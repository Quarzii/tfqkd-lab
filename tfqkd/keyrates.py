"""Stage 2 key-rate models from bertaina2024 and QKD.ipynb."""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math

import numpy as np

from .config import cli_config, load
from .integration import PhaseIntegral, frequency_grid
from .protocol import duty_cycle
from .spectra import components


@dataclass(frozen=True)
class KeyRateParameters:
    clockrate_hz: float
    attenuation_db_per_km: float
    stab_overhead_s: float
    detector_efficiency: float
    detector_dark_count_rate_hz: float
    detector_error: float
    decoy_big: float
    decoy_medium: float
    decoy_mini: float
    f_error: float
    pz_sns: float
    eps_sns_aopp: float
    duty_bb84: float
    pz_cal: float
    nmin_cal: int
    nmax_cal: int
    u_cal: float

    def asdict(self):
        return asdict(self)


def table_ii_sns_pd_parameters():
    # Table II, bertaina2024; QKD.ipynb Cell 23, parameters_source/channels/protocols/SNSPD.
    return KeyRateParameters(**load()["keyrate"])


def h2(probability):
    p = np.asarray(probability, dtype=float)
    out = np.zeros_like(p, dtype=float)
    mask = p >= 1.0e-15
    # Binary Shannon entropy H2, Appendix A Eq. A1 / Appendix B Eq. B2 / Appendix C Eq. C7, bertaina2024; QKD.ipynb Cell 13.
    out[mask] = -p[mask] * np.log2(p[mask]) - (1 - p[mask]) * np.log2(1 - p[mask])
    return float(out) if out.ndim == 0 else out


def etadb(loss_db):
    # Channel transmittance from dB loss, QKD.ipynb Cell 13 etadB.
    return 10 ** (-np.asarray(loss_db, dtype=float) / 10)


def phase_error_exact(sigma_phi):
    # Gaussian phase-error model, main text Eq. 1 discussion and QKD.ipynb Cell 13 EphaseFromSigmaPhi.
    return (1 - np.exp(-(np.asarray(sigma_phi, dtype=float) ** 2) / 2)) / 2


def realistic_secret_key_capacity(loss_db, detector_efficiency):
    # Realistic PLOB bound with detector efficiency, Sec. II and QKD.ipynb Cell 15.
    rate = -np.log2(1 - detector_efficiency * etadb(loss_db))
    return np.where(rate > 1, np.nan, rate)


def decoy_bb84_rate_per_pulse(loss_db, sigma_phi, params: KeyRateParameters):
    # QKD.ipynb Cell 17 RdB_Decoy_err; Appendix C, bertaina2024.
    pdc = params.detector_dark_count_rate_hz / params.clockrate_hz
    # QKD.ipynb Cell 17 RdB_Decoy_err; Appendix C, bertaina2024.
    eta = params.detector_efficiency * etadb(loss_db)
    u, v, w = params.decoy_big, params.decoy_medium, params.decoy_mini
    # Appendix C Eq. C1, bertaina2024; QKD.ipynb Cell 17 GainDecoy_model.
    qu = 1 - (1 - pdc) * np.exp(-u * eta)
    # Appendix C Eq. C1, bertaina2024; QKD.ipynb Cell 17 GainDecoy_model.
    qv = 1 - (1 - pdc) * np.exp(-v * eta)
    # Appendix C Eq. C1, bertaina2024; QKD.ipynb Cell 17 GainDecoy_model.
    qw = 1 - (1 - pdc) * np.exp(-w * eta)
    ephase = phase_error_exact(sigma_phi)
    # QKD.ipynb Cell 17 RdB_Decoy_err; Appendix C, bertaina2024.
    edet = params.detector_error + ephase
    # Appendix C Eq. C2, bertaina2024; QKD.ipynb Cell 17 ErrorDecoy_model.
    eu = (pdc / 2 + (edet - pdc / 2) * (1 - np.exp(-u * eta))) / qu
    # Appendix C Eq. C2, bertaina2024; QKD.ipynb Cell 17 ErrorDecoy_model.
    ev = (pdc / 2 + (edet - pdc / 2) * (1 - np.exp(-v * eta))) / qv
    # Appendix C Eq. C2, bertaina2024; QKD.ipynb Cell 17 ErrorDecoy_model.
    ew = (pdc / 2 + (edet - pdc / 2) * (1 - np.exp(-w * eta))) / qw
    # Appendix C Eq. C3, bertaina2024.
    y0_lower = (v * qw * np.exp(w) - w * qv * np.exp(v)) / (v - w)
    # Appendix C Eq. C4, bertaina2024.
    y1_lower = (u**2 * (qv * np.exp(v) - qw * np.exp(w)) - (v**2 - w**2) * (qu * np.exp(u) - y0_lower)) / (u * (u - v - w) * (v - w))
    # Appendix C Eq. C5, bertaina2024.
    q1_lower = y1_lower * u / np.exp(u)
    # Appendix C Eq. C6, bertaina2024.
    e1_upper = (ev * qv * np.exp(v) - ew * qw * np.exp(w)) / ((v - w) * y1_lower)
    # Appendix C Eq. C7, bertaina2024.
    return params.duty_bb84 * (q1_lower * (1 - h2(e1_upper)) - qu * params.f_error * h2(eu))


def cal_rate_per_pulse(loss_db, sigma_phi, params: KeyRateParameters):
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    pdc = params.detector_dark_count_rate_hz / params.clockrate_hz
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    sqrt_eta = params.detector_efficiency * etadb(loss_db / 2)
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    theta = 2 * np.arcsin(np.sqrt(params.detector_error))
    alpha2 = params.u_cal
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    gamma = sqrt_eta * alpha2
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    omega = np.cos(theta) * np.cos(sigma_phi)
    # Appendix B Eq. B3, bertaina2024; QKD.ipynb Cell 19 pZZ.
    pzz = 0.5 * (1 - pdc) * (np.exp(-gamma * omega) + np.exp(gamma * omega)) * np.exp(-gamma) - (1 - pdc) ** 2 * np.exp(-2 * gamma)
    # Appendix B Eq. B4, bertaina2024; QKD.ipynb Cell 19 eZ.
    ez = (np.exp(-gamma * omega) - (1 - pdc) * np.exp(-gamma)) / (np.exp(-gamma * omega) + np.exp(gamma * omega) - 2 * (1 - pdc) * np.exp(-gamma))
    ez = np.minimum(0.5, ez)

    def summ(x, j, n_min, n_max):
        # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 summ finite S_j truncation.
        return sum(x ** (2 * n + j) / math.sqrt(math.factorial(2 * n + j)) for n in range(n_min, n_max))

    def stot(x, j, n_min, n_max):
        # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 Stot.
        return np.exp(-(x**2)) * summ(x, j, n_min, n_max) ** 2

    def delta_j(x, j, n_min, n_max, xj):
        # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 Deltaj.
        return stot(x, j, n_min, n_max) - xj

    def pxx(qxx, n_a, n_b):
        # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 pXX photon-number gain.
        return (1 - pdc) * (pdc * (1 - sqrt_eta) ** (n_a + n_b) + qxx)

    # Appendix B Eq. B5 with S0={(0,0),(0,1),(1,0),(1,1)}, bertaina2024; QKD.ipynb Cell 19 qXX02.
    qxx02 = sqrt_eta * (1 - sqrt_eta) + sqrt_eta**2 / 4
    # Appendix B Eq. B5 with S1={(0,0)}, bertaina2024; QKD.ipynb Cell 19 qXX11.
    qxx11 = sqrt_eta * (1 - sqrt_eta) / 2 + sqrt_eta * (1 - sqrt_eta) / 2 + sqrt_eta**2 * np.cos(theta) ** 2 / 4
    # Appendix B Eq. B5 with S0={(0,0),(0,1),(1,0),(1,1)}, bertaina2024; QKD.ipynb Cell 19 qXX22.
    qxx22 = 2 * sqrt_eta * (1 - sqrt_eta) ** 3 + (sqrt_eta**2 * (1 - sqrt_eta) ** 2 / 4) * (2 + 4 * (1 + np.cos(theta) ** 2)) + (sqrt_eta**3 * (1 - sqrt_eta) / 8) * 60 + (sqrt_eta**4 / 64) * 384
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    x0 = np.exp(-alpha2) * (1 + np.sqrt(2) * alpha2 + alpha2**2 / 2)
    # QKD.ipynb Cell 19 RdB_CAL_TF_QKD_err; Appendix B, bertaina2024.
    x1 = np.exp(-alpha2) * (alpha2 + 2 / np.sqrt(6) * alpha2**2)
    delta0 = delta_j(np.sqrt(alpha2), 0, params.nmin_cal, params.nmax_cal, x0)
    # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 temp00.
    temp00 = np.exp(-alpha2 / 2) ** 2 * np.sqrt(pxx(0, 0, 0)) + np.sqrt(2) * alpha2 * np.sqrt(pxx(qxx02, 0, 2)) + (alpha2**2 / 2) * np.sqrt(pxx(qxx22, 2, 2)) + delta0
    delta1 = delta_j(np.sqrt(alpha2), 1, params.nmin_cal, params.nmax_cal, x1)
    # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 temp11.
    temp11 = (np.exp(-alpha2 / 2) * np.sqrt(alpha2)) ** 2 * np.sqrt(pxx(qxx11, 1, 1)) + np.exp(-alpha2) * 2 / np.sqrt(6) * alpha2**2 + delta1
    # Appendix B Eq. B5, bertaina2024; QKD.ipynb Cell 19 eX.
    ex = (temp00**2 + temp11**2) / pzz
    # Appendix B Eq. B2 and Eq. B1, bertaina2024; QKD.ipynb Cell 19 Rlow doubles symmetric single-click events.
    return params.pz_cal**2 * 2 * pzz * (1 - params.f_error * h2(ez) - h2(np.minimum(0.5, ex)))


def sns_aopp_rate_per_pulse(loss_db, sigma_phi, params: KeyRateParameters):
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    pdc = params.detector_dark_count_rate_hz / params.clockrate_hz
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    eta_tf = params.detector_efficiency * etadb(loss_db / 2)
    eps = params.eps_sns_aopp
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    s = params.decoy_big / 2
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    n = params.decoy_mini / 2
    decoy2, decoy1, decoy0 = params.decoy_big, params.decoy_medium, params.decoy_mini
    ephase = phase_error_exact(sigma_phi)
    edet = params.detector_error

    def gain(mu):
        # Appendix C Eq. C1 adapted to two TF-QKD detectors, bertaina2024; QKD.ipynb Cell 21 Gain_SNS_TFQKD_model.
        return 1 - (1 - pdc) ** 2 * np.exp(-mu * eta_tf)

    def error(mu, q):
        # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
        pdcma = 2 * pdc
        # Appendix C Eq. C2 adapted to two TF-QKD detectors, bertaina2024; QKD.ipynb Cell 21 Error_SNS_TFQKD_model.
        return (pdcma / 2 + (edet + ephase - pdcma / 2) * (1 - np.exp(-mu * eta_tf))) / q

    q2, q1, q0 = gain(decoy2), gain(decoy1), gain(decoy0)
    e1, e0 = error(decoy1, q1), error(decoy0, q0)
    # Appendix C Eq. C3, bertaina2024; QKD.ipynb Cell 21.
    y0_lower = (decoy1 * q0 * np.exp(decoy0) - decoy0 * q1 * np.exp(decoy1)) / (decoy1 - decoy0)
    # Appendix C Eq. C4, bertaina2024; QKD.ipynb Cell 21.
    y1_lower = (decoy2**2 * (q1 * np.exp(decoy1) - q0 * np.exp(decoy0)) - (decoy1**2 - decoy0**2) * (q2 * np.exp(decoy2) - y0_lower)) / (decoy2 * (decoy2 - decoy1 - decoy0) * (decoy1 - decoy0))
    # Appendix C Eq. C6, bertaina2024; QKD.ipynb Cell 21.
    e1_upper = (e1 * q1 * np.exp(decoy1) - e0 * q0 * np.exp(decoy0)) / ((decoy1 - decoy0) * y1_lower)
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    pz_tilde = eps * (1 - eps)
    # Appendix A Eq. A1 untagged one-photon term after decoy estimate; QKD.ipynb Cell 21.
    n11 = pz_tilde * (s * np.exp(-s) * np.exp(-n) + np.exp(-s) * n * np.exp(-n)) * y1_lower
    n10 = n11
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    qu = gain(s + s)
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    qua = gain(s + n)
    qub = qua
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    qvac = gain(n + n)
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    nc1 = pz_tilde * qua
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    nc0 = pz_tilde * qub
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    nd = eps**2 * qu
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    nv = (1 - eps) ** 2 * qvac
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    n0_total = nd + nc0
    # QKD.ipynb Cell 21 RdB_SNS_AOPP_TF_QKD_err; Appendix A, bertaina2024.
    n1_total = nv + nc1
    # SNS-AOPP odd-pair count min(N0,N1), Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    naopp = np.minimum(n0_total, n1_total)
    # SNS-AOPP retained error-pair contribution, Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    nvddv = (nd / n0_total) * (nv / n1_total) * naopp
    # SNS-AOPP retained correct-pair contribution, Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    nc1c0 = (nc0 / n0_total) * (nc1 / n1_total) * naopp
    # SNS-AOPP retained string length n_t', Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    n_tilde = nc1c0 + nvddv
    # SNS-AOPP bit-flip error E_Z', Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    ez_prime = nvddv / n_tilde
    # SNS-AOPP untagged bits n_1', Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    n1_prime = (n10 / n0_total) * (n11 / n1_total) * naopp
    # SNS-AOPP phase error e_1'^ph, Appendix A text after Eq. A1; QKD.ipynb Cell 21.
    e1_prime = 2 * e1_upper * (1 - e1_upper)
    # Appendix A Eq. A1 with SNS-AOPP primed quantities, bertaina2024; QKD.ipynb Cell 21.
    return params.pz_sns**2 * (n1_prime * (1 - h2(e1_prime)) - params.f_error * n_tilde * h2(ez_prime))


def loss_configuration(config, scenario, loss_db, params: KeyRateParameters, points=None):
    # QKD.ipynb calc_sigma_tau_loss, Cell 3; Eq. 4, bertaina2024.
    total_km = loss_db / params.attenuation_db_per_km
    # QKD.ipynb calc_sigma_tau_loss, Cell 3; Eq. 4, bertaina2024.
    la_km = total_km / 2
    # QKD.ipynb calc_sigma_tau_loss, Cell 3; Eq. 4, bertaina2024.
    lb_km = la_km - scenario["delta_L_km"]
    if lb_km < 0:
        raise ValueError("LB cannot be negative: increase loss or reduce deltaL")
    f = frequency_grid(config, points=points)
    integral = PhaseIntegral(f, components(f, scenario, config, LB_km=lb_km)["total"])
    tau_q, status = integral.operating_time(config["operation"]["sigma_limit_rad"], config["operation"]["tau_max_s"])
    # Eq. 4, bertaina2024: sigma_phi is the square root of integrated phase variance.
    sigma_phi = float(np.sqrt(integral.variance(tau_q)))
    return {"tau_q_s": tau_q, "status": status, "sigma_phi_rad": sigma_phi, "loss_db": float(loss_db), "lb_km": float(lb_km)}


def rates_for_loss(config, scenario, loss_db, params=None, points=None):
    # QKD.ipynb calc_sigma_tau_loss, Cell 3; Eq. 4, bertaina2024.
    p = KeyRateParameters(**config["keyrate"]) if params is None else params
    cfg = loss_configuration(config, scenario, loss_db, p, points=points)
    tau_q = cfg["tau_q_s"]
    sigma_phi = cfg["sigma_phi_rad"]
    d = duty_cycle(tau_q, p.stab_overhead_s)
    return {
        **cfg,
        # Sec. II realistic PLOB bound, bertaina2024; QKD.ipynb Cell 23 plot_panel_scenarios.
        "plob_bps": float(p.clockrate_hz * realistic_secret_key_capacity(loss_db, p.detector_efficiency)),
        # Appendix C Eq. C7, bertaina2024; QKD.ipynb Cell 23 uses SigmaPhi=0 for BB84.
        "bb84_bps": float(p.clockrate_hz * decoy_bb84_rate_per_pulse(loss_db, 0, p)),
        # Appendix A Eq. A1 and SNS-AOPP text, bertaina2024; QKD.ipynb Cell 23 duty multiplication.
        "sns_aopp_bps": float(p.clockrate_hz * d * sns_aopp_rate_per_pulse(loss_db, sigma_phi, p)),
        # Appendix B Eq. B1-B2, bertaina2024; QKD.ipynb Cell 23 duty multiplication.
        "cal_bps": float(p.clockrate_hz * d * cal_rate_per_pulse(loss_db, sigma_phi, p)),
        "duty": float(d),
    }


def figure3_panel_scenarios(config):
    by_name = {scenario["name"]: scenario for scenario in config["scenarios"]}
    return {"a": by_name["4"], "b": by_name["5"], "c": by_name["3"], "d": by_name["6"]}


def config_toml_snippet(params=None):
    p = table_ii_sns_pd_parameters() if params is None else params
    rows = ["[keyrate]"]
    labels = {
        "clockrate_hz": "Hz; Table II, bertaina2024",
        "attenuation_db_per_km": "dB/km; Table II, bertaina2024",
        "stab_overhead_s": "s; Table II, bertaina2024",
        "detector_efficiency": "SNSPD eta_D; Table II, bertaina2024",
        "detector_dark_count_rate_hz": "Hz; SNSPD P_DC, Table II, bertaina2024",
        "detector_error": "e_theta; Table II, bertaina2024",
        "decoy_big": "u; Table II, bertaina2024",
        "decoy_medium": "v; Table II, bertaina2024",
        "decoy_mini": "w; Table II, bertaina2024",
        "f_error": "f_EC; Table II, bertaina2024",
        "pz_sns": "QKD.ipynb Cell 23, asymptotic pZ_SNS",
        "eps_sns_aopp": "epsilon; Table II, bertaina2024",
        "duty_bb84": "QKD.ipynb Cell 23, asymptotic efficient BB84",
        "pz_cal": "QKD.ipynb Cell 23, asymptotic pZ_CAL",
        "nmin_cal": "QKD.ipynb Cell 23, CAL B5 truncation",
        "nmax_cal": "QKD.ipynb Cell 23, CAL B5 truncation",
        "u_cal": "mu_zeta; Table II, bertaina2024",
    }
    for key, value in p.asdict().items():
        rows.append(f"{key} = {value!r} # {labels[key]}")
    return "\n".join(rows)


if __name__ == "__main__":
    config = cli_config()
    params = table_ii_sns_pd_parameters()
    print(config_toml_snippet(params))
    for panel, scenario in figure3_panel_scenarios(config).items():
        row = rates_for_loss(config, scenario, loss_db=40.0, params=params)
        print(panel, scenario["name"], row)
