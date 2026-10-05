"""Phase error and duty cycle from bertaina2024 / the authors' notebook."""

import numpy as np
from .config import cli_config


def phase_error(variance):
    if np.any(np.asarray(variance) < 0):
        raise ValueError("Variance cannot be negative")
    # Eq. 1, bertaina2024: requested small-phase approximation.
    return np.asarray(variance) / 4


def duty_cycle(tau_q_s, tau_ps_s):
    tau_q_s = np.asarray(tau_q_s, dtype=float)
    if np.any(tau_q_s <= 0) or tau_ps_s < 0:
        raise ValueError("Require tau_Q > 0 and tau_PS >= 0")
    # Sec. I, p. 3 (unnumbered duty formula), bertaina2024; QKD.ipynb duty_stabilization.
    return tau_q_s / (tau_q_s + tau_ps_s)


if __name__ == "__main__":
    from .integration import PhaseIntegral, frequency_grid
    from .spectra import components
    c = cli_config()
    f = frequency_grid(c)
    for s in c["scenarios"]:
        integral = PhaseIntegral(f, components(f, s, c)["total"])
        tau, status = integral.operating_time(c["operation"]["sigma_limit_rad"], c["operation"]["tau_max_s"])
        print(s["name"], "e_phi=", float(phase_error(integral.variance(tau))), "d=", duty_cycle(tau, c["operation"]["tau_ps_s"]), status)
