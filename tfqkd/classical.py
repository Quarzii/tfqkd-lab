"""Williams A6/A8 with the explicitly user-authorized ideal controller C(s)=g/s."""
import warnings
import numpy as np
from .config import cli_config


def delay_and_boundary(length_km, physics):
    if length_km <= 0:
        raise ValueError('Length must be positive')
    # Appendix A, Eqs. A1, A6, williams2008: tau is ONE-WAY delay.
    tau = physics['n'] * length_km / physics['c_km_s']
    # Appendix A after Eq. A7, williams2008: servo characteristic, not universal hard stability bound.
    return tau, 2*tau, 1/(4*tau)


def remote_power_ratio(f, tau_oneway, controller, spatial_nodes):
    """controller(s) supplies G0*K*F(s)/s; uniform spatially uncorrelated noise."""
    if controller is None:
        raise ValueError('A controller must be supplied explicitly for Eq. A8')
    f=np.asarray(f,dtype=float)
    if np.any(f<=0) or tau_oneway<0:
        raise ValueError('Require positive frequency and nonnegative delay')
    # Eqs. A3, A6, A8, williams2008; evaluate the distributed-noise spatial integral.
    s=2j*np.pi*f
    C=np.asarray(controller(s))
    G=C*(1+np.exp(-2*s*tau_oneway))
    nodes,weights=np.polynomial.legendre.leggauss(spatial_nodes)
    # Eq. A4, williams2008: uniform uncorrelated noise normalized over x=z/L in [0,1].
    x=(nodes+1)/2
    weights=weights/2
    # Eqs. A3 and A8, williams2008; cancel cosh algebraically to avoid artificial singularity.
    amplitude=np.exp(-s[:,None]*tau_oneway*(1-x))-2*(C/(1+G))[:,None]*np.exp(-2*s[:,None]*tau_oneway)*np.cos(2*np.pi*f[:,None]*tau_oneway*(1-x))
    return np.sum(np.abs(amplitude)**2 * weights, axis=1)


def high_gain_ratio(f,tau_oneway,spatial_nodes):
    # Eq. A9, williams2008, uniform spatial profile; formal high-gain limit, not finite servo response.
    nodes,weights=np.polynomial.legendre.leggauss(spatial_nodes)
    x=(nodes+1)/2
    u=2*np.pi*np.asarray(f)*tau_oneway
    return ((np.sin(u[:,None]*x)/np.cos(u[:,None]))**2) @ (weights/2)


def low_frequency_ratio(f,tau_oneway):
    # Eqs. A10-A11 / Eq. 3.1, williams2008: a=1/3 for uniform spatial noise.
    return (2*np.pi*np.asarray(f)*tau_oneway)**2/3




def integrator_controller(g):
    """User-authorized engineering choice C(s)=g/s; not a fitted Williams controller."""
    if not np.isfinite(g) or g <= 0:
        raise ValueError('Require a finite gain g > 0 in s^-1')
    # Eq. A6, williams2008, with user-authorized G0*K*F(s)=g.
    return lambda s: g / s


def loop_gain(f, tau_roundtrip, g):
    f = np.asarray(f, dtype=float)
    if np.any(f <= 0) or tau_roundtrip < 0:
        raise ValueError('Require f>0 and tau_RT>=0')
    # Eq. A6, williams2008, C(s)=g/s is an engineering assumption authorized by user.
    s = 2j * np.pi * f
    return integrator_controller(g)(s) * (1 + np.exp(-s * tau_roundtrip))


def local_sensitivity(f, tau_roundtrip, g):
    # Eq. A7, williams2008: LOCAL sensitivity; not the remote one-way PSD ratio.
    return np.abs(1 / (1 + loop_gain(f, tau_roundtrip, g))) ** 2


def remote_integrator_ratio(f, tau_oneway, g, spatial_nodes=64):
    """Backward-compatible ideal-integrator spatial integral of Williams A8."""
    return remote_controller_ratio(f, tau_oneway, integrator_controller(g))


def remote_controller_ratio(f, tau_oneway, controller):
    """Exact uniform-profile integral of A3/A8 for an explicitly supplied controller."""
    f = np.atleast_1d(np.asarray(f, dtype=float))
    if np.any(f <= 0) or np.any(np.asarray(tau_oneway) < 0):
        raise ValueError('Require f>0 and nonnegative delay')
    # Eqs. A3,A6,A8, williams2008: D multiplies cos(omega*tau*y), y=1-z/L.
    s = 2j * np.pi * f
    C = controller(s)
    G = C * (1 + np.exp(-2 * s * tau_oneway))
    u = 2 * np.pi * f * tau_oneway
    D = 2 * C / (1 + G) * np.exp(-2 * s * tau_oneway)
    u, D = np.broadcast_arrays(u, D)
    # Eqs. A3,A4,A8, williams2008: analytic integrals of cos² and exp(-i*u*y)*cos(u*y).
    cos2 = (1 + np.sinc(2 * u / np.pi)) / 2
    cross = (1 + np.exp(-1j * u) * np.sinc(u / np.pi)) / 2
    result = 1 + np.abs(D)**2 * cos2 - 2 * np.real(np.conj(D) * cross)
    # Same A8 spatial integral; stable low-u evaluation of integral(sin(u*y)^2, y=0..1).
    # Its Taylor coefficients are algebraic consequences of the exact spatial integral,
    # not added physical coefficients. The stopping rule is machine precision.
    small = np.abs(u) < 1
    if np.any(small):
        us = u[small]
        sin2 = np.zeros_like(us)
        term = us**2 / 3
        order = 1
        while np.any(np.abs(term) > np.finfo(float).eps * np.maximum(sin2, np.finfo(float).tiny)):
            sin2 += term
            # Derived integration of sin² in Eqs. A3,A8, williams2008: successive Taylor terms.
            term *= -(2*us)**2 / ((2*order+2)*(2*order+3))
            order += 1
        # Eqs. A3,A8, williams2008: amplitude=(1-D)*cos(u*y)-i*sin(u*y).
        b = 1-D[small]
        result[small] = np.abs(b)**2*(1-sin2)+sin2-np.imag(b)*us*np.sinc(us/np.pi)**2
    if np.any(result < 0) or np.any(~np.isfinite(result)):
        raise ArithmeticError('Invalid remote PSD ratio')
    return result


def actuator_controller(g, omega_a_rad_s):
    """User-authorized engineering assumption C=g/[s*(1+s/omega_a)]."""
    if not np.isfinite(g) or g <= 0 or not np.isfinite(omega_a_rad_s) or omega_a_rad_s <= 0:
        raise ValueError('Require finite positive g [s^-1] and omega_a [rad/s]')
    # Eq. A6, williams2008; additional actuator pole explicitly requested by user, no source coefficient.
    return lambda s: g / (s * (1 + s / omega_a_rad_s))


def actuator_loop_gain(f, tau_roundtrip, g, omega_a_rad_s):
    # Eq. A6, williams2008, with user-authorized one-pole controller.
    s = 2j * np.pi * np.asarray(f, dtype=float)
    if np.any(np.asarray(f) <= 0) or tau_roundtrip < 0:
        raise ValueError('Require positive frequency and nonnegative delay')
    return actuator_controller(g, omega_a_rad_s)(s) * (1 + np.exp(-s * tau_roundtrip))


def actuator_stability(tau_roundtrip, omega_a_rad_s, crossing_index=0):
    """Numerical imaginary-axis crossing of adopted A6; first crossing is g_crit."""
    from scipy.optimize import brentq
    if tau_roundtrip <= 0 or omega_a_rad_s <= 0 or crossing_index < 0:
        raise ValueError('Require positive delay/pole and nonnegative crossing index')
    # Eq. A6, williams2008 + user-requested actuator pole: b=omega_a*tau_RT, x=omega*tau_RT.
    b = omega_a_rad_s * tau_roundtrip
    # Algebraic imaginary-axis conditions from A6: x*tan(x/2)=b,
    # with sin(x)>0, hence x lies in (2*k*pi,(2*k+1)*pi).
    left = 2 * crossing_index * np.pi
    right = (2 * crossing_index + 1) * np.pi
    x = brentq(lambda value: value * np.tan(value / 2) - b,
               np.nextafter(left, right), np.nextafter(right, left))
    # Eq. A6 characteristic equation: g_crit=omega/sin(x)=(x^2+b^2)/(2*b*tau_RT).
    # This equality proves higher crossing indices require higher gains.
    omega = x / tau_roundtrip
    gain = (x * x + b * b) / (2 * b * tau_roundtrip)
    loop = actuator_loop_gain([omega / (2 * np.pi)], tau_roundtrip, gain, omega_a_rad_s)[0]
    return dict(g_crit_per_s=float(gain), critical_frequency_hz=float(omega / (2 * np.pi)),
                x=float(x), b=float(b), unity_gain_error=float(abs(abs(loop) - 1)),
                negative_real_axis_error=float(abs(loop + 1)))


def require_stable_actuator(tau_roundtrip, g, omega_a_rad_s):
    result = actuator_stability(tau_roundtrip, omega_a_rad_s)
    if not np.isfinite(g) or g <= 0 or g >= result['g_crit_per_s']:
        raise ValueError(f'Unstable/marginal controller: g={g} >= g_crit={result["g_crit_per_s"]} s^-1')
    return result


def actuator_local_sensitivity(f, tau_roundtrip, g, omega_a_rad_s):
    require_stable_actuator(tau_roundtrip, g, omega_a_rad_s)
    # Eq. A7, williams2008, with the user-authorized actuator pole in A6.
    return np.abs(1 / (1 + actuator_loop_gain(f, tau_roundtrip, g, omega_a_rad_s))) ** 2


def actuator_characteristic(z, a, b):
    # Eq. A6 + user-authorized pole: z=s*tau_RT, a=g*tau_RT, b=omega_a*tau_RT.
    return z * (1 + z / b) + a * (1 + np.exp(-z))


def actuator_pole(a, b, initial):
    """Solve an adopted A6 characteristic root; used only for quadrature diagnostics."""
    from scipy.optimize import root
    def fun(vector):
        with np.errstate(over='ignore', invalid='ignore'):
            value = actuator_characteristic(complex(*vector), a, b)
        return [value.real, value.imag]
    def jac(vector):
        # Derivative of the adopted Eq. A6 characteristic equation (numerical root finding).
        with np.errstate(over='ignore', invalid='ignore'):
            derivative = 1 + 2 * complex(*vector) / b - a * np.exp(-complex(*vector))
        return [[derivative.real, -derivative.imag], [derivative.imag, derivative.real]]
    solution = root(fun, [initial.real, initial.imag], jac=jac, tol=1e-11)
    z = complex(*solution.x)
    # Eq. A6: scale the characteristic residual by the magnitudes of its two terms.
    # Eq. A6 root verification: overflowing trial seeds are rejected, not used in the spectrum.
    with np.errstate(over='ignore', invalid='ignore'):
        scale = 1 + abs(z * (1 + z / b)) + abs(a * (1 + np.exp(-z)))
        residual = abs(actuator_characteristic(z, a, b)) / scale
    if not np.isfinite(residual) or residual > 1e-9:
        raise ArithmeticError(f'Characteristic root did not converge: a={a}, b={b}, residual={residual}')
    return z


def actuator_resonance_poles(tau_roundtrip, g, omega_a_rad_s, branches):
    """Find positive-frequency delay roots to resolve narrow A8 servo peaks."""
    stability = require_stable_actuator(tau_roundtrip, g, omega_a_rad_s)
    # Eq. A6 with adopted pole: nondimensional characteristic equation.
    a, b = g * tau_roundtrip, omega_a_rad_s * tau_roundtrip
    seeds = [1j * stability['x']]
    # Low-delay expansion of the same A6 characteristic equation; numerical seeds only.
    seeds.extend(complex(z) for z in np.roots([1 / b, 1 - a, 2 * a]))
    # Newton seeds near zeros of 1+exp(-z), derived from the adopted A6 equation.
    for index in range(branches):
        q = (2 * index + 1) * np.pi
        at_zero = 1j * q
        seeds.append(at_zero - actuator_characteristic(at_zero, a, b) / (1 + a + 2 * at_zero / b))
    if a <= 600:
        seeds.extend(complex(z) for z in dimensionless_poles(a, branches) if z.imag > 0)
    poles = []
    for seed in seeds:
        try:
            z = actuator_pole(a, b, seed)
        except ArithmeticError:
            continue
        if z.imag > 0 and z.real < 0 and not any(abs(z - previous) <= 1e-7 * (1 + abs(z)) for previous in poles):
            poles.append(z)
    return np.asarray(sorted(poles, key=lambda z: z.imag), dtype=complex)


def dimensionless_poles(g_tau, branches):
    """Poles z=s*tau_RT of z+a+a*exp(-z)=0, derived from adopted Eq. A6."""
    from scipy.special import lambertw
    if g_tau <= 0 or g_tau > 600:
        raise ValueError('Pole diagnostic requires 0<g*tau_RT<=600 to avoid exponential overflow')
    # Eq. A6, williams2008, adopted C=g/s: algebraic solution z=W_k(-a*exp(a))-a.
    return lambertw(-g_tau * np.exp(g_tau), np.arange(-branches, branches + 1)) - g_tau


def phase_resonance_frequency(tau_roundtrip):
    """Numerical first delayed-path cancellation; infinite-g pole limit, not stability cutoff."""
    from scipy.optimize import brentq
    if tau_roundtrip <= 0:
        raise ValueError('Positive round-trip delay required')
    # Eq. A6 / Appendix A after A7, williams2008: solve omega*tau_RT=pi numerically.
    omega = brentq(lambda value: value * tau_roundtrip - np.pi, 0, 2 * np.pi / tau_roundtrip)
    return omega / (2 * np.pi)


if __name__ == '__main__':
    c = cli_config()
    tau, rt, _ = delay_and_boundary(c['classical']['length_km'], c['physics'])
    print('C(s)=g/s is a user-authorized engineering assumption; g_crit=+infinity for this ideal loop.')
    print('One-way delay, roundtrip delay [s]:', tau, rt)
    print('Numerical limiting first resonance [Hz]:', phase_resonance_frequency(rt))
    for g in c['classical']['comparison_g_per_s']:
        print('g [s^-1]', g, 'remote ratio at fc1:', remote_integrator_ratio([c['physics']['fc1_hz']], tau, g)[0])
