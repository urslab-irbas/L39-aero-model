"""Comparison of the model with independent L-39 flight-manual data (Table 3 of the article).

Steady-state performance computed from the aerodynamic and engine models
(lift = weight, thrust = drag, rotor speed in steady state, throttle at maximum):
static thrust, maximum level-flight speed, maximum rate of climb, stall speed
(clean configuration) and roll rate at full aileron deflection.
Reference values: Flight manual L-39 (MO SSSR, 1988, pp. 216-217) and
T.O. 1T-L39C-1 (Aero Vodochody), Fig. 6-1 and the roll-rate range 100-150 deg/s.

Usage:  python validation.py

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from citation import handle_help, print_citation
handle_help(__doc__)
import numpy as np
from scipy.optimize import brentq
from l39 import atm, coeffs, engine, S, bref, d2r

G0 = 9.80665


def n1_ss(V, H, thr):
    """Steady-state high-pressure rotor speed (dn1/dt = 0)."""
    return brentq(lambda n: engine(n, V, H, thr)[0], 5000., 25000.)


def thrust_drag(V, H, thr, mass):
    ro, a = atm(H); q = ro * V * V / 2
    cl = mass * G0 / (q * S)                       # lift coefficient for level flight
    cd = coeffs(V / a)['cx1'] + 0.09 * cl ** 2     # drag polar, eq. (18)
    return engine(n1_ss(V, H, thr), V, H, thr)[1], q * S * cd


def v_max(H, mass=4000., thr=1.):
    f = lambda V: np.subtract(*thrust_drag(V, H, thr, mass))
    Vs = np.arange(120., 300., 1.); ex = [f(V) for V in Vs]
    i = max(k for k in range(len(Vs) - 1) if ex[k] > 0 >= ex[k + 1])
    return brentq(f, Vs[i], Vs[i + 1])


def roc_max(H, mass=4300., thr=1.):
    best = max(((lambda P, D: (P - D) * V / (mass * G0))(*thrust_drag(V, H, thr, mass)), V)
               for V in np.arange(70., 260., 0.5))
    return best


def stall_speed(mass, M=0.15):
    F = coeffs(M); cl_max = -0.007 + F['cy1'] * (F['cy6'] - 0.012217)    # eq. (17) at the critical angle
    return np.sqrt(2 * mass * G0 / (1.225 * S * cl_max)), cl_max


def roll_rate(V=380 / 3.6, H=1000., da=17.):
    """Mean roll rate during the first full roll (360 deg of bank) after a full aileron
    deflection from the trimmed flight at V = 380 km/h (simulation of the full model)."""
    import warnings
    from common import trim2, sim
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        x0, u0 = trim2(V, H)

    def ufun(t, x):
        u = u0.copy(); u[1] = da * d2r
        return u

    t, X, U = sim(x0, ufun, 5.0, 0.01)
    k = np.argmax(np.abs(X[:, 7] - X[0, 7]) >= 2 * np.pi)
    return 360. / t[k], t[k]


def main():
    print('Quantity                                   model    flight manual')
    P0 = engine(n1_ss(0.1, 0., 1.), 0.1, 0., 1.)[1]
    print(f'static thrust, kN                          {P0/1000:6.1f}    16.9 (T.O. 1T-L39C-1)')
    for H, ref in ((0., 702.), (5000., 757.), (10000., 737.)):
        print(f'max level speed, H = {H:5.0f} m, km/h (TAS)  {v_max(H)*3.6:6.0f}    {ref:.0f}')
    for H, ref in ((0., 22.), (6000., 10.8), (10000., 3.4)):
        r, V = roc_max(H)
        print(f'max rate of climb, H = {H:5.0f} m, m/s       {r:6.1f}    {ref}   (at {V*3.6:.0f} km/h TAS)')
    for mass, ref in ((3700., 172.), (4500., 190.)):
        Vs, clm = stall_speed(mass)
        print(f'stall speed, {mass:.0f} kg, km/h               {Vs*3.6:6.0f}    {ref:.0f}   (C_L,max = {clm:.2f})')
    p, T = roll_rate()
    print(f'roll rate, full aileron, 380 km/h, deg/s    {p:6.0f}    100-150   (first 360 deg in {T:.2f} s)')
    print_citation()


if __name__ == '__main__':
    main()
