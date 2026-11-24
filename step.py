"""Step deflection of the elevator (Fig. 6 of the article).

Trims the L-39 for straight and level flight at H = 1000 m and V = 190 m/s,
deflects the elevator by -0.3 deg at t = 5 s with a fixed throttle and
integrates the equations of motion for 120 s (Runge-Kutta 4, step 0.01 s).
The result is saved in step.npz (t, X, U, ny, th). The response is then
continued to 800 s from the moment of the elevator step (a) with the throttle
fixed and (b) with the throttle law for reaching the altitude H1 = 1500 m,
    throttle = thr1 + k_h (H1 - h) - k_hd dh/dt,  k_h = 1e-4 1/m, k_hd = 1e-3 s/m,
where thr1 is the throttle that trims the new elevator deflection at H1;
the results are saved in long.npz. figs.py and figs_en.py use both files for Fig. 6.

Usage:  python step.py

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from citation import handle_help, print_citation
handle_help(__doc__)

import numpy as np
from common import *

T_STEP, D_ELEV, T_END = 5.0, -0.3, 120.0
H1, K_H, K_HD, T_LONG = 1500., 1e-4, 1e-3, 800.


def throttle_for(H, de, V_range=(160., 215.)):
    """Throttle and speed of level flight at altitude H trimmed with the elevator deflection de."""
    from scipy.optimize import brentq
    f = lambda V: trim2(V, H)[1][3] - de
    Vs = np.arange(V_range[0], V_range[1] + 1, 5.); fv = [f(V) for V in Vs]
    k = [i for i in range(len(Vs) - 1) if fv[i] * fv[i + 1] < 0][0]
    V1 = brentq(f, Vs[k], Vs[k + 1])
    return trim2(V1, H)[1][0], V1


def main():
    x0, u0 = trim2(190., 1000.)

    def ufun(t, x):
        u = u0.copy()
        if t >= T_STEP:
            u[3] = u0[3] + D_ELEV * d2r
        return u

    t, X, U = sim(x0, ufun, T_END)
    ny = np.empty(len(t)); th = np.empty(len(t))
    for k, (x, u) in enumerate(zip(X, U)):
        _, o = rhs(x, u)
        ro, a = atm(x[12]); q = ro * x[0] ** 2 / 2
        Rb = Qmat(x[2], x[1]) @ (q * S * np.array([-o['cx'], o['cy'], 0.]))
        ny[k] = (Rb[1] / m) / 9.80665                      # normal load factor n_y1
        th[k] = np.arcsin(o['Wn'][2] / np.linalg.norm(o['Wn']))  # flight-path angle
    np.savez('step.npz', t=t, X=X, U=U, ny=ny, th=th)
    print(f'step.npz saved: {len(t)} points, max n_y1 = {ny.max():.3f}')
    import warnings
    de1 = u0[3] + D_ELEV * d2r
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        thr1, V1 = throttle_for(H1, de1)

    def u_fixed(t, x):
        u = u0.copy(); u[3] = de1
        return u

    def u_law(t, x):
        u = u0.copy(); u[3] = de1
        hdot = x[0] * np.sin(x[8] - x[1])          # wings level, no sideslip
        u[0] = min(max(thr1 + K_H * (H1 - x[12]) - K_HD * hdot, 0.55), 1.0)
        return u

    tl, Xf, Uf = sim(x0, u_fixed, T_LONG)
    _, Xl, Ul = sim(x0, u_law, T_LONG)
    d = slice(None, None, 10)                     # every 0.1 s is enough for the figure
    np.savez_compressed('long.npz', t=tl[d], Xf=Xf[d], Uf=Uf[d], Xl=Xl[d], Ul=Ul[d], thr1=thr1, H1=H1)
    hl = Xl[:, 12]; ks = np.where(np.abs(hl - H1) > 15)[0]
    print(f'long.npz saved: fixed throttle h({T_LONG:.0f} s) = {Xf[-1,12]:.0f} m; throttle law thr1 = {thr1*100:.1f} %, '
          f'h({T_LONG:.0f} s) = {hl[-1]:.0f} m, max h = {hl.max():.0f} m, |h - H1| < 15 m after {tl[ks[-1]]:.0f} s, '
          f'throttle {Ul[:,0].min()*100:.0f}...{Ul[:,0].max()*100:.0f} %')
    print_citation()


if __name__ == '__main__':
    main()
