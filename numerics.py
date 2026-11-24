"""Numerical verification: time-step refinement, eigenvalues and long-term response.

1. The step response of Fig. 6 is computed with steps 0.02, 0.01 and 0.0025 s and the
   maximum differences with respect to the 0.0025 s solution are printed.
2. The Jacobian of the right-hand sides at the trim point (H = 1000 m, V = 190 m/s) is
   computed by central differences; its eigenvalues give the time scales of the model
   and the stability margin of the fourth-order Runge-Kutta method (|lambda|*dt < 2.78).
3. The step response is continued to 1500 s to show the long-term behaviour.

Usage:  python numerics.py      (about 2 minutes)

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from citation import handle_help, print_citation
handle_help(__doc__)
import numpy as np
from common import *


def main():
    x0, u0 = trim2(190., 1000.)

    def ufun(t, x):
        u = u0.copy(); u[3] = u0[3] - 0.3 * d2r
        return u

    sol = {dt: sim(x0, ufun, 120., dt) for dt in (0.02, 0.01, 0.0025)}
    tr, Xr, _ = sol[0.0025]
    print('1. Time-step refinement (0...120 s, reference step 0.0025 s)')
    for dt in (0.02, 0.01):
        t, X, _ = sol[dt]; s = int(round(1 / dt)); idx = np.searchsorted(tr, t[::s]); Xi = X[::s]
        print(f'   dt = {dt:5.3f} s: max |dV| = {np.abs(Xi[:,0]-Xr[idx,0]).max():.1e} m/s, '
              f'max |d alpha| = {np.abs(Xi[:,1]-Xr[idx,1]).max()/d2r:.1e} deg, '
              f'max |dh| = {np.abs(Xi[:,12]-Xr[idx,12]).max():.1e} m')
    f = lambda x: rhs(x, u0)[0]
    n = len(x0); J = np.zeros((n, n)); h = np.maximum(1e-6, 1e-7 * np.abs(x0))
    for i in range(n):
        e = np.zeros(n); e[i] = h[i]; J[:, i] = (f(x0 + e) - f(x0 - e)) / (2 * h[i])
    ev = np.linalg.eigvals(J)
    print('2. Eigenvalues of the linearized model at the trim point')
    for l in sorted(ev, key=lambda z: -abs(z)):
        extra = f'  period {2*np.pi/abs(l.imag):6.1f} s, damping ratio {-l.real/abs(l):.3f}' if abs(l.imag) > 1e-9 else ''
        print(f'   {l.real:11.5f} {l.imag:+11.5f}j{extra}')
    lm = max(abs(ev)); print(f'   largest |lambda| = {lm:.1f} 1/s, |lambda|*dt = {lm*0.01:.2f} for dt = 0.01 s (RK4 limit 2.78)')
    t, X, _ = sim(x0, ufun, 1500., 0.01)
    print('3. Long-term response to the sustained elevator step')
    for T in (0, 100, 200, 400, 800, 1200, 1500):
        k = int(T / 0.01)
        print(f'   t = {T:4d} s: V = {X[k,0]:6.1f} m/s, alpha = {X[k,1]/d2r:5.2f} deg, h = {X[k,12]:6.0f} m, '
              f'flight-path angle = {(X[k,8]-X[k,1])/d2r:5.2f} deg')
    print_citation()


if __name__ == '__main__':
    main()
