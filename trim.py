"""Trimming of the L-39 model for level flight (least-squares variant).

Usage:  python trim.py

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
import numpy as np
from scipy.optimize import least_squares
from l39 import *
B0, L0 = 43.45*d2r, 24.53*d2r     # Долна Митрополия
def state(V, al, n1, H, psi=0., ga=0.):
    x = np.array([V, al, 0., 0, 0, 0, psi, ga, al, n1, B0, L0, H])
    return x
def trim(V=150., H=1000., psi=0.):
    def fun(z):
        al, thr, de, n1 = z
        x = state(V, al, n1, H, psi)
        # абсолютни ъглови скорости = преносни (относителните са нула)
        _, o = rhs(x, [thr, 0, 0, de])
        C = Cmat(psi, 0., al); Wn = o['Wn']
        Gm = MA/np.sqrt(1-E2*np.sin(B0)**2); Qm = Gm*(1-E2)/(1-E2*np.sin(B0)**2)
        w = C @ (OM*np.array([0, np.cos(B0), np.sin(B0)]) + np.array([-Wn[1]/(Qm+H), Wn[0]/(Gm+H), Wn[0]*np.tan(B0)/(Gm+H)]))
        x[3:6] = w
        d, _ = rhs(x, [thr, 0, 0, de])
        return [d[0], d[1], d[5]*100, d[9]/1000]
    z = least_squares(fun, [3*d2r, 0.6, -2*d2r, 15500.], bounds=([-5*d2r,0.55,-30*d2r,10000],[15*d2r,1.,30*d2r,17600]), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    al, thr, de, n1 = z.x; z=(z.x,z.cost,z.status)
    x = state(V, al, n1, H, psi)
    return x, np.array([thr, 0., 0., de]), z
if __name__ == '__main__':
    from citation import handle_help, print_citation
    handle_help(__doc__)
    for V in (170., 180., 190., 200.):
        x, u, z = trim(V)
        d, o = rhs(x, u)
        print(f'V={V} al={x[1]/d2r:.3f}deg thr={u[0]:.3f} de={u[3]/d2r:.3f}deg n1={x[9]:.0f} P={o["P"]:.0f}N M={o["M"]:.3f} cy={o["cy"]:.3f} cx={o["cx"]:.4f} ier={z[2]} ', np.round(d[:10],5))
    print_citation()
