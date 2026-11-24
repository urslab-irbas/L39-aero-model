"""Trim values for straight and level flight (Table 2 of the article).

Solves the trim equations of the L-39 model by Newton's method for
H = 1000, 3000, 5000 m and V = 150...190 m/s and saves trim.json with rows
[H (m), V (m/s), M, alpha (deg), throttle (%), elevator (deg), n1 (rpm),
thrust (kN), c_y].

Usage:  python trim_table.py

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from citation import handle_help, print_citation
handle_help(__doc__)

import json
import warnings
from common import *


def main():
    rows = []
    for H in (1000., 3000., 5000.):
        for V in (150., 160., 170., 180., 190.):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                x, u = trim2(V, H)
            _, o = rhs(x, u, clip=False)
            rows.append([H, V, o['M'], x[1] / d2r, u[0] * 100, u[3] / d2r, x[9], o['P'] / 1000, o['cy']])
            print(f'H={H:5.0f} V={V:5.1f}  M={o["M"]:.3f}  alpha={x[1]/d2r:5.2f} deg  '
                  f'throttle={u[0]*100:5.1f} %  elevator={u[3]/d2r:5.2f} deg  n1={x[9]:6.0f}  P={o["P"]/1000:5.2f} kN')
    json.dump(rows, open('trim.json', 'w'))
    print('trim.json saved.')
    print_citation()


if __name__ == '__main__':
    main()
