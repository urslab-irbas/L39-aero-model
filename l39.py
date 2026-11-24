"""L-39 aerodynamic model with polynomial approximation of the aerodynamic
characteristics, engine model and navigation equations on the WGS-84 ellipsoid.

Aerodynamics and engine after the manual of N. Patev (2004) for the TL-39
flight simulator; according to Patev (2004, p. 14) the aerodynamic relations
were derived from the aerodynamic characteristics of the L-39 design
specification and corrected by wind-tunnel and flight-test results (polynomials in the Mach number and in the lift coefficient);
the aerodynamic force is rotated with the matrix Q(beta, alpha); the diurnal
rotation of the Earth is taken into account; navigation equations
dB/dt = W_N/(Q+h), dL/dt = W_E/((G+h)cos B), dh/dt = W_U.

Main functions: rhs(x, u) - right-hand sides of the 13 state equations,
rk4(...) - one Runge-Kutta step, coeffs(M) - aerodynamic polynomials,
engine(...) - thrust model, atm(H) - standard atmosphere.

State x = [V, alpha, beta, wx, wy, wz, psi, gamma, theta, n1, B, L, h]
Control u = [throttle, aileron, rudder, elevator] (rad)

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
Run  python citation.py  for the BibTeX entry.
"""
import numpy as np
d2r = np.pi/180
MA, E2 = 6378137.0, 6.69437999014e-03
OM = 7.292115e-5
Ix, Iy, Iz = 535., 2310., 1890.
m = 4100.; S = 18.8; bref = 9.11; cref = 2.15

def atm(H):
    ro = 1.225*(1-H/44300)**4.256
    a = 20.047*np.sqrt(288.15-6.5e-3*H)          # ISA (в AeroL39.m: 300 K)
    return ro, a

def coeffs(M, cy=None):
    P = np.polyval
    F = dict(
        cx1=P([4.6071,-8.7777,6.1683,-1.942,0.2661,-0.0122,0.023], M),
        cy1=P([427.78,-1112.2,1037.3,-427.7,77.007,-4.6132,4.4158], M),
        cy5=P([37.442,-74.654,44.817,-5.7505,-1.0525,0.1506,0.914], M),
        cy6=P([-0.0388,-3.0388,4.059,-1.252,0.0462,0.3], M),
        my3=P([-0.058718,0.10277,0.033169,0.19524], M),
        mz1=P([-0.36451,0.36082,-0.1099,0.010243,0.042108], M),
        mz5=P([268.56,-858.15,1093.7,-708.80,247.65,-45.782,4.0469,-0.12810,0.32999], M),
        mz11=P([911.16,-2187.4,1861.9,-676.29,104.45,-5.2877,19.039], M),
        mz12=P([2331.7,-5618.5,4904.8,-1904.4,323.05,-19.044,6.1631], M),
        mx1=P([0.504,-1.133,0.7879,0,-0.2676,0.1421,-0.0306,0.0023,0.0015], M),
        mx4=P([0.0015627,-0.0038264,0.0033264,0.0011537], M))
    return F

def cy_fun(M, a):
    F = coeffs(M); fak = 1. if a > F['cy6'] else 0.
    return -0.007 + F['cy1']*(a-0.012217)*(1-fak) + F['cy5']*fak, fak, F

def engine(n1, V, H, thr):
    maxn1, maxP, maxV = 17600., 1800., 277.8
    ro, a = atm(H); M = V/a
    fl5 = 1. if thr > 0.55 else 0.
    Fn1 = (1.351*thr**4-2.1284*thr**3+0.195*thr**2+1.1661*thr+0.418)*maxn1
    FnVol = 6617.6+23.18*V+0.28*H
    flVol = 1. if FnVol > Fn1 else 0.
    tau = 4.
    Fn3 = (-0.85187*(M/0.85)**3+1.7075*(M/0.85)**2+0.14143*(M/0.85)-0.0026990)*640
    Fn4 = (1.1423*thr**3-2.4282*thr**2+0.2863*thr+1.004)*640
    dn1 = (1/tau)*(-n1+(Fn1+0.236*H-Fn3-Fn4+(FnVol-Fn1)*flVol))
    r = n1/maxn1
    # Patev (2004), p. 50, Fig. 45: F_p4 is tabulated as F_p4/640, i.e. its scale is 640 kgf
    # (versions <= 1.1.0 used maxP = 1800 kgf here, which overestimated the thrust by up to 1170 kgf)
    Fp4 = (14.425*r**5-29.718*r**4+22.187*r**3-6.7024*r**2+0.8218*r-0.00005)*640.
    v = V/maxV
    Fp8 = -1. if v > 0.74 else 15.207*v**4-26.538*v**3+13.728*v**2-4.1316*v+0.73637
    Fp8 *= 95
    h = H/18000
    Fp1 = 0. if h > 0.75 else 9.159*h**4-16.24*h**3+8.126*h**2-2.022*h+1.0019
    Fp1 *= 900
    Fp2 = 1.181*v**4-3.382*v**3+3.94*v**2-2.077*v+0.995
    Fp3 = 2.17*r**4-0.102*r**3-1.64*r**2+0.6*r-0.004
    hh = H/12500
    Fp5 = (-1. if hh > 0.8 else 1.007*hh**2-2.4956*hh+0.3635)*455
    Fp6 = 0. if r < 0.8 else 32.32*r**5-67.05*r**4+48.74*r**3-14.42*r**2+1.42*r-0.0001
    Fp7 = 0. if r > 0.8 else 1.
    P = (Fp4+Fp7*Fp8+fl5*(Fp1*Fp2*Fp3+Fp5*Fp6))*9.81
    return dn1, P

def Cmat(psi, th, ga):   # свързана = C * ENU  (формула (4) от статията)
    cp, sp, ct, st, cg, sg = np.cos(psi), np.sin(psi), np.cos(th), np.sin(th), np.cos(ga), np.sin(ga)
    return np.array([[ct*sp, ct*cp, st],
                     [sg*cp-cg*st*sp, -sg*sp-cg*st*cp, cg*ct],
                     [cg*cp+sg*st*sp, -cg*sp+sg*st*cp, -sg*ct]])

def Qmat(b, a):
    return np.array([[np.cos(b)*np.cos(a), np.sin(a), -np.sin(b)*np.cos(a)],
                     [-np.cos(b)*np.sin(a), np.cos(a), np.sin(b)*np.sin(a)],
                     [np.sin(b), 0., np.cos(b)]])

def gravity(B, h):
    return 9.780318*(1+0.0053024*np.sin(B)**2-0.0000058*np.sin(2*B)**2) - 3.086e-6*h

def lim(x, a): return min(max(x, -a), a)

def rhs(x, u, earth=True, clip=True):
    """x = [V, al, be, wx, wy, wz, psi, ga, th, n1, B, L, h]
       u = [throttle, d_ail, d_rud, d_elev] (rad);  w - абсолютни ъглови скорости"""
    V, al, be, wx, wy, wz, psi, ga, th, n1, B, L, h = x
    if clip:
        thr = min(max(u[0], .55), 1.); da = lim(u[1], 17*d2r); dr = lim(u[2], 30*d2r); de = lim(u[3], 30*d2r)
    else:
        thr, da, dr, de = u
    # производните по отклоненията на кормилата (по Пътев) са за 1 градус
    da, dr, de = da/d2r, dr/d2r, de/d2r
    ro, a = atm(h); M = V/a; q = ro*V*V/2
    cy, fak, F = cy_fun(M, al)
    cx = F['cx1'] + 0.09*cy**2
    cz = -0.96*be - dr*(1.4e-3+0.3e-3*M)
    my5 = np.polyval([0.012580,-0.056991,0.0077878,0.14348,0.084738,0.16867], cy)
    mx10 = np.polyval([-0.020544,0.072837,-0.010521,-0.15617,0.026120,0.14757,0.068556,0.070154], cy)
    # страничен канал по Пътев (2004), с. 29-32; m_x^{omega_x} = -0,20 по скоростта на въртене по крен 140 deg/s
    mx5 = np.polyval([0.0094,-0.0299,0.0238,-0.0685], cy)   # Пътев (2004): коефициентите са в единици 10^-2
    mx = da*F['mx1'] - dr*F['mx4'] + be*mx5 + (bref/(2*V))*(-0.20*wx + wy*mx10)
    my = -0.0015*dr - be*F['my3'] + (bref/(2*V))*(my5*wx - 0.27*wy)
    Ra = q*S*np.array([-cx, cy, cz])             # в скоростната к.с. x_a y_a z_a
    Rb = Qmat(be, al) @ Ra                       # Q(beta,alpha): към свързаната к.с.
    dn1, P = engine(n1, V, h, thr)
    C = Cmat(psi, th, ga)
    g = gravity(B, h)
    Vb = V*np.array([np.cos(be)*np.cos(al), -np.cos(be)*np.sin(al), np.sin(be)])
    Wn = C.T @ Vb                                # пътна скорост в ENU (J1 J2 J3)
    Gm = MA/np.sqrt(1-E2*np.sin(B)**2); Qm = Gm*(1-E2)/(1-E2*np.sin(B)**2)
    Gh, Qh = Gm+h, Qm+h
    if earth:
        w_ie = OM*np.array([0., np.cos(B), np.sin(B)])
        w_en = np.array([-Wn[1]/Qh, Wn[0]/Gh, Wn[0]*np.tan(B)/Gh])
    else:
        w_ie = w_en = np.zeros(3)
    w_in_b = C @ (w_ie + w_en)
    w_b = np.array([wx, wy, wz])                 # абсолютни ъглови скорости
    fb = (Rb + np.array([P, 0., 0.]))/m + C @ np.array([0., 0., -g])
    # уравнения на относителното движение: dV/dt = f + g - (w + w_ie)xV
    dVb = fb - np.cross(w_b + C @ w_ie, Vb)
    dV = Vb @ dVb / V
    Vx, Vy, Vz = Vb
    dal = (Vy*dVb[0]-Vx*dVb[1])/(Vx*Vx+Vy*Vy)
    dbe = (V*dVb[2]-Vz*dV)/(V*V*np.cos(be))
    mz = F['mz1'] + cy*(0.048-F['mz5']) - 0.0154*de - (cref/V)*(F['mz11']*wz + F['mz12']*dal) - 0.1*fak
    Mx, My, Mz = q*S*bref*mx, q*S*bref*my, q*S*cref*mz
    dwx = (Mx-(Iz-Iy)*wy*wz)/Ix; dwy = (My-(Ix-Iz)*wz*wx)/Iy; dwz = (Mz-(Iy-Ix)*wx*wy)/Iz
    # кинематика на ойлеровите ъгли: относителни ъглови скорости w - C(w_ie+w_en)
    rx, ry, rz = w_b - w_in_b
    dpsi = (rz*np.sin(ga) - ry*np.cos(ga))/np.cos(th)
    dga = rx + np.sin(th)*dpsi
    dth = ry*np.sin(ga) + rz*np.cos(ga)
    # навигационни уравнения (Release5, BLortgeo.m): dB/dt=W_N/Q, dL/dt=W_E/(G cosB)
    dB = Wn[1]/Qh; dL = Wn[0]/(Gh*np.cos(B)); dh = Wn[2]
    return np.array([dV, dal, dbe, dwx, dwy, dwz, dpsi, dga, dth, dn1, dB, dL, dh]), \
        dict(P=P, cy=cy, cx=cx, M=M, Wn=Wn, q=q, nb=fb/9.80665 - C @ np.array([0,0,-g])/9.80665 + C@np.array([0,0,1.]) )

def rk4(f, x, u, dt, **k):
    k1, _ = f(x, u, **k); k2, _ = f(x+dt/2*k1, u, **k); k3, _ = f(x+dt/2*k2, u, **k); k4, _ = f(x+dt*k3, u, **k)
    return x + dt/6*(k1+2*k2+2*k3+k4)
