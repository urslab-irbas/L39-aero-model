"""Trimming by Newton's method (trim2) and Runge-Kutta 4 simulation (sim).

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from trim import *
from scipy.optimize import fsolve
def trim2(V=190., H=1000., psi=0.):
    def res(z):
        al,thr,de,n1=z
        x=state(V,al,n1,H,psi); d,o=rhs(x,[thr,0,0,de],clip=False); return [d[0],d[1]*100,d[5]*100,d[9]/100]
    z=fsolve(res,[2*d2r,0.6,1*d2r,15500])
    return state(V,z[0],z[3],H,psi), np.array([z[1],0.,0.,z[2]])
def sim(x0, ufun, T, dt=0.01, earth=True):
    n=int(round(T/dt)); X=np.zeros((n+1,len(x0))); U=np.zeros((n+1,4)); X[0]=x0; x=x0.copy()
    U[0]=ufun(0.,x)
    for k in range(n):
        t=k*dt; u=ufun(t,x); U[k]=u
        x=rk4(rhs,x,u,dt,earth=earth); X[k+1]=x
    U[n]=ufun(n*dt,x)
    return np.arange(n+1)*dt, X, U
