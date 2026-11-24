"""Figures 4-6 of the article with English labels.

Fig. 4 - polynomial approximation of the aerodynamic characteristics vs. tabular data;
Fig. 5 - lift coefficient vs. angle of attack and polars; Fig. 6 - step response
(needs step.npz and long.npz, produced by step.py).

Usage:  python figs_en.py

If you use this code, please cite:
    Madzharov, A. (2026). A six-degree-of-freedom flight-dynamics model of the L-39 aircraft with
    polynomial approximation of the aerodynamic characteristics. Aerospace Research in Bulgaria, XX
    (submitted).
"""
from citation import handle_help, print_citation
handle_help(__doc__)
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Liberation Serif','DejaVu Serif'],'font.size':10,'axes.grid':True,'grid.alpha':.35,'lines.linewidth':1.3,'figure.dpi':150})
P=np.polyval
M=np.linspace(0,0.86,200)
# --- Фиг. 4: полиноми по М спрямо табличните стойности [2]
pan=[('$F_{D1}(M)=C_{D0}$',[0,0.4,0.6,0.7,0.75,0.8,0.85],np.array([2.3,2.3,2.4,2.6,2.9,3.8,7.5])*0.01,[4.6071,-8.7777,6.1683,-1.942,0.2661,-0.0122,0.023]),
('$F_{L1}(M)=C_L^{\\alpha}$, 1/rad',[0,0.188,0.47,0.564,0.61,0.658,0.752,0.85],[4.41,4.44,4.5,4.76,4.96,4.96,4.18,2.78],[427.78,-1112.2,1037.3,-427.7,77.007,-4.6132,4.4158]),
('$F_{m1}(M)=C_{m0}$',[0,0.2,0.477,0.608,0.712,0.85],np.array([42.1,42.1,42.1,39.33,30.03,2.77])*0.001,[-0.36451,0.36082,-0.1099,0.010243,0.042108]),
('$F_{m5}(M)=\\bar x_{ac}$',[0,0.07,0.15,0.282,0.37,0.47,0.584,0.658,0.78,0.85],[0.33,0.33,0.33,0.33,0.33,0.32,0.315,0.326,0.35,0.373],[268.56,-858.15,1093.7,-708.80,247.65,-45.782,4.0469,-0.12810,0.32999]),
('$F_{m11}(M)=-C_m^{\\bar\\omega_z}$',[0,0.282,0.47,0.564,0.658,0.752,0.85],[19.02,19.02,19.85,20.66,21.50,20.57,19.67],[911.16,-2187.4,1861.9,-676.29,104.45,-5.2877,19.039]),
('$F_{m12}(M)=-C_m^{\\dot{\\bar\\alpha}}$',[0,0.47,0.564,0.658,0.752,0.85],[6.047,6.047,7.42,6.65,3.62,0.6],[2331.7,-5618.5,4904.8,-1904.4,323.05,-19.044,6.1631])]
fig,ax=plt.subplots(3,2,figsize=(6.3,6.6))
for a,(tt,x,y,p) in zip(ax.flat,pan):
    a.plot(M,P(p,M),'k-',label='polynomial'); a.plot(x,y,'ko',mfc='w',ms=4,label='tabular data')
    a.set_title(tt,fontsize=10); a.set_xlim(0,0.86)
for a in ax[-1]: a.set_xlabel('M')
ax[0,0].legend(fontsize=8,loc='upper left')
fig.tight_layout(); fig.savefig('fig4_en.png',dpi=300); plt.close(fig)
# --- Фиг. 5: c_y(alpha), c_x(alpha) за различни М
al=np.linspace(-10,24,300)*d2r
fig,ax=plt.subplots(1,2,figsize=(6.3,2.7))
for Mi,ls in zip((0.3,0.5,0.6,0.7,0.8),('-','--','-.',':','-')):
    cy=np.array([cy_fun(Mi,a)[0] for a in al]); cx=coeffs(Mi)['cx1']+0.09*cy**2
    lw=0.8 if Mi==0.8 else 1.3
    ax[0].plot(al/d2r,cy,'k',ls=ls,lw=lw,label=f'M={Mi}'); ax[1].plot(cy,cx,'k',ls=ls,lw=lw)
ax[0].set_xlabel('α, °'); ax[0].set_ylabel('$C_L$'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('$C_L$'); ax[1].set_ylabel('$C_D$'); ax[1].set_title('drag polar',fontsize=9)
fig.tight_layout(); fig.savefig('fig5_en.png',dpi=300); plt.close(fig)
# --- Fig. 6: elevator step; a-c: first 120 s / 4 s with the orthogonal-basis model, d-f: 800 s (fixed throttle / throttle law)
d=np.load('step.npz'); t=d['t']; X=d['X']; ny=d['ny']; th=d['th']
k0=int(round(5/0.01)); t=t[k0:]-t[k0]; X=X[k0:]; ny=ny[k0:]; th=th[k0:]
W0=X[0,0]; Tth=25.; thz=3*d2r
thk=thz*(1-np.exp(-t/Tth)); dthk=(thz-thk)/Tth; nk=np.cos(thk)+W0*dthk/9.80665   # orthogonal basis
L=np.load('long.npz'); tl=L['t']; Xf=L['Xf']; Xl=L['Xl']; Uf=L['Uf']; Ul=L['Ul']
fig,ax=plt.subplots(3,2,figsize=(6.3,6.6))
a=ax[0,0]; a.plot(t,X[:,1]/d2r,'k'); a.set_ylabel('α, °'); a.set_title('a) angle of attack',fontsize=9); a.set_xlim(0,t[-1])
a=ax[1,0]; a.plot(t,X[:,8]/d2r,'k',label='ϑ (flight dynamics)'); a.plot(t,th/d2r,'k--',label='θ (flight dynamics)'); a.plot(t,thk/d2r,'k:',lw=1.6,label='θ = ϑ (orthogonal basis)')
a.set_ylabel('°'); a.legend(fontsize=6.5); a.set_title('b) pitch and flight-path angle',fontsize=9); a.set_xlim(0,t[-1])
a=ax[2,0]; z=slice(0,int(4/0.01)+1); a.plot(t[z],ny[z],'k',label='flight dynamics'); a.plot(t[z],nk[z],'k:',lw=1.6,label='orthogonal basis')
a.set_ylabel('$n_{y1}$'); a.legend(fontsize=7); a.set_title('c) $n_{y1}$, t = 0…4 s',fontsize=9); a.set_xlim(0,4); a.set_xlabel('t, s')
a=ax[0,1]; a.plot(tl,Xf[:,0],'k',label='fixed throttle'); a.plot(tl,Xl[:,0],'k--',label='throttle law'); a.set_ylabel('V, m/s'); a.legend(fontsize=7); a.set_title('d) airspeed',fontsize=9)
a=ax[1,1]; a.plot(tl,Xf[:,12],'k'); a.plot(tl,Xl[:,12],'k--'); a.set_ylabel('h, m'); a.set_title('e) altitude',fontsize=9)
a=ax[2,1]; a.plot(tl,Uf[:,0]*100,'k'); a.plot(tl,Ul[:,0]*100,'k--'); a.set_ylabel('$\\delta_T$, %'); a.set_title('f) throttle',fontsize=9); a.set_xlabel('t, s')
for a in ax[:,1]: a.set_xlim(0,tl[-1])
ax[1,0].set_xlabel(''); ax[0,0].set_xlabel('')
fig.tight_layout(); fig.savefig('fig6_en.png',dpi=300); plt.close(fig)
print('ok')
print_citation()
