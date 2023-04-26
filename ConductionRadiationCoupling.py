/********************************************************
 * Condution-Radiation Coupling  March-2023
 * Authors: David Becerril Rodriguez & Raul Esquivel-Sirvent
 * Affilitation: CNR-ISM & UNAM
 *
 * Solves the thermal conduction equations in planar bodies
 * coupled to the thermal radiation of the bodies 
 *************************************************************/
 
 
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.interpolate import interp1d
from scipy import interpolate 
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import os
from scipy.integrate import quad
import time 

c    = 2.99792458e8     # Velocidad de la luz en el vacío en m/s
hbar = 6.62606957e-34/(2.0*np.pi)      #Cte de Planck J*s dividida por 2Pi.
kB   = 1.380658e-23     #Cte de Boltzmann en J/K
eps0 = 1.0          #permitividad del vacío
eVtoRads = 1.5193e15  #Factor de conversión de eV a rad/s*)
cminvtoRads = 1.88363e11 #Factor de conversión de cm-1 a rad/s*)
w0=1.0e14
c_cm = c*1e2


## Transform from wavlength in cm 
## to angular frequency. 
def wtaf(wvl_cm):
    return 2*np.pi*c_cm/wvl_cm
    
def wtf(wvl_cm):
    return c_cm/wvl_cm


def wl_t_wn(wvl_cm):
    return 1/wvl_cm


# Just to test for now
def SiC(wvl_cm):
    omegaL = 1.82652e14#969
    omegaT = 1.49477e14#793
    gamma  = 8.9724e11#4.76
    Einf   = 6.7
    wn     = wtaf(wvl_cm) 
    return Einf*(1.0 + (omegaL**2 - omegaT**2) /(omegaT**2-wn**2 - 1j*gamma*wn) )
    #return Einf*(omegaL**2-wn**2-1j*wn*gamma)/(omegaT**2-wn**2-1j*wn*gamma)
    
    
Diel  = SiC
###########################
## Distance is in cm
## 1e-4 is 1 micron,
## 0.1e-4 are 100 nm, etc.
## 1e-7 is a nanometer
### slab width in cm


#specific heat
kappa = 120*1.0e-2  #(W/m*K)
#kappa  = 490e-3
cv   =  3200*600*(1.0e-6) #J/(m².K) #2330e-3*710      #J⋅K−1⋅kg−1
alpha = kappa/cv

#Initial temperatures
T10 = 400
T20 = 300


def Rp(wvl_cm,beta,dwidth):
    d = dwidth
    k0 = 2*np.pi/wvl_cm      #Componente perpendicular del vector de onda en la Película de YBCO (temperatura T1) para pol. p.
    kz0 = np.sqrt(k0**2 - np.abs(beta)**2 + 0*1j)
    kz1 = np.sqrt(Diel(wvl_cm)*k0**2 - np.abs(beta)**2 + 0*1j)
  
    if (kz1.imag)<0.0e0:
        kz1 = -kz1
    else:
        kz1 = kz1
    
    rp0P = (Diel(wvl_cm)*kz0- kz1)/(Diel(wvl_cm)*kz0 + kz1)            #Coeficiente de reflexión en la interfaz 1 (medio_eps0|eps_ybco T1) para la pol. P
    rpPS = (-Diel(wvl_cm)*kz0 + kz1)/(Diel(wvl_cm)*kz0 + kz1)  
        #Coeficiente de reflexión en la interfaz 2 (eps_ybco T1|sustrato_fepsS) para la pol. P
    return (rp0P+rpPS*np.exp(2.0*1j*kz1*d))/(1.0e0+rp0P*rpPS*np.exp(2.0*1j*kz1*d))              #Coeficiente de reflexión de de sistema eps0/pelicula_T1/sustrato


def Rs(wvl_cm,beta,dwidth):
    d  = dwidth
    k0 = 2*np.pi/wvl_cm      #Componente perpendicular del vector de onda en la Película de YBCO (temperatura T1) para pol. p.
    kz0 = np.sqrt(k0**2 - np.abs(beta)**2 + 0*1j)
    kz1 = np.sqrt(Diel(wvl_cm)*k0**2 - np.abs(beta)**2 + 0*1j)
  
    if (kz1.imag)<0.0e0:
        kz1 = -kz1
    else:
        kz1 = kz1
  
    rs0P = (kz0-kz1)/(kz0+kz1)     #Coeficiente de reflexión en la interfaz 1 (medio_eps0|eps_ybco T1) para la pol. s
    rsPS = (-kz0+kz1)/(kz0+kz1)      #Coeficiente de reflexión en la interfaz 2 (eps_ybco T1|sustrato_fepsS) para la pol. s
    return (rs0P+rsPS*np.exp(2.0*1j*kz1*d))/(1.0e0+rs0P*rsPS*np.exp(2.0*1j*kz1*d))        #Coeficiente de reflexión de de sistema eps0/pelicula_T1/sustrato

def tau_prop(wvl_cm,beta,alphap,dwidth,dgap):
    rr = 0
    L0 = dgap
    k0 = 2*np.pi/wvl_cm

    kz0 = np.sqrt(k0**2 - np.abs(beta)**2 + 0*1j)

    if np.imag(kz0) != 0:
        return 0
    else:
        if alphap == 's':
            rr = Rs(wvl_cm,beta,dwidth)
        elif alphap == 'p':
            rr = Rp(wvl_cm,beta,dwidth)

        return (1-np.abs(rr)**2)**2/np.abs(1-rr*rr*np.exp(2*1j*kz0*L0))**2

        
def tau_evan(wvl_cm,beta,alphap,dwidth,dgap):
    L0 = dgap
    rr = 0

    k0 = 2*np.pi/wvl_cm

    kz0 = np.sqrt(k0**2 - np.abs(beta)**2 + 0*1j)
    if np.imag(kz0) == 0:

        return 0.0
    else:

        if alphap == 's':
            rr = Rs(wvl_cm,beta,dwidth)
        elif alphap == 'p':
            rr = Rp(wvl_cm,beta,dwidth)
        
        return (4*rr.imag**2*np.exp(-2*np.abs(kz0)*L0))/np.abs(1.0-rr*rr*np.exp(-2.0*np.abs(kz0)*L0))**2


def ThetaT(T,wvl_cm):
    w = wtaf(wvl_cm)
    x = hbar*w
    return x/(np.exp( x/(kB*T) )-1.0 )

def dThetaT(T,wvl_cm):
    w = wtaf(wvl_cm)
    x = kB*T
    cc = hbar*w
    xx = cc/x
    return kB*xx**2*np.exp(xx)/(np.exp(xx)-1)**2



### Integral over the wavenumbers.
def sw(wvl_cm,ki,kf,dwidth,dgap):
    
    def fkernel(beta,wvl_cm,dwidth,dgap):
        out = (beta/(2*np.pi)**2)*(tau_evan(wvl_cm,beta,'p',dwidth,dgap)  + tau_prop(wvl_cm,beta,'p',dwidth,dgap))
        out += (beta/(2*np.pi)**2)*(tau_evan(wvl_cm,beta,'s',dwidth,dgap)  + tau_prop(wvl_cm,beta,'s',dwidth,dgap))
        
        return out
    
    return quad(fkernel, ki, kf, args=(wvl_cm,dwidth,dgap))[0]

def makeKlist(x_angfreq,angfreq_range):
    kfinal = 1e4
    if angfreq_range[0] < x_angfreq < angfreq_range[1]:
        klist0 = np.arange(0.0,10.0*x_angfreq/c_cm,angfreq_range[2]*x_angfreq/c_cm)
        return np.append(klist0,[10*x_angfreq/c_cm,kfinal*x_angfreq/c_cm] )
    else:
        return np.array([0.0,kfinal*x_angfreq/c_cm])
    
def swConverged(wvl_cm,dwidth,dgap,angfreq_range):
    
    def fkernel(beta,wvl_cm,dwidth,dgap):
        out = (beta/(2*np.pi)**2)*(tau_evan(wvl_cm,beta,'p',dwidth,dgap)  + tau_prop(wvl_cm,beta,'p',dwidth,dgap))
        out += (beta/(2*np.pi)**2)*(tau_evan(wvl_cm,beta,'s',dwidth,dgap)  + tau_prop(wvl_cm,beta,'s',dwidth,dgap))
        
        return out
    
    out   = 0.0
    omega = wtaf(wvl_cm)
    klist = makeKlist(omega,angfreq_range)
    
    for i,ki in enumerate(klist[:-1]):
        out += quad(fkernel, klist[i], klist[i+1], args=(wvl_cm,dwidth,dgap),limit = angfreq_range[3])[0]
    
    return out



def G_Rad(T,omegai,omegaf,dwidth,dgap):

    def ff(x_angfreq,dwidth,dgap):
        kf = 1e4*x_angfreq/c_cm
        return sw(wtaf(x_angfreq),0.0,kf,dwidth,dgap)*dThetaT(T,wtaf(x_angfreq))
        #return sw(wtaf(x_angfreq),ki,kf)*( ThetaT(T10,wtaf(x_angfreq)) - ThetaT(T20,wtaf(x_angfreq))  )

    return   quad(ff, omegai, omegaf, args = (dwidth,dgap) )[0]


def G_RadConverged(T,omegai,omegaf,dwidth,dgap,angfreq_range):

    def ff(x_angfreq,dwidth,dgap,T,angfreq_range):
        return swConverged(wtaf(x_angfreq),dwidth,dgap,angfreq_range)*dThetaT(T,wtaf(x_angfreq))

    return   quad(ff, omegai, omegaf, args = (dwidth,dgap,T,angfreq_range) )[0]

def GRadSpline(*args):
    
    omegai = args[0]
    omegaf = args[1]
    Nomega = args[2]
    
    dwidth = args[3]
    dgap   = args[4]
    
    convg_omegai = args[5]
    convg_omegaf = args[6]
    steps        = args[7]
    limitN       = args[8]
    showplt      = args[9]
    
    angfreq_convg = [convg_omegai,convg_omegaf,steps,limitN]
    omegas = np.linspace(omegai,omegaf,Nomega)
    swlist = np.zeros(Nomega)
    
    for i,omegaii in enumerate(omegas):
        lambdai = wtaf(omegaii)
        kfi     = 1e4*omegaii/c_cm
        swlist[i] = swConverged(lambdai,dwidth,dgap,angfreq_convg)*dThetaT(T20,lambdai)  
    
    spl = interpolate.interp1d(omegas,swlist)
 
    if showplt == True:
        fig, ax = plt.subplots()
        ax.plot(omegas, swlist,'-o', lw = 2, label = "values", )
        ax.plot(omegas,spl(omegas), label = "spline")
        ax.set(xlabel='Ang. Frequency ](Hz)', ylabel='Sw',yscale = 'log')
        plt.legend(loc='best')
        plt.show()

    integral1 = quad(spl,omegai,omegaf,limit = limitN)[0]

    return  integral1

def showTaus(lambdai,lambdaf,betai,betaf,dwidth,dgap):
    Nk = 80
    Nw = 100
    
    lambdas = np.linspace(lambdai,lambdaf,Nw)
    omegas = np.linspace(wtaf(lambdaf),wtaf(lambdai),Nw)
    ks      =np.linspace(betai,betaf,Nk)
    out = np.zeros((Nw,Nk))

    for i,wi in enumerate(omegas):
        for j,ki in enumerate(ks):
            out[i,j] =  tau_evan(wtaf(wi),ki,'p',dwidth,dgap)   + tau_prop(wtaf(wi),ki,'p',dwidth,dgap)


    fig, ax = plt.subplots()
    y = omegas
    x = ks
    z = out
    cc = ax.pcolormesh( ks,omegas, z, cmap='viridis', vmin=z.min(), vmax=z.max())
    ax.set(xlabel='wavenumber(cm^-1)', ylabel='Angular Frequency (Hz)')

    ax.set_title('Taus')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(cc, ax=ax)
    plt.yscale('log')
    plt.show()

def showSw(lambdai,lambdaf,ki,kf,units,dwidth,dgap):
    Nw      = 100
    Nk      = 80
    swdata  = np.zeros(Nw)

 
    fig, ax = plt.subplots()
    
    if units == "af":
        omegas = np.linspace(wtaf(lambdaf),wtaf(lambdai),Nw)
        for i in range(Nw):
            lambdai = wtaf(omegas[i])
            kfi = 1e3*omegas[i]/c_cm
            swdata[i] = sw(lambdai,0,kfi)*( dThetaT(T20,lambdai)  )
            #swdata[i] = sw(lambdai,0,kfi,dwidth,dgap)#*( ThetaT(T10,lambdai) - ThetaT(T20,lambdai)  )

        ax.plot(omegas, swdata/1e2,'-o')
        plt.yscale('log')
        ax.set(xlabel='Ang. Frequency (Hz)', ylabel='Sw')
    elif units == "wl":
        lambdas = np.linspace(lambdai,lambdaf,Nw)

        for i in range(Nw):
            #swdata[i] = sw(lambdas[i],ki,kf)*( dThetaT(T20,lambdas[i])  )
            swdata[i] = sw(lambdas[i],ki,kf,dwidth,dgap)#*( ThetaT(T10,lambdas[i]) - ThetaT(T20,lambdas[i])  )

    
        ax.plot(lambdas, swdata,'-o')
    ax.grid()
    plt.show()
    
def get_rsFourier(xnvec,dwidth,GG):
    return xnvec*alpha/dwidth**2


def TempDistF(z,t,xns,dwidth):
    
    xn    = xns
    nterms = len(xn)
    dT    = T10-T20


    T1n = 0.0 + 1j*0.0
    T2n = 0.0 + 1j*0.0
 

    for i in range(nterms):

        coef0 = 8*dT    
        
        T1n +=  coef0*np.sin(xn[i])*np.cos(xn[i])**2/(4*xn[i]+np.sin(4*xn[i]) )*np.exp(-(xn[i]**2*alpha/dwidth**2)*t)*np.cos(xn[i]*(z+dwidth)/dwidth) 
        T2n += -coef0*np.sin(xn[i])**2*np.cos(xn[i])/(4*xn[i]+np.sin(4*xn[i]) )*np.exp(-(xn[i]**2*alpha/dwidth**2)*t)*np.sin(xn[i]*z/dwidth) 

    return T20+ T1n.real,T20 + T2n.real
def get_xn(Nroots,gg,dwidth):    
    f0 = lambda x:x*math.tan(2*x) - (2*dwidth*gg/kappa)
    derr = 1e-12
    out = []
    for ni in range(Nroots):
        xi = 0.5*ni*np.pi 
        out.append( brentq(f0,xi, xi + 0.25*np.pi-derr )  )
     
    return  np.array(out)



def get_rs(xnvec,tau,dwidth,GG):
    v = np.sqrt(alpha/tau)
    
    out = 1j*np.zeros( [len(xnvec),2] )
    for i in range(len(xnvec) ):
        gamman = (xnvec[i]/dwidth)**2 
        
        rm = -v**2/(2*alpha)*(1.0 + np.sqrt(1.0-4.0*gamman*(alpha/v)**2+0j ) ) 
        rp = -v**2/(2*alpha)*(1.0 - np.sqrt(1.0-4.0*gamman*(alpha/v)**2 +0j) )

        out[i,0] = rp
        out[i,1] = rm
    
    return out


def get_un(xnvec,dwidth):
    return np.sin(xnvec)/xnvec

def get_vn(xnvec,dwidth):
    return (4*xnvec + np.sin(4*xnvec) )/(8*xnvec*np.cos(xnvec)**2)

def get_wn(xnvec,dwidth):
    return -np.cos(xnvec) + np.tan(xnvec)*np.sin(xnvec)


def TempDist(z,t,xns,tau,GG,dwidth,dgap):
    
    xn    = xns
    nterms = len(xn)
    dT    = T10-T20
    rs = get_rs(xn,tau,dwidth,GG)
    r1 = rs[:,0]
    r2 = rs[:,1]
    
    un = get_un(xn,dwidth)
    vn = get_vn(xn,dwidth)
    wn = get_wn(xn,dwidth)

    T1n = 0.0 + 1j*0.0
    T2n = 0.0 + 1j*0.0
 
    an = 0.0
    bn = 0.0
    
    for i in range(nterms):
        an = -r2[i]*(un[i]/vn[i]) + (GG*wn[i])/(dwidth*cv*vn[i])
        bn = r1[i]*(un[i]/vn[i]) - (GG*wn[i])/(dwidth*cv*vn[i])
        coef0 = dT/(r1[i]-r2[i])       
        
        T1n += coef0*(an*np.exp(r1[i]*t) + bn*np.exp(r2[i]*t) )*np.cos(xn[i]*(z+dwidth)/dwidth)
        T2n += -coef0*np.tan(xn[i])*(an*np.exp(r1[i]*t) + bn*np.exp(r2[i]*t) )*np.sin(xn[i]*(z-dwidth-dgap)/dwidth)

    return T20+ T1n.real,T20 + T2n.real



