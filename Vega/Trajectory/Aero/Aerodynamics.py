from math import sqrt,pi
import numpy as np
from numpy import log10

import gekko as GEKKO



def Drag_Coeff(Mach):
    CD=-3*10**-6*Mach**6 + 0.0002*Mach**5 - 0.0046*Mach**4 + 0.053 * Mach**3 - 0.2806*Mach**2 + 0.6211*Mach + 0.0568
    return CD


"""
def Drag_force(rho,Sref,V,Mach,Knudsen):
    Kn=Knudsen
    Kn_inf = 0.01
    Kn_sup = 10
    if Kn<Kn_inf:
        CD=Drag_Coeff(Mach)
        D = 1/2*rho*CD*Sref*V**2
    elif Kn>Kn_sup:
        CD=Drag_Coeff(Mach)
        D = 1/2*rho*CD*Sref*V**2
    else:
        CD=Drag_Coeff(Mach)
        D = 1/2*rho*CD*Sref*V**2
    return D

"""
	


def Dynamic_pressure(rho,V):
    q=1/2*rho*V**2
    return q

def Knudsen(L,rho,T,m):
    dynamic_viscosity=1.716e-5*(T/273.15)**(3/2)*(273.15+110.4)/(T+110.4)
    boltzmann=1.3806485279e-23
    molecular_mass=28.97/6.023e26
    kn=dynamic_viscosity/(L*rho)*m.sqrt((pi*molecular_mass)/(2*T*boltzmann))
    return kn






def Drag_force(rho,T,Sref,V,L,m):

    if type(rho) is np.float64:
        kn=Knudsen(L,rho,T,m)
    else:
        kn=Knudsen(L,rho,T,m)
        print("Deu")
        print(kn.value)
    kn_inf=0.0146241
    kn_sup=1000*kn_inf
    A=2
    B=0.5113
    CDfm=1.75+sqrt(np.pi)/2*(V/sqrt(2*287.053*T))


    a=sqrt(1.4*287.053*T)
    Mach=V/a
    CD=Drag_Coeff(Mach)

    if kn.value<kn_inf:
        D = 1/2*rho*CD*Sref*V**2
    elif kn>kn_sup:
        D = 1/2*rho*CDfm*Sref*V**2
    else:
        CD = CD + (CDfm-CD)*(1/3*log10(A*kn)+B)
        D = 1/2*rho*CD*Sref*V**2
        
    return D





