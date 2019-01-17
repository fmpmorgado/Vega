#First 80 km -> US standard atmosphere 1976
#Above 80 km -> MSISE-90 Model of Earth's Upper Atmosphere mean solar activity



import numpy as np
from numpy import polyfit,poly1d
import math


def interp(htop,hbottom,h,Ttop,Tbottom,Ptop,Pbottom,rhotop,rhobottom):
    P=Pbottom+(Ptop-Pbottom)/(htop-hbottom)*(h-hbottom)
    T=Tbottom+(Ttop-Tbottom)/(htop-hbottom)*(h-hbottom)
    rho=rhobottom+(rhotop-rhobottom)/(htop-hbottom)*(h-hbottom)

    return T,P,rho


def atmosphere(h):
    
    z=[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
        12000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
        60000, 70000, 80000, 100000, 120000, 140000, 160000, 180000,
        200000, 220000, 240000, 260000, 280000, 300000, 320000, 340000,
        360000, 380000, 400000, 420000, 440000, 460000, 480000, 50000000]

    T=[288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65,
        236.15, 229.65, 223.15, 216.65, 216.65, 216.65, 221.65, 226.65,
        237.05, 251.05, 265.05, 270.65, 245.45, 217.45, 196.65, 184.0160,
        374.9715, 635.5703, 787.5532, 877.6729, 931.2806, 963.2701,
        982.4191, 993.9173, 1000.8427, 1005.0267, 1007.5620, 1009.1030,
        1010.0423, 1010.6166, 1010.9688, 1011.1853, 1011.3190,
        1011.4014, 1011.4526, 1011.4845]

    P=[1.01325*10**+5, 8.98746*10**+4, 7.94952*10**+4, 7.01085*10**+4, 6.16402*10**+4, 5.40199*10**+4,
        4.71810*10**+4, 4.10607*10**+4, 3.55998*10**+4, 3.07425*10**+4, 2.64363*10**+4, 1.93304*10**+4,
        1.20446*10**+4, 5.47489*10**+3, 2.51102*10**+3, 1.17187*10**+3, 5.58924*10**+2, 2.77522*10**+2,
        1.43135*10**+2, 7.59448*10**+1, 2.03143*10**+1, 4.63422, 8.86280*10**-1, 2.81*10**-2,
        2.17*10**-3, 7.03*10**-4, 3.31*10**-4, 1.80*10**-4, 1.05*10**-4, 6.44*10**-5, 4.09*10**-5,
        2.66*10**-5, 1.77*10**-5, 1.20*10**-5, 8.20*10**-6, 5.69*10**-6, 3.98*10**-6, 2.81*10**-6,
        2.01*10**-6, 1.44*10**-6, 1.04*10**-6, 7.55*10**-7, 5.53*10**-7, 4.07*10**-7]

    rho=[1.22500, 1.11164, 1.00649, 9.09122*10**-1, 8.19129*10**-1, 7.36116*10**-1, 6.59697*10**-1, 5.89501*10**-1, 5.25168*10**-1,
         4.66348*10**-1, 4.12707*10**-1, 3.10828*10**-1, 1.93674*10**-1, 8.80349*10**-2, 3.94658*10**-2, 1.80119*10**-2, 8.21392*10**-3,
         3.85101*10**-3, 1.88129*10**-3, 9.77525*10**-4, 2.88321*10**-4, 7.42430*10**-5, 1.57005*10**-5, 5.08*10**-7, 1.80*10**-8,
         3.26*10**-9, 1.18*10**-9, 5.51*10**-10, 2.91*10**-10, 1.66*10**-10, 9.91*10**-11, 6.16*10**-11, 3.94*10**-11, 2.58*10**-11,
         1.72*10**-11, 1.16*10**-11, 7.99*10**-12, 5.55*10**-12, 3.89*10**-12, 2.75*10**-12, 1.96*10**-12, 1.40*10**-12, 1.01*10**-12, 7.30*10**-13]

    for i in range(0,len(z)-1):
        if h>=z[i] and h<=z[i+1]:
            Data=interp(z[i+1],z[i],h,T[i+1],T[i],P[i+1],P[i],rho[i+1],rho[i])  

    return Data

def Temp_func(h):

    z=[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
            12000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
            60000, 70000, 80000, 100010]

    T=[288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65,
            236.15, 229.65, 223.15, 216.65, 216.65, 216.65, 221.65, 226.65,
            237.05, 251.05, 265.05, 270.65, 245.45, 217.45, 196.65, 184.0160]

    A=polyfit(z,T,5) #8
    p = poly1d(A)

    return p(h)

def rho_func(h):

    z=[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
            12000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
            60000, 70000, 80000, 100010]

    rho=[1.22500, 1.11164, 1.00649, 9.09122*10**-1, 8.19129*10**-1, 7.36116*10**-1, 6.59697*10**-1, 5.89501*10**-1, 5.25168*10**-1,
         4.66348*10**-1, 4.12707*10**-1, 3.10828*10**-1, 1.93674*10**-1, 8.80349*10**-2, 3.94658*10**-2, 1.80119*10**-2, 8.21392*10**-3,
         3.85101*10**-3, 1.88129*10**-3, 9.77525*10**-4, 2.88321*10**-4, 7.42430*10**-5, 1.57005*10**-5, 5.08*10**-7]


    A=polyfit(z,rho,8)
    p = poly1d(A)



    return p(h)

def Pfunc(h):

    z=[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
        12000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
        60000, 70000, 80000, 100010]

    P=[1.01325*10**+5, 8.98746*10**+4, 7.94952*10**+4, 7.01085*10**+4, 6.16402*10**+4, 5.40199*10**+4,
        4.71810*10**+4, 4.10607*10**+4, 3.55998*10**+4, 3.07425*10**+4, 2.64363*10**+4, 1.93304*10**+4,
        1.20446*10**+4, 5.47489*10**+3, 2.51102*10**+3, 1.17187*10**+3, 5.58924*10**+2, 2.77522*10**+2,
        1.43135*10**+2, 7.59448*10**+1, 2.03143*10**+1, 4.63422, 8.86280*10**-1, 2.81*10**-2]

    A=polyfit(z,P,8)
    p = poly1d(A)

    return p(h)


"""


import matplotlib.pyplot as plt
import numpy as np


z=[0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000,
            12000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
            60000, 70000, 80000, 100010]

T=[288.15, 281.65, 275.15, 268.65, 262.15, 255.65, 249.15, 242.65,
            236.15, 229.65, 223.15, 216.65, 216.65, 216.65, 221.65, 226.65,
            237.05, 251.05, 265.05, 270.65, 245.45, 217.45, 196.65, 184.0160]

A=polyfit(z,T,5) #8
p = poly1d(A)

aux=np.linspace(0,100010,100000)

plt.plot(z, T, '.', aux,p(aux), '-')
plt.show()

"""
