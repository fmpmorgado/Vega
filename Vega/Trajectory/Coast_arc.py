import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def coast_arc(mission,rocket,stage,trajectory,general,optimization,Data):
    g0=9.810665
    Re=6371000

    def coast(t, y):
        D=0
        return [y[2]*cos(y[3]),y[2]*sin(y[3]),T/y[4]-D/y[4]-g0*sin(y[3])*(Re/(Re+y[1]))**2, (y[2]/(Re+y[1])-g0/y[2]*(Re/(Re+y[1]))**2)*cos(y[3]),-T/ISP/g0, g0*sin(y[3])*(Re/(Re+y[1]))**2, D/y[4]]

    T=0
    ISP=stage[-1].Isp
    tb=30

    sol = solve_ivp(coast, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1]] , t_eval=np.linspace(0,tb,1000))
      
    tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
    Data[0].extend(tm[1:])
    Data[1].extend(sol.y[0][1:])
    Data[2].extend(sol.y[1][1:])
    Data[3].extend(sol.y[2][1:])
    Data[4].extend(sol.y[3][1:])
    Data[5].extend(sol.y[4][1:])
    Data[6].extend(np.zeros(len(sol.t)-1))
    Data[7].extend(sol.y[5][1:])
    Data[8].extend(sol.y[6][1:])

    return Data
