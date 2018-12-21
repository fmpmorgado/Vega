import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func, pressure_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure



def Vertical_flight(mission,rocket,stage,trajectory,general,optimization,Data):

    
    g0=9.810665
    Re=6371000

    T=stage[0].thrust
    Mass=rocket.mass
    ISP=stage[0].Isp
    Area=stage[0].diameter**2/4*np.pi

    dm=-T/ISP/g0


    T=stage[0].thrust
    #+pressure_func(0)*stage[0].A


    def Vertical(t,y):
        D=1/2*rho_func(y[1])*Area*Drag_Coeff(y[2]/sqrt(1.4*287.053*Temp_func(y[1])))*y[2]**2

        if pressure_func(y[1])<0: Pa=0
        else: Pa=pressure_func(y[1])
        
        T_aux=T
        #-stage[0].A*Pa
        

        return [y[2]*cos(y[3]),y[2]*sin(y[3]),T_aux/y[4]-D/y[4]-g0*sin(y[3])*(Re/(Re+y[1]))**2, 0,dm, g0*sin(y[3])*(Re/(Re+y[1]))**2, D/y[4]]


    def alt_reach(t,y):
        return y[1]-trajectory.vertical_altitude

    alt_reach.terminal=True
    sol = solve_ivp(Vertical, [0, trajectory.vertical_time], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1]], events=[alt_reach], t_eval=np.linspace(0,trajectory.vertical_time,1000))

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

#    Data[4][-1]=trajectory.pitch

    return Data
