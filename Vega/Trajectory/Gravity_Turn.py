import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func,pressure_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def Gravity_turn(mission,rocket,stage,trajectory,general,optimization,Data):
    g0=9.810665
    Re=6371000

    def Gravity_flight(t, y):
        if y[1]<100000 and y[1]>0:
            D=1/2*1.225*exp(-y[1]/8440)*Area*Drag_Coeff(y[2]/sqrt(1.4*287.053*Temp_func(y[1])))*y[2]**2
        else:
            D=0

        if pressure_func(y[1])<0: Pa=0
        else: Pa=pressure_func(y[1])
        
        if rocket.actual_stage==0: T_aux=T
        #-stage[0].A*Pa
        else: T_aux=T

        
        return [y[2]*cos(y[3]),y[2]*sin(y[3]),T_aux/y[4]-D/y[4]-g0*sin(y[3])*(Re/(Re+y[1]))**2, (y[2]/(Re+y[1])-g0/y[2]*(Re/(Re+y[1]))**2)*cos(y[3]),dm, g0*sin(y[3])*(Re/(Re+y[1]))**2, D/y[4]]
    
    def gamma_reach(t,y):
        return y[3]

    gamma_reach.terminal=True

    def mass_reach(t,y):
        return y[4]-(stage[0].structural_mass+stage[0].payload)

    mass_reach.terminal = True
    
    def alt_reach(t,y):
        return y[1]-100000

    alt_reach.terminal = True

    DV=0

    for i in range(0,rocket.num_stages):

        if i==rocket.num_stages-1:break

        T=stage[i].thrust
        ISP=stage[i].Isp
        tb=stage[i].propellant_mass/(T)*g0*ISP
        Area=stage[i].diameter**2/4*np.pi

        dm=-T/ISP/g0
        DV=DV+stage[i].DV

        if rocket.actual_stage==0:
            T=stage[0].thrust
            #+pressure_func(Data[2][-1])*stage[0].A

        
        
        def vel_reach(t,y):
            return y[2]+y[5]+y[6]-DV         

        vel_reach.terminal = True
     


        #sol = solve_ivp(Gravity_flight, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1]], events=[vel_reach, gamma_reach] , t_eval=np.linspace(0,tb,1000))
        sol =solve_ivp(Gravity_flight, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1]], events=[mass_reach,alt_reach] , t_eval=np.linspace(0,tb,1000))

        
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

        print(Data[2][-1], Data[0][-1],Data[5][-1])

        
        if Data[2][-1]>99000:
            break

        if rocket.actual_stage==0: rocket.actual_DV=Data[7][-1]+Data[8][-1]+Data[2][-1]

        
        rocket.actual_stage=rocket.actual_stage+1
        rocket.actual_DV=rocket.actual_DV+stage[rocket.actual_stage].DV

        #print(Data[2][-1],rocket.actual_stage)        
        Data[5][-1]=stage[i].payload
        
    
    return Data
        
