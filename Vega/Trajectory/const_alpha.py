import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp,log
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def const_alpha(mission,rocket,stage,trajectory,general,optimization,Data):
    g0=9.810665
    Re=6371000    

    tb=0

    Data[5][-1]=Data[5][-1]-550

    #print("payload", stage[-1].payload)
    #for i in range(rocket.actual_stage,rocket.num_stages):
        #stage[i].payload=stage[i].payload-450

    
    def const(t, y):
        D=0
        return [y[2]*cos(y[3]),
                y[2]*sin(y[3]),
                T*cos(alpha)/y[4]-D/y[4]-g0*sin(y[3])*(Re/(Re+y[1]))**2,
                T/(y[4]*y[2])*sin(alpha)+(y[2]/(Re+y[1])-g0/y[2]*(Re/(Re+y[1]))**2)*cos(y[3]),
                -T/ISP/g0,
                g0*sin(y[3])*(Re/(Re+y[1]))**2,
                D/y[4],
                #T/y[4]*sin(alpha)]
                T/y[4]-T/y[4]*cos(alpha)]
    vA=0
    #trajectory.alpha=[0,0]

    for i in range(0,rocket.num_stages-rocket.actual_stage-1):

        T=stage[i+rocket.actual_stage].thrust
        ISP=stage[i+rocket.actual_stage].Isp


        tb=(stage[i+rocket.actual_stage].propellant_mass)/(T)*g0*ISP
#        print(Data[0][-1],tb)
       
        alpha=trajectory.alpha[i]



        def vel_reach(t,y):
            return y[2]+y[5]+y[6]+y[7]-rocket.actual_DV        
        vel_reach.terminal = True


        sol = solve_ivp(const, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],vA] , events=[vel_reach], t_eval=np.linspace(0,tb,1000))

        
        tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
        Data[0].extend(tm[1:])
        Data[1].extend(sol.y[0][1:])
        Data[2].extend(sol.y[1][1:])
        Data[3].extend(sol.y[2][1:])
        Data[4].extend(sol.y[3][1:])
        Data[5].extend(sol.y[4][1:])
        Data[6].extend(np.ones(len(sol.t)-1)*alpha)
        Data[7].extend(sol.y[5][1:])
        Data[8].extend(sol.y[6][1:])

        #print(Data[0][-1])

        vA=sol.y[7][-1]

        #print(Data[3][-1],Data[7][-1],Data[8][-1],vA,rocket.actual_DV)
        Data[5][-1]=stage[i+rocket.actual_stage].payload-550

        #print((stage[i+rocket.actual_stage+1].total_mass+stage[i+rocket.actual_stage+1].payload),((stage[i+rocket.actual_stage+1].structural_mass+stage[i+rocket.actual_stage+1].payload)))

        stage[i+rocket.actual_stage+1].DV=log((stage[i+rocket.actual_stage+1].total_mass+stage[i+rocket.actual_stage+1].payload-550)/((stage[i+rocket.actual_stage+1].structural_mass+stage[i+rocket.actual_stage+1].payload-550)))*stage[i+rocket.actual_stage+1].Isp*g0
     
        rocket.actual_DV=rocket.actual_DV+stage[i+rocket.actual_stage+1].DV
        #print(rocket.actual_stage, rocket.actual_DV)
        

        
    if trajectory.coast==True:

        T=0
        tb=trajectory.coast_time
        alpha=0
        sol = solve_ivp(const, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],0] , t_eval=np.linspace(0,tb,100))


        tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
        Data[0].extend(tm[1:])
        Data[1].extend(sol.y[0][1:])
        Data[2].extend(sol.y[1][1:])
        Data[3].extend(sol.y[2][1:])
        Data[4].extend(sol.y[3][1:])
        Data[5].extend(sol.y[4][1:])
        Data[6].extend(np.ones(len(sol.t)-1)*alpha)
        Data[7].extend(sol.y[5][1:])
        Data[8].extend(sol.y[6][1:])

    
    
        #print(Data[3][-1]+Data[7][-1]+Data[8][-1]+vA-rocket.actual_DV+stage[i+rocket.actual_stage].DV)

    return Data,vA
