import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def pitch(mission,rocket,stage,trajectory,general,optimization,Data):

    alpha_v=20

    g0=9.810665
    Re=6371000

    T=stage[0].thrust
    Mass=rocket.mass
    ISP=stage[0].Isp
    Area=stage[0].diameter**2/4*np.pi

    #The pitch is also considered only in the first stage

    m = GEKKO(remote=False)

    m.time = np.linspace(0,trajectory.pitch_time,101)

    final=np.zeros(len(m.time))
    final[-1]=1
    final=m.Param(value=final)

    tf=m.FV(value=1,lb=0.1,ub=100)
    tf.STATUS=0

    alpha=m.MV(value=0,lb=-0.2,ub=0.1)
    alpha.STATUS=1


    alpha.DMAX=alpha_v*np.pi/180*trajectory.pitch_time/101



    x=m.Var(value=Data[1][-1], lb=0) 
    y=m.Var(value=Data[2][-1], lb=0) 
    v=m.Var(value=Data[3][-1], lb=0)
    phi=m.Var(value=Data[4][-1], ub=np.pi/2,lb=0)
    mass=m.Var(value=Data[5][-1])
    vG=m.Var(value=Data[7][-1])
    vD=m.Var(value=Data[8][-1])


    
    rho=m.Intermediate(rho_func(y))
    cD=m.Intermediate(Drag_Coeff(v/(m.sqrt(1.4*287.053*Temp_func(y)))))
    D=m.Intermediate(1/2*rho*Area*cD*v**2)
    #D=m.Intermediate(Drag_force(rho_func(y),Temp_func(y),Area,v,1,n))
    

    m.Equations([x.dt()/tf==v*m.cos(phi),
                 y.dt()/tf==v*m.sin(phi),
                 v.dt()/tf==T/mass*m.cos(alpha)-D/mass-g0*m.sin(phi)*(Re/(Re+y))**2,
                 ISP*g0*mass.dt()/tf==-T,
                 v*phi.dt()/tf==(v**2/(Re+y)-g0*(Re/(Re+y))**2)*m.cos(phi)+T*m.sin(alpha)/mass,
                 vG.dt()/tf==g0*m.sin(phi)*(Re/(Re+y))**2,
                 vD.dt()/tf==D/mass])

    #m.fix(alpha,100,0)
    m.Obj(final*(phi-trajectory.pitch)**2)
    m.Obj(1e-4*alpha**2)

    m.options.IMODE=6
    m.options.SOLVER=3
    m.options.MAX_ITER=500
    m.solve(disp=False)

    ###Save Data
    tm=np.linspace(Data[0][-1],m.time[-1]+Data[0][-1],101)

    Data[0].extend(tm[1:])
    Data[1].extend(x.value[1:])
    Data[2].extend(y.value[1:])
    Data[3].extend(v.value[1:])
    Data[4].extend(phi.value[1:])
    Data[5].extend(mass.value[1:])
    Data[6].extend(alpha.value[1:])
    Data[7].extend(vG.value[1:])
    Data[8].extend(vD.value[1:])

    print(Data[4][-1])


    return Data
