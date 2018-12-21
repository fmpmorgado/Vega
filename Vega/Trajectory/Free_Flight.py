import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def free_flight(vA,mission,rocket,stage,trajectory,general,optimization,Data):
    aux_vA=vA
    g0=9.810665
    Re=6371000
    
    tb=trajectory.coast_time
    for i in range(0,rocket.num_stages):
        T=stage[i].thrust
        ISP=stage[i].Isp 
        tb=tb+(stage[i].propellant_mass)/(T)*g0*ISP

    tb=Data[0][-1]+(stage[-1].propellant_mass)/(stage[-1].thrust)*g0*stage[-1].Isp

    
    ##############Ccntrol law for no coast ARC########################
    #### Only for last stage
        
    Obj=1000

        #Boundary Defined, for better convergence in Free flight
    t_lb=-(tb-Data[0][-1])*0.5
    t_ub=0
    aux3=0

    m = GEKKO(remote=False)

    m.time=np.linspace(0,1,100)

    final=np.zeros(len(m.time))
    final[-1]=1
    final=m.Param(value=final)


    #Lagrange multipliers with l1=0

    l2=m.FV(-0.01,lb=trajectory.l2_lb, ub=trajectory.l2_ub)
    l2.STATUS=1

    l3=m.FV(-1, lb=trajectory.l3_lb, ub=trajectory.l3_ub)
    l3.Status=1

    l4=m.Var(value=0)



    rate=m.Const(value=T/ISP/g0)
    #value=(tb-Data[0][-1]+(t_ub+t_lb)/2)
    tf=m.FV(value=(tb-Data[0][-1]+t_lb), ub=tb-Data[0][-1]+t_ub, lb=tb-Data[0][-1]+t_lb)
    tf.STATUS=1


    aux=m.FV(0, lb=-5, ub=5)
    aux.STATUS=1

    
    vD=m.FV(value=Data[8][-1])
    
    #Initial Conditions

    x=m.Var(value=Data[1][-1])
    y=m.Var(value=Data[2][-1])
    vx=m.Var(value=Data[3][-1]*cos(Data[4][-1]))
    vy=m.Var(value=Data[3][-1]*sin(Data[4][-1]), lb=0)
    vG=m.Var(value=Data[7][-1])
    vA=m.Var(value=0)
    gamma=m.Var()
    alpha=m.Var()
    
    t=m.Var(value=0)
    

    m.Equations([l4.dt()/tf==-l2])

    m.Equation(t.dt()/tf==1)

    mrt = m.Intermediate(Data[5][-1]-rate*t)
    srt = m.Intermediate(m.sqrt(l3**2+(l4-aux)**2))


    m.Equations([x.dt()/tf==vx,
                   y.dt()/tf==vy,
                   vx.dt()/tf==(T/mrt*(-l3)/srt-(vx**2+vy**2)/(Re+y)*m.sin(gamma)),
                   vy.dt()/tf==(T/mrt*(-(l4-aux)/srt)+(vx**2+vy**2)/(Re+y)*m.cos(gamma)-g0 *(Re/(Re+y))**2),
                   vG.dt()/tf==g0 *(Re/(Re+y))**2*m.sin(m.atan(vy/vx))])

    m.Equation(gamma==m.atan(vy/vx))
    m.Equation(alpha==m.atan((l4-aux)/l3))

                    
    m.Equation(vA.dt()/tf==T/mrt-T/mrt*m.cos((alpha-gamma)))
    
    #m.Obj(final*vA)

    # Soft constraints
    m.Obj(final*(y-mission.final_altitude)**2)
    m.Obj(final*100*(vx-mission.final_velocity*cos(mission.final_flight_angle))**2)
    m.Obj(final*100*(vy-mission.final_velocity*sin(mission.final_flight_angle))**2)
    #m.Obj(100*tf)
    m.Obj(final*(l2*vy+
                 l3*(T/(Data[5][-1]-rate*t)*(-l3/(m.sqrt(l3**2+(l4-aux)**2)))-(vx**2+vy**2)/(Re+y)*m.sin(gamma))
                 +(l4-aux)*(T/(Data[5][-1]-rate*t)*(-(l4-aux)/(m.sqrt(l3**2+(l4-aux)**2)))+(vx**2+vy**2)/(Re+y)*m.cos(gamma)-g0*(Re/(Re+y))**2)
                 +1)**2)
    
    
    #Options
    m.options.IMODE=6
    m.options.SOLVER=1
    m.options.NODES=3
    m.options.MAX_MEMORY=10

    m.options.MAX_ITER=500
    m.options.COLDSTART=0
    m.options.REDUCE=0
    m.solve(disp=False)



    print("vA",vA.value[-1]+aux_vA)
    print()
    
    """
            
        #print("Obj= ",Obj)
    print("altitude: ", y.value[-1], vx.value[-1], vy.value[-1])
    print("Final mass=",Data[5][-1]-rate.value*t.value[-1])
    print("Gravity= ",vG.value[-1], " Drag= ",vD.value[-1], " alpha= ", vA.value[-1]," Velocity= ", sqrt(vx.value[-1]**2+vy.value[-1]**2)  )
    print("Flight angle=",gamma.value[-1])
    print("")
    print("l2:",l2.value[-1],"l3:", l3.value[-1], "aux:", aux.value[-1])
    print("tf:", tf.value[-1], "tb:",tb-Data[0][-1])
    print("Obj:", Obj)
   # if t_lb==-10.0 or t_ub==10.0:
   #     break               


    """
        


    tm=np.linspace(Data[0][-1],tf.value[-1]+Data[0][-1],100)

    mass=Data[5][-1]


    Data[0].extend(tm[1:])
    Data[1].extend(x.value[1:])
    Data[2].extend(y.value[1:])

    for i in range(1,100):
        Data[3].extend([sqrt(vx.value[i]**2+vy.value[i]**2)])
        Data[4].extend([atan(vy.value[i]/vx.value[i])])
        Data[5].extend([mass-rate.value*t.value[i]])
        Data[6].extend([atan((l4.value[i]-aux.value[i])/l3.value[i])])

    Data[7].extend(vG.value[1:])
    Data[8].extend(vD.value[1:])    

    DV_required=mission.final_velocity+Data[7][-1]+Data[8][-1]+aux_vA+vA.value[-1]


    Obj=m.options.objfcnval+DV_required


    return Data,DV_required,Obj
