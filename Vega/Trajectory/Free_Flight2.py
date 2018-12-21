import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp,log,atan2
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure

def free_flight2(vA,mission,rocket,stage,trajectory,general,optimization,Data):

    aux_vA=vA
    aux_g=Data[7][-1]
    
    g0=9.810665
    Re=6371000

    Ue=398604.3*10**9


    #print(vA) 
    T=stage[-1].thrust
    ISP=stage[-1].Isp
    tb=(stage[-1].propellant_mass)/(T)*g0*ISP
    #print(tb)


#    print(Data[0][-1],tb)
    def vel_reach(t,y):
        return y[2]+y[5]+y[6]+y[7]-rocket.actual_DV        
    vel_reach.terminal = True    


    
    def traj2(t,y):
        """

        #gamma=atan(y[3]/y[2])
        alpha=atan(y[10]/y[9])
        #V=sqrt((T/y[4]*cos(alpha-y[3])-(y[2])/(Re+y[1])*sin(y[3]))**2+(T/y[4]*sin(alpha-y[3])+(y[2])/(Re+y[1])*cos(y[3])-g0 *(Re/(Re+y[1]))**2)**2)
        V=T*cos(alpha)/y[4]-g0*sin(y[3])*(Re/(Re+y[1]))**2
        #gamma=T/(y[4]*y[2])*sin(alpha-y[3])+(y[2]/(Re+y[1])-g0/y[2]*(Re/(Re+y[1]))**2)*cos(y[3])
        gamma=T/(y[4]*y[2])*sin(alpha)+(y[2]/(Re+y[1])-g0/y[2]*(Re/(Re+y[1]))**2)*cos(y[3])
        return [y[2]*cos(y[3]),
                y[2]*sin(y[3]),
                #T/y[4]*(-y[9])/(sqrt(y[9]**2+y[10]**2))-(y[2]**2+y[3]**2)/(Re+y[1])*sin(gamma),
                #T/y[4]*(-y[10])/(sqrt(y[9]**2+y[10]**2))+(y[2]**2+y[3]**2)/(Re+y[1])*cos(gamma)-g0 *(Re/(Re+y[1]))**2,
                V,
                gamma,
                -T/ISP/g0,
                g0 *(Re/(Re+y[1]))**2*sin(y[3]),
                0,
                T/y[4]-T/y[4]*cos(alpha),
                0,
                0,
                -y[8]]


        """


        x=y[0]
        ay=y[1]
        V=y[2]
        gamma=y[3]
        m=y[4]

        lr=y[8]
        lg=y[9]
        lv=y[10]


        sa=-lg/(V*sqrt((lg/V)**2+lv**2))
        ca=-lv/sqrt((lg/V)**2+lv**2)


        g=Ue/(ay+Re)**2
        TM=T/m
        r=ay+Re
        
        dx=V*cos(gamma)
        dy=V*sin(gamma)
        dV=TM*ca-g*sin(gamma)
        dgamma=TM*sa/V+(V/r-g/V)*cos(gamma)
        dm=-T/ISP/g0
        dvg=g*sin(gamma)
        dD=0
        dva=TM-TM*ca

        dlr=V*lg*cos(gamma)/r**2-(2*Ue*lv*sin(gamma)+2*Ue*lg*cos(gamma)/V)/r**3
        dlg=-lr*cos(gamma)*V+Ue*lv*cos(gamma)/r**2+lg*sin(gamma)*(V/r-Ue/(r**2*V))
        dlv=-lr*sin(gamma)-lg*(cos(gamma)*(1/r+Ue/((r**2)*(V**2)))-TM/V**2*sa)

        #print(dlr,dlv,dlg)
        
        return [dx,dy,dV,dgamma,dm,dvg,dD,dva,dlr,dlg,dlv]
    #x,y,v,gamma,m,vg,D,va,l2,l3,l4
    tf=trajectory.aux*tb

    #print(tf)
    sol = solve_ivp(traj2, [0, tf], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],vA,trajectory.l[0],trajectory.l[1],trajectory.l[2]] , events=[vel_reach],t_eval=np.linspace(0,tf,100))

    
    vA=sol.y[7][-1]



    #print(vA)
    tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
    Data[0].extend(tm[1:])
    Data[1].extend(sol.y[0][1:])
    Data[2].extend(sol.y[1][1:])


    for i in range(1,len(sol.y[0])):
        
        sa=-sol.y[9][i]/(sol.y[2][i]*sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
        ca=-sol.y[10][i]/(sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))

        
       # if atan2(ca,sa)<-np.pi/2: Data[6].extend([atan2(ca,sa)+np.pi])
       # elif atan2(ca,sa)>np.pi/2: Data[6].extend([atan2(ca,sa)-np.pi])
       # else: Data[6].extend([atan2(ca,sa)])

        Data[6].extend([atan2(ca,sa)])
        
        #print(atan2(ca,sa),sol.y[9][i],sol.y[10][i])
    Data[3].extend(sol.y[2][1:])
    Data[4].extend(sol.y[3][1:])
    Data[5].extend(sol.y[4][1:])
    Data[7].extend(sol.y[5][1:])
    Data[8].extend(sol.y[6][1:])

    #print(rocket.actual_DV)
    #print(Data[3][-1]+Data[8][-1]+Data[7][-1]+vA-rocket.actual_DV)

    """
    y0=X
    y1=y
    y2=vx
    y3=vy
    y4=mass
    y5=vG
    y6=Drag
    y7=vA
    y8=l2
    y9=l3
    y10=l4


    """


    #print(Data[0][-1])
    DV_required=0
    Obj=0
    return Data,DV_required,Obj
