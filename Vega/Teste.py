from scipy.integrate import solve_ivp
import pygmo as pg
from math import sin,cos,atan,pi,sqrt
import numpy as np
import matplotlib.pyplot as plt



TM=30
Ue=398604.3*10**9
    
g0=9.810665
Re=6371000

R0=Re
V0=30
Gamma0=0.5

Robj=500000+Re
Gammaobj=0.3
Vobj=5000

g0=9.810665


x0=0
y0=0
V0=100
ga0=np.pi/2


def traj(t,y):
    """
    if y[5]==0: y[5]=0.00000001

    alpha=atan(y[4]/y[5]/y[2])

    r=y[2]*sin(y[1])
    gamma=TM*sin(alpha)/y[2]+(y[0]*y[2]**2-Ue)/(y[0]**2*y[2])*cos(y[1])
    V=TM*cos(alpha)-(Ue/y[0]**2)*sin(y[1])

    l1=-y[4]*(-y[0]*y[2]**2+2*Ue)/(y[0]**3*y[2])*cos(y[1])-y[5]*(2*Ue/(y[0]**3))*sin(y[1])
    l2=-y[3]*y[2]*cos(y[1])-y[4]*Ue*sin(y[1])/(y[0]**2*y[2])+y[4]*y[2]/y[0]*sin(y[1])+y[5]*Ue/y[0]**2*cos(y[1])
    l3=-y[3]*sin(y[1])-y[4]*Ue*cos(y[1])/(y[0]**2*y[2]**2)-y[4]*cos(y[1])/y[0]+y[4]*TM*sin(alpha)/y[2]**2



    return [r,gamma,V,l1,l2,l3]

    """

    x=y[0]
    ay=y[1]
    V=y[2]
    g=y[3]

    lx=y[4]
    ly=y[5]
    lv=y[6]
    lg=y[7]


    if lv==0: lv=0.00000000000000001
    a=atan(lg/(V*lv))
    

    sa=-lg/(V*sqrt((lg/V)**2+lv**2))
    ca=-lv/sqrt((lg/V)**2+lv**2)

    dx=V*cos(g)
    dy=V*sin(g)
    dV=TM*ca
    dg=TM*sa/V

    dlx=0
    dly=0
    dlv=-lx*cos(g)-ly*sin(g)+lg*(TM*sa)/V**2
    dlg=lx*V*sin(g)-ly*cos(g)*V
    
    return[dx,dy,dV,dg,dlx,dly,dlv,dlg]

def model(x):

    """
    P_l1=x[0]
    P_l2=x[1]
    P_l3=x[2]
    tf=x[3]
    
    
    
    sol = solve_ivp(traj, [0, tf], [R0,Gamma0,V0,P_l1,P_l2,P_l3] , events=[],t_eval=np.linspace(0,tf,1000))
    

    Rf=sol.y[0][-1]
    Vf=sol.y[2][-1]
    Gammaf=sol.y[1][-1]


    l1f=sol.y[3][-1]
    l2f=sol.y[4][-1]
    l3f=sol.y[5][-1]

    alphaf=atan(l2f/(Vf*l3f))


    Hf=l1f*(Vf*sin(Gammaf))+l2f*(TM*sin(alphaf)/Vf+(Vf**2*Rf-Ue)/(Rf**2*Vf)*cos(Gammaf))+l3f*(TM*cos(alphaf)-Ue/Rf**2*sin(Gammaf))


    Obj=10*abs(Rf-Robj)+10*3*abs(Vf-Vobj)+10**6*abs(Gammaf-Gammaobj)
    #+abs(Hf+1)

    print("h=",Rf-Re,"V=",Vf,"Gamma",Gammaf,"Obj=", Obj)


    

    """
    P_lx=0
    P_ly=x[0]
    P_lv=x[1]
    P_lg=x[2]
    tf=x[3]
    

    sol = solve_ivp(traj, [0, tf], [x0,y0,V0,ga0,P_lx,P_ly,P_lv,P_lg] , events=[],t_eval=np.linspace(0,tf,1000))



    xf=sol.y[0][-1]
    yf=sol.y[1][-1]
    Vf=sol.y[2][-1]
    gf=sol.y[3][-1]

    lxf=sol.y[4][-1]
    lyf=sol.y[5][-1]
    lvf=sol.y[6][-1]
    lgf=sol.y[7][-1]



    af=atan(lgf/lvf/Vf)

    saf=-lgf/(Vf*sqrt((lgf/Vf)**2+lvf**2))
    caf=-lvf/sqrt((lgf/Vf)**2+lvf**2)


    Hf=lxf*Vf*cos(gf)+lyf*Vf*sin(gf)+lvf*TM*caf+lgf*TM*saf/Vf

    Obj=10**6*abs(gf)+abs(yf-(Robj-Re))

    #print(Obj)
    if abs(Obj)<10000: print(Obj,tf,gf,Vf,yf,Hf)

    print(Obj,gf)
    return abs(Obj)



class Optimization:
    def fitness(self, x):
        fit=model(x)
        return [fit]

    def get_bounds(self):
        return ([-1,-1,-1,250],[1,1,1,400]) 
        




prob = pg.problem(Optimization())
algo = pg.algorithm(pg.pso_gen(gen=300))          #Choose of heuristic algorithm and number of generations
pop = pg.population(prob,100)                        #Choose number of individuals
pop = algo.evolve(pop)                              #Evolve the population
y=pop.champion_x                                    #Extract the best DV division solution to minimize mass



P_lx=0
P_ly=y[0]
P_lv=y[1]
P_lg=y[2]
tf=y[3]


sol = solve_ivp(traj, [0, tf], [x0,y0,V0,ga0,P_lx,P_ly,P_lv,P_lg] , events=[],t_eval=np.linspace(0,tf,1000))



xf=sol.y[0][-1]
yf=sol.y[1][-1]
Vf=sol.y[2][-1]
gf=sol.y[3][-1]

lxf=sol.y[4][-1]
lyf=sol.y[5][-1]
lvf=sol.y[6][-1]
lgf=sol.y[7][-1]



af=atan(lgf/lvf/Vf)

saf=-lgf/(Vf*sqrt((lgf/Vf)**2+lvf**2))
caf=-lvf/sqrt((lgf/Vf)**2+lvf**2)

Hf=lxf*Vf*cos(gf)+lyf*Vf*sin(gf)+lvf*TM*caf+lgf*TM*saf/Vf

Obj=10**3*abs(gf)+10*abs(Vf-Vobj)+abs(yf-(Robj-Re))

print(Obj,yf,Vf,gf,tf)


alpha=[]
for i in range(0,len(sol.y[0])):    
    alpha.extend([atan((sol.y[7][i])/sol.y[6][i]/sol.y[2][i])])
    
t=np.linspace(0,tf,1000)

plt.figure()
plt.subplot(4,1,1)
plt.plot(t,sol.y[2],'b--',LineWidth=2)
plt.legend(['v'],loc='best')
plt.ylabel('Velocity')
plt.subplot(4,1,2)
plt.plot(t,sol.y[1],'g:',LineWidth=2)
plt.legend(['y'],loc='best')
plt.ylabel('Position')
plt.subplot(4,1,3)
plt.plot(t,sol.y[3],'k.-',LineWidth=2)
plt.legend(['b'],loc='best')
plt.ylabel('Flight path angle')
plt.subplot(4,1,4)
plt.plot(t,alpha,'k.-',LineWidth=2)
plt.legend(['b'],loc='best')
plt.ylabel('angle of attach')


plt.show()

"""
sol = solve_ivp(traj, [0, tf], [R0,Gamma0,V0,P_l1,P_l2,P_l3] , events=[],t_eval=np.linspace(0,tf,1000))



Rf=sol.y[0][-1]
Vf=sol.y[1][-1]
Gammaf=sol.y[2][-1]


l1f=sol.y[3][-1]
l2f=sol.y[4][-1]
l3f=sol.y[5][-1]

alphaf=atan(l2f/(Vf*l3f))


Hf=l1f*(Vf*sin(Gammaf))+l2f*(TM*sin(alphaf)/Vf+(Vf**2*Rf-Ue)/(Rf**2*Vf)*cos(Gammaf))+l3f*(TM*cos(alphaf)-Ue/Rf**2*sin(Gammaf))


Obj=10*abs(Rf-Robj)+abs(Vf-Vobj)+abs(Gammaf-Gammaobj)+abs(Hf+1)

print(Obj,Rf-Re,Vf,Gammaf)

alpha=[]
for i in range(0,len(sol.y[0])):    
    alpha.extend([atan((sol.y[4][i])/sol.y[5][i]/sol.y[2][i])])


t=np.linspace(0,tf,1000)

plt.figure()
plt.subplot(4,1,1)
plt.plot(t,sol.y[2],'b--',LineWidth=2)
plt.legend(['v'],loc='best')
plt.ylabel('Velocity')
plt.subplot(4,1,2)
plt.plot(t,sol.y[0],'g:',LineWidth=2)
plt.legend(['y'],loc='best')
plt.ylabel('Position')
plt.subplot(4,1,3)
plt.plot(t,sol.y[1],'k.-',LineWidth=2)
plt.legend(['b'],loc='best')
plt.ylabel('Flight path angle')
plt.subplot(4,1,4)
plt.plot(t,alpha,'k.-',LineWidth=2)
plt.legend(['b'],loc='best')
plt.ylabel('angle of attach')


plt.figure()
plt.subplot(3,1,1)
plt.plot(t,sol.y[3],'b--',LineWidth=2)
plt.legend(['l1'],loc='best')
plt.ylabel('Velocity')
plt.subplot(3,1,2)
plt.plot(t,sol.y[4],'g:',LineWidth=2)
plt.legend(['l2'],loc='best')
plt.ylabel('Position')
plt.subplot(3,1,3)
plt.plot(t,sol.y[5],'k.-',LineWidth=2)
plt.legend(['l3'],loc='best')
plt.ylabel('Flight path angle')

plt.show()

#r,gamma,V,l1,l2,l3,
"""
