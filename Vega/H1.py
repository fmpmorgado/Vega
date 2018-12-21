from scipy.integrate import solve_ivp
import pygmo as pg
from math import sin,cos,atan,pi
import numpy as np
import matplotlib.pyplot as plt



TM=20
Ue=398604.3*10**9
    
g0=9.810665
Re=6371000

R0=Re
V0=30
Gamma0=0.5

Robj=800000+Re
Gammaobj=0.3
Vobj=7100

g0=9.810665

V=2

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
    if y[2]==0: y[2]=0.0000000000000001
    gamma=atan(y[3]/y[2])

    x=V*cos(gamma)
    y=V*sin(gamma)

    l1=0
    l2=0
    return[x,y,l1,l2]

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
    P_l1=x[0]
    P_l2=x[1]
    tf=x[2]
    

    sol = solve_ivp(traj, [0, tf], [0,0,P_l1,P_l2] , events=[],t_eval=np.linspace(0,tf,1000))



    xf=sol.y[0][-1]
    yf=sol.y[1][-1]

    l1f=sol.y[2][-1]
    l2f=sol.y[3][-1]

    gammaf=atan(l2f/l1f)

    Hf=l1f*V*cos(gammaf)+l2f*V*sin(gammaf)

    Obj=abs(1-xf)+abs(1-yf)+abs(1+Hf)

    print(Obj)
    return Obj



class Optimization:
    def fitness(self, x):
        fit=model(x)
        return [fit]

    def get_bounds(self):
        return ([-1,-1,0.1],[1,1,10])
        




prob = pg.problem(Optimization())
algo = pg.algorithm(pg.pso_gen(gen=100))          #Choose of heuristic algorithm and number of generations
pop = pg.population(prob,100)                        #Choose number of individuals
pop = algo.evolve(pop)                              #Evolve the population
y=pop.champion_x                                    #Extract the best DV division solution to minimize mass




P_l1=y[0]
P_l2=y[1]
tf=y[2]

    
sol = solve_ivp(traj, [0, tf], [0,0,P_l1,P_l2] , events=[],t_eval=np.linspace(0,tf,1000))



xf=sol.y[0][-1]
yf=sol.y[1][-1]


l1f=sol.y[2][-1]
l2f=sol.y[3][-1]

gammaf=atan(l2f/l1f)

Hf=l1f*V*cos(gammaf)+l2f*V*sin(gammaf)

Obj=abs(1-xf)+abs(1-yf)+abs(Hf+1)


print(y, Obj)


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
