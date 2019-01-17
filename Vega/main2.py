import math
import sys


sys.path.append('./Configuration')
from Config_B import Rocket_initialize

sys.path.append('./Output')
from Output import Graphs,Information


import matplotlib.pyplot as plt

import time
import pygmo as pg

sys.path.append('./Build')
sys.path.append('./Build/Mass_model')
from Build import Rocket_build


sys.path.append('./Trajectory')
sys.path.append('./Trajectory/Aero')

from Trajectory import init_trajectory

g0=9.810665

optim=True

now = time.clock()

from Atmosphere import atmosphere
from Aerodynamics import Knudsen,Drag_force,Drag_Coeff,Dynamic_pressure
from math import exp,pi

mission,rocket,stage,trajectory,general=Rocket_initialize()

#print(rocket.DV)
Information(mission,rocket,stage)




general.TWR=False



old=0
new=rocket.DV    
Obj=10000

optim=True

rocket.gen=0
rocket.optim=True
rocket.optim2=True
rocket.count=0

rocket,stage=Rocket_build(mission,rocket,stage, general)


print(rocket.DV)

Information(mission,rocket,stage)



trajectory.alpha=[0,0]
trajectory.l=[0,0,0]


optim = True


trajectory.vertical_time=trajectory.vertical_time+trajectory.pitch_time



def traj_model(x,rocket):
    trajectory.alpha[0]=x[0]
    trajectory.alpha[1]=x[1]
    trajectory.coast_time=x[2]
    trajectory.pitch=x[3]
    trajectory.l[0]=x[4]
    trajectory.l[1]=x[5]
    trajectory.l[2]=x[6]
    trajectory.aux=x[7]

    if rocket.optim==True:

        Data,new,Obj=init_trajectory(mission,rocket,stage,trajectory,general)
        rocket.gen=rocket.gen+1/100
        print(Data[2][-1],Data[3][-1],Data[4][-1])
    else: return 1e20



    if abs(Data[2][-1]-mission.final_altitude)>10 or abs(Data[3][-1]-mission.final_velocity)>5 or abs(Data[4][-1])>0.01:
        return Obj
    else:
         rocket.optim=False
         return 0
    
    





def traj_model2(x,rocket):
    trajectory.alpha[0]=x[0]
    trajectory.alpha[1]=x[1]
    trajectory.coast_time=x[2]
    trajectory.pitch=x[3]
    trajectory.l[0]=x[4]
    trajectory.l[1]=x[5]
    trajectory.l[2]=x[6]
    trajectory.aux=x[7]

    if rocket.optim2==False:
        return 0

    Data,new,Obj=init_trajectory(mission,rocket,stage,trajectory,general)
    #print(Data[2][-1],Data[3][-1],Data[4][-1],Data[0][-1])
        

    
    if abs(Data[2][-1]-mission.final_altitude)>10000 or abs(Data[3][-1]-mission.final_velocity)>100 or abs(Data[4][-1])>0.15:
        #print("A")
        return 1e20

    else:
        #print("B")
        print(Data[2][-1],Data[3][-1],Data[4][-1])

        rocket.optim2=False
        return 0




def design_model(x,rocket,stage):
    DV_sum=0
            
    for i in range(0,len(stage)):
        #stage[i].TWR=x[i]
        stage[i].DV_Ratio=x[i]
        DV_sum=DV_sum+stage[i].DV_Ratio
    
    rocket,stage=Rocket_build(mission,rocket,stage, general)
    #Information(mission,rocket,stage)
    #print(rocket.mass)
    rocket.count=rocket.count+1

    prob = pg.problem(Traj_Optimization2())
    algo = pg.algorithm(pg.pso_gen(gen=50))          #Choose of heuristic algorithm and number of generations
    pop = pg.population(prob,100)                        #Choose number of individuals
    pop = algo.evolve(pop)                              #Evolve the population
    y=pop.champion_x                                    #Extract the best DV division solution to minimize mass

    trajectory.alpha[0]=y[0]
    trajectory.alpha[1]=y[1]
    trajectory.coast_time=y[2]
    trajectory.pitch=y[3]
    trajectory.l[0]=y[4]
    trajectory.l[1]=y[5]
    trajectory.l[2]=y[6]
    trajectory.aux=y[7]



    Data,new,Obj=init_trajectory(mission,rocket,stage,trajectory,general)


    if general.constraints:
        for i in range(0,len(Data[0])):
            if Data[12][i]/Data[5][i]/g0>5 or Dynamic_pressure(1.225*exp(-Data[2][i]/8440),Data[3][i])>50000: return 1e20

    rocket.optim2=True
        
    Obj=pop.champion_f

    #print(y, Obj)
    
    if Obj==0:
        print(rocket.count, rocket.mass, "Worked")
        print(stage[0].TWR,stage[1].TWR,stage[2].TWR,stage[0].DV_Ratio,stage[1].DV_Ratio,stage[2].DV_Ratio)
        return rocket.mass

    else:
        print(rocket.count, rocket.mass)
        return rocket.mass*10



class Traj_Optimization:
    def fitness(self, x):
        fit=traj_model(x,rocket)
        return [fit]

    def get_bounds(self):
        return ([-0.8,-0.8,1,1.47,-0.999,-0.999,-0.999,0.01],[0,0,350,1.56,0.999999,0.999999,0.999999,0.999])


class Traj_Optimization2:
    def fitness(self, x):
        fit=traj_model2(x,rocket)
        return [fit]

    def get_bounds(self):
        return ([-0.8,-0.8,250,1.45,-0.999,-0.999,-0.999,0.50],[0,0,350,1.56,0.999999,0.999999,0.999999,0.999])

        

class Design_Optimization:
    def fitness(self, x):
        fit=design_model(x,rocket,stage)
        return [fit]

    def get_bounds(self):
      #  return ([stage[0].TWR*0.75,stage[1].TWR*0.75,stage[2].TWR*0.75,stage[3].TWR*0.75,stage[0].DV_Ratio*0.9,stage[1].DV_Ratio*0.9,stage[2].DV_Ratio*0.9,stage[3].DV_Ratio*0.9],
    #         [stage[0].TWR*1.25,stage[1].TWR*1.25,stage[2].TWR*1.25,stage[3].TWR*1.25,stage[0].DV_Ratio*1.05,stage[1].DV_Ratio*1.05,stage[2].DV_Ratio*1.05,stage[3].DV_Ratio*1.05])
        return([stage[0].DV_Ratio*0.9,stage[1].DV_Ratio*0.9,stage[2].DV_Ratio*0.9],[stage[0].DV_Ratio*1.05,stage[1].DV_Ratio*1.05,stage[2].DV_Ratio*1.05])        




if general.Design_optimization:
    prob = pg.problem(Design_Optimization())
    algo = pg.algorithm(pg.sga(gen=1))          #Choose of heuristic algorithm and number of generations
    pop = pg.population(prob,125)                        #Choose number of individuals
    pop = algo.evolve(pop)                              #Evolve the population
    y=pop.champion_x                               #Extract the best DV division solution to minimize mass


    #print(y)

    DV_sum=0
    for i in range(0,len(stage)):
    #    stage[i].TWR=y[i]
        stage[i].DV_Ratio=y[i]
        DV_sum=DV_sum+stage[i].DV_Ratio

    rocket,stage=Rocket_build(mission,rocket,stage, general)
    Information(mission,rocket,stage)
    print(y)

    
    print("Time=", time.clock()-now)



"""

stage[0].TWR=1.4696
stage[1].TWR=3.3822
stage[2].TWR=1.933
stage[3].TWR=0.0853

stage[0].DV_Ratio=0.2244
stage[1].DV_Ratio=0.20514
stage[2].DV_Ratio=0.3714
stage[3].DV_Ratio=0.058


rocket,stage=Rocket_build(mission,rocket,stage, general)
Information(mission,rocket,stage)

"""


if rocket.num_boosters !=0:
    time_burn=stage[0].propellant_mass*stage[0].Isp*g0/stage[0].thrust
    stage[0].propellant_mass=stage[0].propellant_mass+stage[1].thrust/stage[1].Isp/g0*time_burn
    stage[1].propellant_mass=stage[1].propellant_mass-stage[1].thrust/stage[1].Isp/g0*time_burn
    stage[0].payload=stage[0].payload-stage[1].thrust/stage[1].Isp/g0*time_burn

if general.Trajectory:
    prob = pg.problem(Traj_Optimization())
    algo = pg.algorithm(pg.pso_gen(gen=500))          #Choose of heuristic algorithm and number of generations
    pop = pg.population(prob,250)                        #Choose number of individuals
    pop = algo.evolve(pop)                              #Evolve the population
    y=pop.champion_x                                    #Extract the best DV division solution to minimize mass


    trajectory.alpha[0]=y[0]
    trajectory.alpha[1]=y[1]
    trajectory.coast_time=y[2]
    trajectory.pitch=y[3]
    trajectory.l[0]=y[4]
    trajectory.l[1]=y[5]
    trajectory.l[2]=y[6]
    trajectory.aux=y[7]


    print()

"""
trajectory.coast_time=2.77594163e+02
trajectory.pitch=1.47575848e+00
trajectory.l[0]=-2.10763106e-03
trajectory.l[1]=7.12057203e-01 
trajectory.l[2]=-8.85632384e-01
trajectory.aux=9.37180755e-01
"""



Data,new,Obj=init_trajectory(mission,rocket,stage,trajectory,general)
print(Data[2][-1], time.clock()-now)

    
#print("Time=", time.clock()-now, y)

kn=[]
D=[]
a=[]
q=[]
cd=[]

"""
for i in range(0,len(Data[0])):
    T=atmosphere(Data[2][i])[0]
    rho=1.225*exp(-Data[2][i]/8440)
    L=0.1
    V=Data[3][i]
    Sref=4
         
    kn.extend([Knudsen(L,rho,T)])
    D.extend([1/2*1.225*exp(-Data[2][i]/8440)*Sref*Drag_Coeff(V/math.sqrt(1.4*287.053*T))*V**2])
    a.extend([Data[12][i]/Data[5][i]/g0])
    cd.extend([Drag_Coeff(V/math.sqrt(1.4*287.053*T))])    
    q.extend([Dynamic_pressure(rho,V)])
    

plt.plot(Data[2][0:650], kn[0:650])
plt.show()

plt.plot(Data[0],a)
plt.show()

plt.plot(Data[0][0:650],D[0:650])
plt.plot(Data[0][0:650],D2[0:650])

plt.show()
"""
    

Graphs(Data)


