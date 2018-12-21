import pygmo as pg
import numpy as np
import math
import Massa2

################ Constants and variables initialization ################

#Gravity constants
g0=9.810665



#stage wall thickness considered
th=0.01

#Material constants --- Criar dicionário depois
Ftu=310*10**6
rhomat=2700 #al alloy

###### Mudar depois para Config de material




def Rocket_build(mission,rocket,stage,general):
#def Rocket_build(mission,rocket,stage,TWR_active=0, optimize=0):


    DV_lb=[0.2,0.2,0.2,0.05]             #DV Boundaries for the PSO algorithm
    DV_ub=[0.6,0.6,0.6,0.2]



#Model used to find the fitness in the PSO algorithm

    def model(y):
        DV_sum=0

        rocket.mass=mission.payload+rocket.fairing

        for x in range(len(stage)-1,-1,-1):
            stage[x].payload=rocket.mass
            stage[x].DV_Ratio=y[x]
            DV_sum=DV_sum+y[x]

            while True:
                aux=stage[x].structural_factor
                
                Kn=math.exp(stage[x].DV/(g0*stage[x].Isp))
                stage[x].propellant_mass = (Kn-1)*(stage[x].payload+stage[x].structural_mass)
                stage[x].structural_mass = Massa2.stage_mass(stage[x],rocket.payload,th,Ftu,rhomat)
                stage[x].structural_factor = stage[x].structural_mass/(stage[x].structural_mass+stage[x].propellant_mass)

                stage[x].total_mass=stage[x].structural_mass+stage[x].propellant_mass+stage[x].payload
                if general.TWR:
                    stage[x].thrust=stage[x].total_mass*g0*stage[x].TWR

                if (abs(stage[x].structural_factor-aux)/stage[x].structural_factor)<10**-6:
                    stage[x].total_mass=stage[x].structural_mass+stage[x].propellant_mass
                    break
                
            rocket.mass=rocket.mass+stage[x].total_mass

        return 10*rocket.mass+1e6*abs(DV_sum-1)




    class DV_Optimization:
        def fitness(self, x):
            return [model(x)]

        def get_bounds(self):
            return (DV_lb,DV_ub)



    """    
    if general.Design_optimization:
        prob = pg.problem(DV_Optimization())
        algo = pg.algorithm(pg.pso_gen(gen = 300))          #Choose of heuristic algorithm and number of generations
        pop = pg.population(prob,300)                     #Choose number of individuals
        pop = algo.evolve(pop)                              #Evolve the population
        y=pop.champion_x                                    #Extract the best DV division solution to minimize mass


        
        for x in range(len(stage)-1,-1,-1):
            stage[x].DV_Ratio=y[x]                              #Replace the DV in the stages
            stage[x].DV=stage[x].DV_Ratio*rocket.DV

    """


    rocket.mass=mission.payload+rocket.fairing
    


                
    for x in range(rocket.num_stages-1,-1,-1):              #Building the rocket, taking in consideration the mass equations in the module Mass 
        stage[x].payload=rocket.mass                        #The function starts Building the rocket beginning with the last stage. Each stage iteration
                                                            #ends when the structural factor converges
        while True:
            aux=stage[x].structural_factor
            stage[x].DV=stage[x].DV_Ratio*rocket.DV
            
            Kn=math.exp(stage[x].DV/(g0*stage[x].Isp))
            stage[x].propellant_mass = (Kn-1)*(stage[x].payload+stage[x].structural_mass)
            stage[x].structural_mass = Massa2.stage_mass(stage[x],rocket.payload,th,Ftu,rhomat)
            stage[x].structural_factor = stage[x].structural_mass/(stage[x].structural_mass+stage[x].propellant_mass)

            stage[x].total_mass=stage[x].structural_mass+stage[x].propellant_mass+stage[x].payload
            if general.TWR:
                stage[x].thrust=stage[x].total_mass*g0*stage[x].TWR

            if (abs(stage[x].structural_factor-aux)/stage[x].structural_factor)<10**-6:
                stage[x].total_mass=stage[x].structural_mass+stage[x].propellant_mass
                break



        
        rocket.mass=rocket.mass+stage[x].total_mass


    #for x in range(rocket.num_stages-1,-1,-1):
    #   stage[x].DV=g0*stage[x].Isp*math.log((stage[x].payload+stage[x].total_mass)/(stage[x].payload+stage[x].structural_mass))


    return rocket,stage




