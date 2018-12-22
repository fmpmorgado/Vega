import math
import numpy as np
from Classes import Trajectory, Rocket, Mission, Stage, General, Optimization

######################Propellant Dictionary #########################
#Mudar depois

Propellant=[{'LOX/H2':('Liquid',462,3.8,1142,71)},
            {'LOX/N2H4':('Liquid',363,1.2,1142,1010)},
            {'LOX/RP1':('Liquid',347,2.27,1142,810)},
            {'N2O4/RP1':('Liquid',328,3.51,1440,810)},
            {'N2O4/HTPB':('Liquid',297,3.17,1440,1810)},
            {'LOX/HTPB':('Liquid',317,2.04,1142,1810)},
            {'F2/H2':('Liquid',441,4.26,1509,71)},
            {'N2O4/MMH':('Liquid',318,1,1440,878)},
            {'N2O4/UDMH':('Liquid',313,2.21,1440,1010)},
            {'N2O4/N2H4':('Liquid',309,2.21,1440,1010)},
            {'HTPB/AP':('Solid',265,1880)},
            {'DB-H;X/AP':('Solid',260,1854)}]


            



#Universal constants used in this module

g0=9.810665  #Earth gravity at sea level



########################Auxiliary functions ##############################
# This functions are to find the payload and the stages DV, using the stages mass.
# The functions are only used in the case the user wants to use a known rocket and doesn't
# want to deduce the payload and DV of each stage by hand


def DV_budget(stage,total_mass,num_stages):
    payload=total_mass
    DV=[0]*num_stages

    for x in range(0,num_stages):
        payload=payload-stage[x].total_mass

    for x in range(num_stages-1,-1,-1):
        DV[x]=g0*stage[x].Isp*math.log((payload+stage[x].total_mass)/(payload+stage[x].structural_mass))
        payload=payload+stage[x].total_mass
    return DV
    
def Payload(stage,total_mass,num_stages):
    payload=total_mass
    
    for x in range(0,num_stages):
        payload=payload-stage[x].total_mass

    return payload





def Rocket_initialize():
    
# (Comparison) - In the case of trying to use a pre existing rocket as input, to compare the program output with the rocket
# (*) - Not mandatory to fill, only in the case the user wants to use a pre existing rocket to simulate a trajectory
# (Mandatory) - Need to be filled

######################################### Setup ########################################



    ########### General Initialization #####################

    general= General()

    general.Trajectory=True                 #In order to perform a trajectory simulation    
    general.New_design=False             #If the user wants to use a pre existing rocket              
    general.TWR=True             #If the user wants to use Thrust to weight ratio for each iteration, instead of defining the initial Thrust
    general.Design_optimization=True #Design and DV optimization

    #############  Trajectory Initialization#####################################

    trajectory=Trajectory() 

    trajectory.delta_pitch=0.1    


    trajectory.pitch_time=8    #Pitch maneuver duration 
    trajectory.pitch=1.50        #Final flight path angle of Pitch maneuver

    
    trajectory.vertical_altitude=50    #Altitude to end vertical maneuver and start pitch
    trajectory.vertical_time=3         #Otherwise, the user can specify the time to end the maneuver
                                        #Specifying both


    trajectory.coast=True             #Coast in the last stage
    trajectory.coast_time=50            #Coast time

    #############  Mission setup and initialization #########################
        
    mission=Mission()

    mission.initial_altitude=0          #Initial Rocket Altitude
    mission.initial_velocity=0.1          #Initial Rocket Velocity
    mission.initial_flight_angle=np.pi/2 #Initial Flight Path Angle


    mission.final_altitude=700000       #Final Rocket Altitude
    mission.final_velocity=7500       #Final Rocket Velocity
    mission.final_flight_angle=0.0        #Final Flight Path Angle
    
    mission.payload=1430+900                 #Mission Payload to carry


    #############  Rocket setup and initialization ##################

    rocket=Rocket()

    rocket.num_stages=3               #Number of Rocket stages
    rocket.num_boosters=0               #Number of Rocket boosters

    rocket.DV=10000                     #Initial guess of Rocket DV, which has to englobe thrust vectoring, gravity and drag losses
    rocket.mass=137000                  #Initial Rocket mass (comparison only)

    rocket.fairing=550

    #############  Stage setup and initialization ######################

    #### Parameters to include in the simulation, by stage order, i.e.,
    #the array starts with the specifications of the first stage and so on


    ###     SIZING      ####
    
    Diameter=[3,1.9,1.9,1.9]      #Stage diameter (mandatory)
    length=[11.7,8.39,4.12,1.7]     #Stage length (*)



    ###     PROPULSION  ####

    #The user can choose to run the program with a Thrust specified for each stage
    #or can choose to pick a Thrust to Weight Ratio, where the Thrust will depend
    # on the Rocket mass of the current iteration
     
#    Thrust=[2261000,1196000,260000,2420]*rocket.num_stages      #Stage Thrust (depends on the user choice)
    Thrust=[2261000,1196000,260000]*rocket.num_stages      #Stage Thrust (depends on the user choice)


    TWR=[2]*rocket.num_stages                       #Thrust to Weight ratio (depends on the user choice)
    
    Isp=[280,289,296]                   #Specific impulse of each stage (mandatory)

    Num_engines=[1,1,1,1]             #Number of engines of each stage (mandatory)
    stage_type=["solid","solid","solid","liquid"] #Type of propellant (mandatory)

    OF_ratio=[2.61,2.61,2.61,2.61]       #Propellant Oxidizer Fuel Ratio. (Only in the presence of liquid stages)
                                    #The array needs to be of the same size as the number of stages included

    Oxidizer_density=[1800,1800,1800,1142]*rocket.num_stages       #The oxidizer density for liquid stages only
    Fuel_density=[810]*rocket.num_stages            #The fuel density for liquid stages only
    
    #Falta trocar por um dicionário

    nozzle_area=[7]*rocket.num_stages

    #######  MASS   ######
    
 #   total_mass=[88365+7431,23906+1845,10115+833,550+418]*rocket.num_stages           #Total Mass of the stage (comparison only)
 #   propellant_mass=[88365,23906,10115,550]*rocket.num_stages                        #Propellant Mass of the stage (comparison only)

    total_mass=[88365+7431,23906+1845,10115+833]           #Total Mass of the stage (comparison only)
    propellant_mass=[88365,23906,10115]                        #Propellant Mass of the stage (comparison only)



    Structural_Factor=[0.0]*rocket.num_stages                                        #Structural Factor (*)



    #### VELOCITY BUDGET ##

    DV=[2950]*rocket.num_stages                 #DV of each stage (*)
    
    DV_Ratio=[0.3]*rocket.num_stages                          #DV ratio to divide DV for each rocket stage. It can be chosen by the user
                                                    #or can be selected using a PSO algorithm to reduce the rocket mass, respecting the
                                                    #DV chosen in the rocket.DV, to fulldill the mission


    ####### Stage initialization ########

    stage = [Stage() for _ in range(rocket.num_stages)]


    for x in range(0,rocket.num_stages):
        stage[x].stage_num=x+1
        stage[x].diameter=Diameter[x]
        stage[x].thrust=Thrust[x]
        stage[x].TWR=TWR[x]
        stage[x].Isp=Isp[x]
        stage[x].num_engines=Num_engines[x]
        stage[x].OF_ratio=OF_ratio[x]
        stage[x].stage_type=stage_type[x]
        stage[x].DV=DV[x]
        
        stage[x].payload=mission.payload
        stage[x].propellant_mass=propellant_mass[x]
        stage[x].structural_mass=total_mass[x]-propellant_mass[x]
        stage[x].total_mass=total_mass[x]
        stage[x].length=length[x]
        stage[x].oxidizer_density= Oxidizer_density[x]
        stage[x].fuel_density= Fuel_density[x]
        stage[x].structural_factor=Structural_Factor[x]
        stage[x].DV_Ratio=DV_Ratio[x]

        stage[x].nozzle_area=nozzle_area[x]


    # The code below is only used if the user wants to use a pre existing rocket and wants to compare.
    
    if general.New_design==False:
        rocket.DV=0
        payload=Payload(stage,rocket.mass,rocket.num_stages)  #Payload given by inputs above
        rocket.mass=rocket.mass-payload+mission.payload+rocket.fairing    #Replacing payload and updating Rocket.mass
        DV=DV_budget(stage,rocket.mass,rocket.num_stages)
        stage[rocket.num_stages-1].payload=mission.payload

        for x in range(0,rocket.num_stages):
            stage[x].DV=DV[x]
            rocket.DV=rocket.DV+DV[x]

        for x in range(0,rocket.num_stages):
            stage[x].DV_Ratio=stage[x].DV/rocket.DV
        
        payload=mission.payload+rocket.fairing
        for x in range(rocket.num_stages-1,-1,-1):
            stage[x].payload=payload
            stage[x].TWR=stage[x].thrust/(stage[x].total_mass+stage[x].payload)/g0
            payload=payload+stage[x].total_mass

    else:
        rocket.mass=stage[0].payload+stage[0].total_mass
                
    return mission,rocket,stage,trajectory,general
