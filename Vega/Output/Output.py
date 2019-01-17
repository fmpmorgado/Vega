import math
import matplotlib.pyplot as plt

def Graphs(Data):
    #Data plotting

    plt.figure()
    plt.subplot(4,1,1)
    plt.plot(Data[0],Data[6],'r-',LineWidth=2)
    plt.ylabel('alpha')
    plt.subplot(4,1,2)
    plt.plot(Data[0],Data[7],'r-',LineWidth=2)
    plt.ylabel('vG')
    plt.subplot(4,1,3)
    plt.plot(Data[0],Data[8],'r-',LineWidth=2)
    plt.ylabel('vD')
    plt.subplot(4,1,4)
    plt.plot(Data[0],Data[9],'r-',LineWidth=2)
    plt.ylabel('vA')
    plt.xlabel('Time')



    plt.figure()
    plt.subplot(4,1,1)
    plt.plot(Data[0],Data[3],'b--',LineWidth=2)
    plt.legend(['v'],loc='best')
    plt.ylabel('Velocity')
    plt.subplot(4,1,2)
    plt.plot(Data[0],Data[2],'g:',LineWidth=2)
    plt.legend(['y'],loc='best')
    plt.ylabel('Position')
    plt.subplot(4,1,3)
    plt.plot(Data[0],Data[4],'k.-',LineWidth=2)
    plt.legend(['b'],loc='best')
    plt.ylabel('Flight path angle')
    plt.subplot(4,1,4)
    plt.plot(Data[0],Data[5],'k.-',LineWidth=2)
    plt.legend(['b'],loc='best')
    plt.ylabel('Mass')

    plt.xlabel('Time')
    plt.show()

    return

def Information(mission,rocket,stage):
    
        ########################## Design Output #################################

    print("**************************************************************************************")
    print("**************************************************************************************")
    print("*                                  Mission Profile                                   *")
    print("**************************************************************************************")
    print("**************************************************************************************")

    print("*****                           Altitude = {0:0.2f} km ".format(mission.final_altitude*10**-3))
    print("*****                     Final Velocity = {0:0.2f} m/s ".format(mission.final_velocity))
    print("*****                            Payload = {0:0.1f} kg ".format(mission.payload))

    print("**************************************************************************************")
    print("**************************************************************************************")

    print("\n")
    print("\n")
    print("\n")

    print("**************************************************************************************")
    print("**************************************************************************************")
    print("*                                Rocket Specifications                               *")
    print("**************************************************************************************")
    print("**************************************************************************************")

    print("*****                       Stage number = {0:0d} ".format(rocket.num_stages))
    print("*****                     Booster number = {0:0d} ".format(rocket.num_boosters))
    print("*****     Total mass (Excluding payload) = {0:0.2f} kg ".format(rocket.mass))

    print("**************************************************************************************")
    print("**************************************************************************************")

    print("\n")
    print("\n")
    print("\n")

    print("**************************************************************************************")
    print("**************************************************************************************")
    print("*                                Stages Specifications                               *")



    for x in range(0,len(stage)):
        print("**************************************************************************************")
        print("**************************************************************************************")   
        print("*****                                 Stage {0:0d}                                    *****".format(stage[x].stage_num))
        print("*****                                                                            *****")
        print("*****                                DV = {0:0.0f}     ".format(stage[x].DV))
        print("*****                          Payload  =  {0:0.2f}    ".format(stage[x].payload))
        print("*****                                                                            *****")
        print("*****                                 Sizing                                     *****")
        print("*****                            Length = {0:0.2f}".format(stage[x].length))
        print("*****                          Diameter = {0:0.2f}".format(stage[x].diameter))
        print("*****                            Volume = {0:0.2f}".format(stage[x].diameter*math.pi*stage[x].length))
        print("*****                                                                            *****")
        print("*****                                  Mass                                      *****")
        print("*****                   Propellant Mass = {0:0.2f}".format(stage[x].propellant_mass))
        print("*****                   Structural Mass = {0:0.2f}".format(stage[x].structural_mass))
        print("*****                        Total Mass = {0:0.2f}".format(stage[x].total_mass))
        print("*****         Total Mass (With Payload) = {0:0.2f}".format(stage[x].total_mass+stage[x].payload))
        print("*****                                                                            *****")
        print("*****                                Propulsion                                  *****")
        print("*****                              Type = ",stage[x].stage_type)
        print("*****                               Isp = {0:0d} s".format(stage[x].Isp))
        print("*****                            Thrust = {0:0.0f} kN".format(stage[x].thrust*10**-3))
        print("*****                                                                            *****")
        print("*****                                Propellant                                  *****")    
        print("*****                              Name =                                    ") #Após criar dicionário de Fuel
        if(stage[x].stage_type=="liquid"):
            print("*****                         O/F ratio = {0:0.2f}".format(stage[x].OF_ratio))
            print("*****                  Oxidizer density = {0:0d} kg/m^3".format(stage[x].oxidizer_density))
            print("*****                      fuel density = {0:0d} kg/m^3".format(stage[x].fuel_density))
        elif(stage[x].stage_type=="solid"):
            print("*****                Propellant density = {0:0d} kg/m^3".format(stage[x].oxidizer_density)) #Depois mudar
            

    print("**************************************************************************************")
    print("**************************************************************************************")     

    return

