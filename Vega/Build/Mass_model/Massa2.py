import math
    
    

g0=9.80665
K_RL=1


######################Solid propulsion system ##########################


#engine

def Motor(M_Prop,Thrust):
    M_M=0.677165543*10**-3*Thrust+0.076377593*M_Prop
    M_M=0.91*M_M
    M_M=0.0713984*M_Prop+0.31895*10**-3*Thrust+92.68

    #M_M=81.7+0.081509*M_Prop
    #M_M=0.08424*M_Prop+95.56
    return M_M


#########sizing

def Tank_solid_length(rho,diameter,Mass):
    Vtank=Mass/rho
    Length=(Vtank)/((diameter/2)**2*math.pi-(diameter*0.5/2)**2*math.pi)
    return Length*1.2



def Tank_solid_length(rho,diameter,Mass):
    Vtank=0.000693062*Mass+ 1.805
    #Vtank=0.000778164*Mass
    Length=(Vtank)/((diameter/2)**2*math.pi)
    return Length*1.275
######################Liquid Propulsion system #######################



#Thrust chamber

def Thrust_chamber(T):
    M_TC=T/(g0*(25.2*math.log10(T)-80.7))
    return M_TC

def Thrust_chamber_length(T):
    L_TC=(3.042*T*10**-5+327.7)*10**-3
    return L_TC

#Support structure

def Support_structure(T):
    M_ST=0.88*10**-3*(0.225*T)**1.0687
    return M_ST

#Tank

def Fuel_tank(rho,Mass,Ratio,D,Fmat,rhomat):

    Dtank=0.9*D
    Fuel_Mass=Mass/(Ratio+1)   
    Vtank=Fuel_Mass/rho
    
    Pb=2*1.2*10**(-0.10688*(math.log(Vtank)-0.2588))*10**6


    
    tc=0.5*Pb*Dtank/Fmat
    ts=0.25*Pb*Dtank/Fmat

    lc=(Vtank-4/3*math.pi*(Dtank/2)**3)/(math.pi*(Dtank/2)**2)

    if lc<0:
        rtank=(Vtank/4*3/math.pi)**(1/3)
        TankMass=(ts*4*math.pi*(rtank)**2)*rhomat
        return TankMass

    TankMass=(tc*lc*2*math.pi*Dtank/2+ts*4*math.pi*(Dtank/2)**2)*rhomat

    return TankMass   

def Oxidizer_tank(rho,Mass,Ratio,D,Fmat,rhomat):


    Dtank=0.9*D
    Fuel_Mass=Mass*Ratio/(Ratio+1)   
    Vtank=Fuel_Mass/rho

    if Vtank<=0: Vtank=100000
    Pb=2*1.2*10**(-0.10688*(math.log(Vtank)-0.2588))*10**6


    tc=0.5*Pb*Dtank/Fmat
    ts=0.25*Pb*Dtank/Fmat

    lc=(Vtank-4/3*math.pi*(Dtank/2)**3)/(math.pi*(Dtank/2)**2)

    if lc<0:
        rtank=(Vtank/4*3/math.pi)**(1/3)
        TankMass=(ts*4*math.pi*(rtank)**2)*rhomat
        return TankMass


    TankMass=(tc*lc*2*math.pi*Dtank/2+ts*4*math.pi*(Dtank/2)**2)*rhomat

    #M_Insulation=(lc*2*math.pi*Dtank/2+4*math.pi*(Dtank/2)**2)*1.123
    return TankMass   

def Tank_oxidizer_length(rho,Mass,Ratio,D):
    Dtank=1.00*D
    Fuel_Mass=Mass*Ratio/(Ratio+1)   
    Vtank=Fuel_Mass/rho

    lc=(Vtank-4/3*math.pi*(Dtank/2)**3)/(math.pi*(Dtank/2)**2)
    length=lc+Dtank
    return length
    
def Tank_fuel_length(rho,Mass,Ratio,D):
    Dtank=1.00*D
    Fuel_Mass=Mass/(Ratio+1)   
    Vtank=Fuel_Mass/rho

    lc=(Vtank-4/3*math.pi*(Dtank/2)**3)/(math.pi*(Dtank/2)**2)
    length=lc+Dtank
    return length




#####################Avionics and Power System #####################

#Avionics Akin

def Massa_avionics(M_Stage):
    M_Avionics=10*math.pow(M_Stage,0.361)
    return M_Avionics

#Avionics

def Avionics(S_tot):
    M_A=0.25*K_RL*(246.76+1.3183*S_tot)
    return M_A

#Power system

def Power_system(M_A):
    M_PS=0.82*K_RL*0.405*M_A
    return M_PS











####################### Other components ############################


def Common_formula(k1,k2,K_SM,D,S):
    M_OC=K_SM*K1*S*(D*3.2808)**k2
    return M_OC


#Intertanks Lower stage
    
def Intertanks_lower(D,S,K_SM):
    k1 = 5.4015
    k2 = 0.5169
    M_ITL=Common_formula(k1,k2,K_SM,D,S)
    return M_ITL

#Intertanks Upper stage and boosters

def Intertanks_upper(D,S,K_SM):
    k1 = 3.8664
    k2 = 0.6025
    M_ITU=Common_formula(k1,k2,K_SM,D,S)
    return M_ITU

# Interstage Lower stage

def Interstage_lower(D,S,K_SM):
    k1=7.7165
    k2=0.4856
    M_ISL=Common_formula(k1,k2,K_SM,D,S)
    return M_ISL

#Interstage Upper stage

def Interstage_upper(D,S,K_SM):
    k1=5.5234
    k2=0.5210
    M_ISU=Common_formula(k1,k2,K_SM,D,S)
    return M_ISU

#Pad interface

def Pad_interface(D,S,K_SM):
    k1=25.763
    k2=0.5498
    M_PI=Common_formula(k1,k2,K_SM,D,S)
    return M_PI




















####################################   Payload adapter
def payload_adapter_mass(payload):
    mass=0.0477536*payload**1.01317
    return mass


####################################   Payload Fairing
def payload_fairing_mass(diameter, length):
    mass= 15.4717+30.9005*diameter*length
    return mass

def payload_fairing_length(diameter, payload):
    length = 0.566993+4.12*10**-4*payload+1.733*diameter
    return length



####################################   Stage mass

def stage_mass(stage,rocket,payload,th,Ftu,rhomat):
    Total=0

    #Total=Total+Motor_case(M_prop)
    if(stage.stage_type=="solid"):
        if rocket.num_boosters>0 and stage.stage_num==1:
            Total=Total+Motor(stage.propellant_mass/rocket.num_boosters,stage.thrust/rocket.num_boosters)*rocket.num_boosters
            th=0.015
            Total=Total+math.pi*stage.diameter*stage.length*th*rhomat
        else:
            Total=Total+Motor(stage.propellant_mass,stage.thrust)
    if(stage.stage_type=="liquid"):
        Total=Total+Fuel_tank(stage.fuel_density,stage.propellant_mass,stage.OF_ratio,stage.diameter,Ftu,rhomat)
        Total=Total+Oxidizer_tank(stage.oxidizer_density,stage.propellant_mass,stage.OF_ratio,stage.diameter,Ftu,rhomat)
        Total=Total+Thrust_chamber(stage.thrust/stage.num_engines)*stage.num_engines
        Total=Total+Support_structure(stage.thrust/stage.num_engines)*stage.num_engines

        print(Fuel_tank(stage.fuel_density,stage.propellant_mass,stage.OF_ratio,stage.diameter,Ftu,rhomat),
              Oxidizer_tank(stage.oxidizer_density,stage.propellant_mass,stage.OF_ratio,stage.diameter,Ftu,rhomat),
              Thrust_chamber(stage.thrust/stage.num_engines)*stage.num_engines,
              Support_structure(stage.thrust/stage.num_engines)*stage.num_engines,
              math.pi*stage.diameter*stage.length*th*rhomat)  

    #Out structure
    stage.length=stage_length(stage)

    if(stage.stage_type=="liquid"): 
        Total=Total+math.pi*stage.diameter*stage.length*th*rhomat


   
    #Avionics and Power system
    #Total=Total+Avionics(stage.length*(stage.diameter/2)**2*math.pi)
    #Total=Total+Power_system(Massa_avionics(stage.total_mass))

    return Total

####################################   Stage sizing

def stage_length(stage):
    Length=0
    
    if(stage.stage_type=="liquid"):
        #if stage.stage_num!=1:
        Length=Length + Tank_fuel_length(stage.fuel_density,stage.propellant_mass,stage.OF_ratio,stage.diameter)
        Length=Length + Tank_oxidizer_length(stage.oxidizer_density,stage.propellant_mass,stage.OF_ratio,stage.diameter)
        Length=Length + Thrust_chamber_length(stage.thrust)
        Length=1.00*Length




    if (stage.stage_type=="solid"):
        Length=Length + Tank_solid_length(stage.oxidizer_density,stage.diameter,stage.propellant_mass)
    



    return Length    
