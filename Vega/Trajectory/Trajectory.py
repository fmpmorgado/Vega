import numpy as np
from gekko import GEKKO
import matplotlib.pyplot as plt
from math import sin,cos,sqrt,atan,exp,atan2
from scipy.integrate import solve_ivp
from scipy import polyfit, poly1d
from Atmosphere import interp, atmosphere,Temp_func,rho_func,pressure_func
from Aerodynamics import Drag_force,Drag_Coeff,Dynamic_pressure,Knudsen
from Vertical import Vertical_flight
from Gravity_Turn import Gravity_turn
from Free_Flight import free_flight
from Free_Flight2 import free_flight2
from Pitch import pitch
from Coast_arc import coast_arc
from const_alpha import const_alpha
from Output import Graphs,Information



def init_trajectory(mission,rocket,stage,trajectory,general):
    Traj=True

    

    
    #Constants Definitions
    g0=9.810665
    Re=6371000
    Ue=398604.3*10**9

    #Initialize variables for plotting
    Time=[0]
    x=[0]
    y=[mission.initial_altitude]
    v=[mission.initial_velocity]
    phi=[mission.initial_flight_angle]
    mass=[rocket.mass]
    alpha=[0]
    vG=[0]
    vD=[0]
    vA=[0]

    Data=[Time,x,y,v,phi,mass,alpha,vG,vD,vA,0,0,[0]]

    Hi=[]
    Hf=[]


#    print(Data)
    L=0.1

    def Hamiltonian(t,y):
#Temp_func(y[1])

        
        if y[1]>0 and y[1]<5000000:
            D=Drag_force(1.225*exp(-y[1]/8440),atmosphere(y[1])[0],Area,y[2],L)
        else:
            D=0

        ax=y[0]
        ay=y[1]
        V=y[2]
        gamma=y[3]
        m=y[4]

        lr=y[8]
        lg=y[9]
        lv=y[10]



        if x==0 and ay>=0 and ay<5000000: T=T_aux-stage[0].nozzle_area*atmosphere(ay)[1]
        else: T=T_aux


        if Hamilton==True:
            sa=-lg/(V*sqrt((lg/V)**2+lv**2))
            ca=-lv/sqrt((lg/V)**2+lv**2)
        else:
            sa=sin(alpha)
            ca=cos(alpha)


        g=Ue/(ay+Re)**2
        TM=T/m
        r=ay+Re
        
        dx=V*cos(gamma)
        dy=V*sin(gamma)
        dV=TM*ca-g*sin(gamma)-D/m
        dgamma=TM*sa/V+(V/r-g/V)*cos(gamma)
        dvg=g*sin(gamma)
        dD=D/m
        dva=TM-TM*ca

                
        #print(TM,T,m,dva,ca,sa)        



        if Hamilton==True:
            dlr=V*lg*cos(gamma)/r**2-(2*Ue*lv*sin(gamma)+2*Ue*lg*cos(gamma)/V)/r**3
            dlg=-lr*cos(gamma)*V+Ue*lv*cos(gamma)/r**2+lg*sin(gamma)*(V/r-Ue/(r**2*V))
            dlv=-lr*sin(gamma)-lg*(cos(gamma)*(1/r+Ue/((r**2)*(V**2)))-TM/V**2*sa)

        else:
            dlr=0
            dlg=0
            dlv=0
            
        
        #print(ay)
        
        return [dx,dy,dV,dgamma,dm,dvg,dD,dva,dlr,dlg,dlv]        


    if Traj==True:

        GT=True
        Pitch=True
        Hamilton=False
        Coast=True        

        

        x=0

        if rocket.num_boosters==0:

            while x!=-1:

                aux=x
                T=stage[x].thrust
                ISP=stage[x].Isp
                Area=stage[x].diameter**2/4*np.pi

                dm=-T/ISP/g0
                #print(dm)
                tb=-stage[x].propellant_mass/dm
                #Stage termination
                kn=Knudsen(L,1.225*exp(-Data[2][-1]/8440),atmosphere(Data[2][-1])[0])

                #Not knowing Tvac


                if GT==True: payload=stage[x].payload
                else: payload=stage[x].payload-rocket.fairing

                def stage_termination(t,y):

                    return payload+stage[x].structural_mass-y[4]+0.001

                stage_termination.terminal=True


                def gravity_turn_termination(t,y):
                    if y[1]<0: return 1
                    if Knudsen(L,1.225*exp(-y[1]/8440),atmosphere(y[1])[0])>10:
                        return 0
                    else:
                        return 1

                    #return y[1]-100000

                gravity_turn_termination.terminal=True

                if kn>10:  gravity_turn_termination.terminal=False


                            
                #Vertical condition
                if Data[2][-1]==0:
                    tb=trajectory.vertical_time #+trajectory.pitch_time #Mudar depois

                alpha=0

                if GT==False and x<len(stage)-1:
                    if alpha_aux+1>=2:
                        break
                    alpha= trajectory.alpha[alpha_aux];
                    alpha_aux=alpha_aux+1
         
                if x==len(stage)-1 and Coast==True: T=0; dm=0; tb=trajectory.coast_time; Coast=False


                if x==len(stage)-1 and T!=0: Hamilton=True; tb=tb*trajectory.aux


                if x==0: T_aux=T+stage[0].nozzle_area*pressure_func(0)
                else:
                    T_aux=T
                    
                #print(x,Data[2][-1],Data[5][-1])

                Hi.extend([
                    trajectory.l[0]*Data[3][-1]*sin(Data[4][-1])
                    +trajectory.l[1]*(T_aux/Data[5][-1]*sin(Data[4][-1])/Data[3][-1]
                                      +(Data[3][-1]/(Data[2][-1]+Re)-Ue/((Data[2][-1]+Re)**2*Data[3][-1]))*cos(Data[4][-1]))
                    +trajectory.l[2]*(T_aux/Data[5][-1]*cos(Data[4][-1])+(Data[3][-1]/(Data[2][-1]+Re)**2))*sin(Data[4][-1])])
                
                sol = solve_ivp(Hamiltonian, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],Data[9][-1],trajectory.l[0],trajectory.l[1],trajectory.l[2]], events=[stage_termination,gravity_turn_termination], t_eval=np.linspace(0,tb,250))





                tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
                #print(sol.y[1][-1])

                Data[0].extend(tm[1:])
                Data[1].extend(sol.y[0][1:])
                Data[2].extend(sol.y[1][1:])

                if Data[2][-1]<0: return Data,0,1e20 


                Data[3].extend(sol.y[2][1:])
                Data[4].extend(sol.y[3][1:])
                Data[5].extend(sol.y[4][1:])
                

                for i in range(1,len(sol.y[0])):
                    if x==0:
                        Data[12].extend([T_aux-stage[x].nozzle_area*atmosphere(Data[2][-len(sol.y[0])+i])[1]])
                    else:
                        Data[12].extend([T_aux])


                if Hamilton==False:
                    Data[6].extend(np.ones(len(sol.t)-1)*alpha)
                else:
                    for i in range(1,len(sol.y[0])):

                        sa=-sol.y[9][i]/(sol.y[2][i]*sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                        ca=-sol.y[10][i]/(sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                        Data[6].extend([atan2(sa,ca)])
                    
                Data[7].extend(sol.y[5][1:])
                Data[8].extend(sol.y[6][1:])
                Data[9].extend(sol.y[7][1:])

                trajectory.l[0]=sol.y[8][-1]
                trajectory.l[1]=sol.y[9][-1]
                trajectory.l[2]=sol.y[10][-1]

                
                Hf.extend(
                    [trajectory.l[0]*Data[3][-1]*sin(Data[4][-1])
                    +trajectory.l[1]*(T_aux/Data[5][-1]*sin(Data[4][-1])/Data[3][-1]
                                      +(Data[3][-1]/(Data[2][-1]+Re)-Ue/((Data[2][-1]+Re)**2*Data[3][-1]))*cos(Data[4][-1]))
                    +trajectory.l[2]*(T_aux/Data[5][-1]*cos(Data[4][-1])+(Data[3][-1]/(Data[2][-1]+Re)**2))*sin(Data[4][-1])])


                kn=Knudsen(L,1.225*exp(-Data[2][-1]/8440),atmosphere(Data[2][-1])[0])


                #print(x,Data[2][-1],Data[5][-1],Data[0][-1])
                #print()

                
                if Pitch==True: Data[4][-1]=trajectory.pitch; x=x-1; Pitch=False

                if kn>=10 and GT == True:
                    #print(Data[0][-1])
                    payload=stage[x].payload-rocket.fairing
                    Data[5][-1]=Data[5][-1]-rocket.fairing
                    x=x-1;
                    GT=False;
                    Hamilton=True;
                    alpha_aux=0;


                if sol.status!=-1 and x!=len(stage)-1: x=x+1


                if x!=aux:
                    Data[5][-1]=payload
                    if x<=len(stage)-1:
                        Data[12][-1]=stage[x].thrust

                #print(x, Data[2][-1], T, Hamilton)
                #print(x,Data[2][-1],T)
                

                #print(Hi,Hf)
                if Hamilton==True and x==len(stage)-1 and Coast==False and T!=0: break
                if Data[2][-1]<0: break


















        else:

            #Change propellant mass of the fictional stages to ease the algorithm

#            print(time_burn)

            while x!=-1:

              #  if x==0:
                

                aux=x
                T=stage[x].thrust
                ISP=stage[x].Isp
                Area=stage[x].diameter**2/4*np.pi

                dm=-T/ISP/g0
                
                if x==0:
                    T=stage[0].thrust+stage[1].thrust
                    dm=-stage[0].thrust/stage[0].Isp/g0-stage[1].thrust/stage[1].Isp/g0
                    Area=stage[0].diameter**2/4*np.pi*rocket.num_boosters+stage[1].diameter**2/4*np.pi

                
                tb=-stage[x].propellant_mass/dm
                #print(T,tb)

                                
                #Stage termination
                kn=Knudsen(L,1.225*exp(-Data[2][-1]/8440),atmosphere(Data[2][-1])[0])

                #Not knowing Tvac


                if GT==True: payload=stage[x].payload
                else: payload=stage[x].payload-rocket.fairing

                def stage_termination(t,y):
                    #print(y[4],y[1],t)
                    return payload+stage[x].structural_mass-y[4]+0.001

                stage_termination.terminal=True


                def gravity_turn_termination(t,y):
                    if y[1]<0 or y[1] > 5000000: return 1
                   # print(Knudsen(L,1.225*exp(-y[1]/8440),atmosphere(y[1])[0]),y[1])

                    if Knudsen(L,1.225*exp(-y[1]/8440),atmosphere(y[1])[0])>10:
                        return 0
                    else:
                        return 1

                    #return y[1]-100000

                gravity_turn_termination.terminal=True

                if kn>10:  gravity_turn_termination.terminal=False


                            
                #Vertical condition
                if Data[2][-1]==0:
                    tb=trajectory.vertical_time #+trajectory.pitch_time #Mudar depois

                alpha=0

                if GT==False and x<len(stage)-1:
                    if alpha_aux+1>=2:
                        break
                    alpha= trajectory.alpha[alpha_aux];
                    alpha_aux=alpha_aux+1
         
                if x==len(stage)-1 and Coast==True: T=0; dm=0; tb=trajectory.coast_time; Coast=False


                if x==len(stage)-1 and T!=0: Hamilton=True; tb=tb*trajectory.aux


                if x==0: T_aux=T+stage[0].nozzle_area*pressure_func(0)
                else:
                    T_aux=T
                    
                #print(x,Data[2][-1],Data[5][-1])

                Hi.extend([
                    trajectory.l[0]*Data[3][-1]*sin(Data[4][-1])
                    +trajectory.l[1]*(T_aux/Data[5][-1]*sin(Data[4][-1])/Data[3][-1]
                                      +(Data[3][-1]/(Data[2][-1]+Re)-Ue/((Data[2][-1]+Re)**2*Data[3][-1]))*cos(Data[4][-1]))
                    +trajectory.l[2]*(T_aux/Data[5][-1]*cos(Data[4][-1])+(Data[3][-1]/(Data[2][-1]+Re)**2))*sin(Data[4][-1])])
                
                sol = solve_ivp(Hamiltonian, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],Data[9][-1],trajectory.l[0],trajectory.l[1],trajectory.l[2]], events=[stage_termination,gravity_turn_termination], t_eval=np.linspace(0,tb,250))





                tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
                #print(sol.y[1][-1])

                Data[0].extend(tm[1:])
                Data[1].extend(sol.y[0][1:])
                Data[2].extend(sol.y[1][1:])

 


                Data[3].extend(sol.y[2][1:])
                Data[4].extend(sol.y[3][1:])
                Data[5].extend(sol.y[4][1:])
                



                if Hamilton==False:
                    Data[6].extend(np.ones(len(sol.t)-1)*alpha)
                else:
                    for i in range(1,len(sol.y[0])):

                        sa=-sol.y[9][i]/(sol.y[2][i]*sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                        ca=-sol.y[10][i]/(sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                        Data[6].extend([atan2(sa,ca)])
                    
                Data[7].extend(sol.y[5][1:])
                Data[8].extend(sol.y[6][1:])
                Data[9].extend(sol.y[7][1:])

                if Data[2][-1]<0 or Data[2][-1]>5000000: break

                for i in range(1,len(sol.y[0])):
                    if x==0 :
                        Data[12].extend([T_aux-stage[x].nozzle_area*atmosphere(Data[2][-len(sol.y[0])+i])[1]])
                    else:
                        Data[12].extend([T_aux])


                trajectory.l[0]=sol.y[8][-1]
                trajectory.l[1]=sol.y[9][-1]
                trajectory.l[2]=sol.y[10][-1]

                
                Hf.extend(
                    [trajectory.l[0]*Data[3][-1]*sin(Data[4][-1])
                    +trajectory.l[1]*(T_aux/Data[5][-1]*sin(Data[4][-1])/Data[3][-1]
                                      +(Data[3][-1]/(Data[2][-1]+Re)-Ue/((Data[2][-1]+Re)**2*Data[3][-1]))*cos(Data[4][-1]))
                    +trajectory.l[2]*(T_aux/Data[5][-1]*cos(Data[4][-1])+(Data[3][-1]/(Data[2][-1]+Re)**2))*sin(Data[4][-1])])


                if Data[2][-1]>0 and Data[2][-1]<5000000: kn=Knudsen(L,1.225*exp(-Data[2][-1]/8440),atmosphere(Data[2][-1])[0])


                #print(x,Data[2][-1],Data[5][-1],Data[0][-1])
                #print()

                
                if Pitch==True: Data[4][-1]=trajectory.pitch; x=x-1; Pitch=False

                if kn>=10 and GT == True:
                    #print("Payload dropout")
                    payload=stage[x].payload-rocket.fairing
                    Data[5][-1]=Data[5][-1]-rocket.fairing
                    x=x-1;
                    GT=False;
                    Hamilton=True;
                    alpha_aux=0;


                if sol.status!=-1 and x!=len(stage)-1: x=x+1


                if x!=aux:
                    Data[5][-1]=payload
                    if x<=len(stage)-1:
                        Data[12][-1]=stage[x].thrust

                #print(x, Data[2][-1], T, Hamilton)

               # print(x,Data[0][-1],Data[2][-1],T,kn)
                

                #print(Hi,Hf)
#                if x==2 and Coast==False: break
                if Hamilton==True and x==len(stage)-1 and Coast==False and T!=0: break

    
    DV_required=0
    Obj=1e20
    if len(Hf)>=2:
        Obj=abs(Data[2][-1]-mission.final_altitude)**2+1e4*abs(Data[3][-1]-mission.final_velocity)**2+10e10*abs(Data[4][-1]-mission.final_flight_angle)**2+1e4*abs(Hf[-2]-Hi[-1]+Hf[-1])
    #Obj=abs(Data[2][-1]-600000)**2+100*(Data[3][-1]-5800)**2+10e7*(Data[4][-1])**2
    #print(Data[2][-1],Data[3][-1],Data[4][-1],Hf[-2]-Hi[-1]+Hf[-1], Hf[-1])
        Data[10]=Hf[-2]-Hi[-1]+Hf[-1]
        Data[11]=Hf[-1]
    #if Data[4][-1]<-0.05: Obj=1e20

#    Graphs(Data)    
    return Data,DV_required,Obj


"""

def init_trajectory(mission,rocket,stage,trajectory,general):
    Traj=True

    
    #Constants Definitions
    g0=9.810665
    Re=6371000
    Ue=398604.3*10**9

    #Initialize variables for plotting
    Time=[0]
    x=[0]
    y=[mission.initial_altitude]
    v=[mission.initial_velocity]
    phi=[mission.initial_flight_angle]
    mass=[rocket.mass]
    alpha=[0]
    vG=[0]
    vD=[0]
    vA=[0]

    Data=[Time,x,y,v,phi,mass,alpha,vG,vD,vA]





    def Hamiltonian(t,y):

        if y[1]<100000 and y[1]>0:
            D=1/2*1.225*exp(-y[1]/8440)*Area*Drag_Coeff(y[2]/sqrt(1.4*287.053*Temp_func(y[1])))*y[2]**2
        else:
            D=0


        ax=y[0]
        ay=y[1]
        V=y[2]
        gamma=y[3]
        m=y[4]

        lr=y[8]
        lg=y[9]
        lv=y[10]



        if x==0: T=T_aux-stage[0].nozzle_area*pressure_func(ay)
        else: T=T_aux


        if Hamilton==True:
            sa=-lg/(V*sqrt((lg/V)**2+lv**2))
            ca=-lv/sqrt((lg/V)**2+lv**2)
        else:
            sa=sin(alpha)
            ca=cos(alpha)


        g=Ue/(ay+Re)**2
        TM=T/m
        r=ay+Re
        
        dx=V*cos(gamma)
        dy=V*sin(gamma)
        dV=TM*ca-g*sin(gamma)-D/m
        dgamma=TM*sa/V+(V/r-g/V)*cos(gamma)
        dvg=g*sin(gamma)
        dD=D/m
        dva=TM-TM*ca

        
       # print(TM,T,m)        



        if Hamilton==True:
            dlr=V*lg*cos(gamma)/r**2-(2*Ue*lv*sin(gamma)+2*Ue*lg*cos(gamma)/V)/r**3
            dlg=-lr*cos(gamma)*V+Ue*lv*cos(gamma)/r**2+lg*sin(gamma)*(V/r-Ue/(r**2*V))
            dlv=-lr*sin(gamma)-lg*(cos(gamma)*(1/r+Ue/((r**2)*(V**2)))-TM/V**2*sa)

        else:
            dlr=0
            dlg=0
            dlv=0
            

        #print(ay)
        
        return [dx,dy,dV,dgamma,dm,dvg,dD,dva,dlr,dlg,dlv]        


    if Traj==True:

        GT=True
        Pitch=True
        Hamilton=False
        Coast=True        


        x=0
        while x!=-1:
            aux=x
            T=stage[x].thrust
            ISP=stage[x].Isp
            Area=stage[x].diameter**2/4*np.pi

            dm=-T/ISP/g0
            #print(dm)
            tb=-stage[x].propellant_mass/dm
            #Stage termination

            #Not knowing Tvac


            if GT==True: payload=stage[x].payload
            else: payload=stage[x].payload-550

            def stage_termination(t,y):

                return payload+stage[x].structural_mass-y[4]+0.001

            stage_termination.terminal=True


            def gravity_turn_termination(t,y):
                return y[1]-100000

            gravity_turn_termination.terminal=True

            if Data[2][-1]>99000: gravity_turn_termination.terminal=False


                        
            #Vertical condition
            if Data[2][-1]==0:
                tb=trajectory.vertical_time #+trajectory.pitch_time #Mudar depois

            alpha=0

            if GT==False and x<len(stage)-1:
                if alpha_aux+1>=2:
                    break
                alpha= trajectory.alpha[alpha_aux];
                alpha_aux=alpha_aux+1
     
            if x==len(stage)-1 and Coast==True: T=0; dm=0; tb=trajectory.coast_time; Coast=False


            if x==len(stage)-1 and T!=0: Hamilton=True; tb=tb*trajectory.aux


            if x==0: T_aux=T+stage[0].nozzle_area*pressure_func(0)
            else:
                T_aux=T
                
            #print(x,Data[2][-1],Data[5][-1])

            sol = solve_ivp(Hamiltonian, [0, tb], [Data[1][-1], Data[2][-1],Data[3][-1],Data[4][-1],Data[5][-1],Data[7][-1], Data[8][-1],Data[9][-1],trajectory.l[0],trajectory.l[1],trajectory.l[2]], events=[stage_termination,gravity_turn_termination], t_eval=np.linspace(0,tb,250))





            tm=np.linspace(Data[0][-1],sol.t[-1]+Data[0][-1],len(sol.t))
            #print(sol.y[1][-1])

            Data[0].extend(tm[1:])
            Data[1].extend(sol.y[0][1:])
            Data[2].extend(sol.y[1][1:])
            Data[3].extend(sol.y[2][1:])
            Data[4].extend(sol.y[3][1:])
            Data[5].extend(sol.y[4][1:])
            
            if Hamilton==False:
                Data[6].extend(np.ones(len(sol.t)-1)*alpha)
            else:
                for i in range(1,len(sol.y[0])):

                    sa=-sol.y[9][i]/(sol.y[2][i]*sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                    ca=-sol.y[10][i]/(sqrt((sol.y[9][i]/sol.y[2][i])**2+sol.y[10][i]**2))
                    Data[6].extend([atan2(ca,sa)])
                
            Data[7].extend(sol.y[5][1:])
            Data[8].extend(sol.y[6][1:])
            Data[9].extend(sol.y[7][1:])


            #print(x,Data[2][-1],Data[5][-1],Data[0][-1])
            #print()

            
            if Pitch==True: Data[4][-1]=trajectory.pitch; x=x-1; Pitch=False

            if Data[2][-1]>99000 and GT == True:
                payload=stage[x].payload-550
                Data[5][-1]=Data[5][-1]-550
                x=x-1;
                GT=False;
                alpha_aux=0;


            if sol.status!=-1 and x!=len(stage)-1: x=x+1


            if x!=aux: Data[5][-1]=payload

            #print(x, Data[2][-1], T, Hamilton)
                        
            if Hamilton==True: break
            if Data[2][-1]<0: break



    
    DV_required=0
    Obj=abs(Data[2][-1]-mission.final_altitude)**2+1e4*abs(Data[3][-1]-mission.final_velocity)**2+10e10*abs(Data[4][-1]-mission.final_flight_angle)**2
    #Obj=abs(Data[2][-1]-600000)**2+100*(Data[3][-1]-5800)**2+10e7*(Data[4][-1])**2
        
    #if Data[4][-1]<-0.05: Obj=1e20
    
    return Data,DV_required,Obj

"""
