################ Classes ##################
class a:
    def d(self):
        pass
b=a()

#Trajectory Class
class Trajectory:
    def vertical_time(self):
        pass

    def vertical_altitude(self):
        pass


    def pitch_time(self):
        pass

    def pitch(self):
        pass

    def delta_pitch(self):
        pass


    #Define boundaries as vector to save time.

    def l2_ub(self):
        pass

    def l2_lb(self):
        pass

    def l3_ub(self):
        pass

    def l3_lb(self):
        pass



#Optimization class

class Optimization:
    def DV_Optimization(self):
        pass

    def Trajectory_Optimization(self):
        pass



#General Class
class General:

    def Trajectory(self):
        pass

    def Comparison(self):
        pass


    def TWR(self):
        pass

#Mission class
class Mission:
    def final_altitude(self):
        pass

    def final_velocity(self):
        pass

    def initial_altitude(self):
        pass

    def initial_velocity(self):
        pass

    def initial_flight_angle(self):
        pass

    def final_flight_angle(self):
        pass

    def payload(self):
        pass

    def objective(self):
        pass



#Rocket class

class Rocket:
    def num_stages(self):
        pass

    def num_boosters(self):
        pass

    def payload(self):
        pass

    def mass(self):
        pass

    def DV(self):
        pass

    def DV_History(self):
        pass
    
#Stage class

# ATENÇÂO     Não inclui nozzle Area nem Ratio Area

class Stage:
    def stage_num(self):
        pass
    
    def type(self):
        pass

    def diameter(self):
        pass

    def length(self):
        pass

    def num_engines(self):
        pass

    def thrust(self):
        pass

    def Isp(self):
        pass

    def structural_mass(self):
        pass

    def propellant_mass(self):
        pass

    def total_mass(self):
        pass
    
    def payload(self):
        pass

    def OF_ratio(self):
        pass

    def fuel_density(self):
        pass

    def oxidizer_density(self):
        pass

    def DV(self):
        pass

    def TWR(self):
        pass

    def DV_Ratio(self):
        pass
