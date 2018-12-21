import math

def Grav_acceleration(h):
    g=6.67408*10**-11*5.97219*10**24/(h+6371000)**2
    return g

def Grav_Force(h,Mass,Gamma):
    Fg=Grav_acceleration(h)*Mass*math.sin(Gamma)
    return Fg

