'''
Defined code units and physical quantities
Sam Geen, February 2018
'''

# Physical quantities (base units in cgs)
pc = 3.086e+18
g = 1.66e-24
year = 3.154e+7
kB = 1.3806485279e-16 # in cgs

X = 0.74

# Units (only really used internally to make VHone not complain about numerical accuracy)
distance = pc # in cm
density = g/1.0**3.0 # = g/cm^3
time = 1e6*year
# Derived units
velocity = distance / time
mass = density*distance**3.0
pressure = density * velocity**2.0
energy = mass*velocity**2.0