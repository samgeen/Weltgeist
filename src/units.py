'''
Defined code units and physical quantities
Sam Geen, February 2018
'''

# Physical quantities (base units in cgs)
pc = 3.086e+18
mH = 1.66e-24
g = mH
year = 3.154e+7
kB = 1.3806485279e-16 # in cgs
G = 6.67428e-8
X = 0.74
mp = mH / X

# Units (only really used internally to make VHone not complain about numerical accuracy)
distance = pc # in cm
density = g # = 1 g/cm^3
time = 1e6*year
# Derived units
velocity = distance / time
mass = density*distance**3.0
pressure = density * velocity**2.0
energy = mass*velocity**2.0
# Note: this is acceleration! In the code (e.g. forces.f90), grav = v*v/r
# e.g. 2*GM/r = v_esc^2, so g=GM/r^2=0.5*v_esc^2/r
gravity = velocity*velocity/distance 