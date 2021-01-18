"""
Defined code units and physical quantities
The Python parts of Weltgeist use cgs
VH1 uses units defined below
Sam Geen, February 2018
"""

import numpy as np

# Physical quantities (base units in cgs)
pc = 3.086e+18
mH = 1.66e-24
year = 3.154e+7
Myr = 1e6*year
kB = 1.3806485279e-16 # in cgs
G = 6.67428e-8
X = 0.74
mp = mH / X
c = 2.998e+10
eV = 1.60217662e-12 # in ergs
Msun = 1.9891e33 # g

# Code units
# Used by VH1 - the Python parts of Weltgeist use cgs
distance = pc # in cm
density = mH # 1 g/cm^3
time = 1.0 / np.sqrt(G*density) # sets G=1 in VH1 (not super important here, though)
# Derived units
velocity = distance / time
mass = density*distance**3.0
pressure = density * velocity**2.0
energy = mass*velocity**2.0
# Note: this is acceleration! In the code (e.g. forces.f90), grav = v*v/r
# e.g. 2*GM/r = v_esc^2, so g=GM/r^2=0.5*v_esc^2/r
gravity = G*mass/distance**2 # velocity*velocity/distance 
