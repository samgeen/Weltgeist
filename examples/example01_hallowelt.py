"""
Example 1 - Hallo, Welt!
Get the code running for one timestep

@author: samgeen
"""

# This piece of code basically adds the parent directory to PYTHONPATH
import os, sys
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)

# Import numpy and weltgeist
import numpy as np
import weltgeist
import weltgeist.units as wunits # make this easier to type

def run_example():
    # First we need the integrator object
    # This controls the simulation setup and running
    # weltgeist is the package name (folder weltgeist/)
    # integrator is the module name (file integrator.py)
    # Integrator() is the function that returns the object
    # Since f2py doesn't like instancing multiple versions of Fortan
    #  modules, we're stuck with the singleton pattern
    # (https://en.wikipedia.org/wiki/Singleton_pattern)
    integrator = weltgeist.integrator.Integrator()
    # Let's go ahead and set up a simulation volume
    # This runs the Setup method of the integrator object we got from 
    #  the weltgeist package
    integrator.Setup(ncells = 256, # 256 cells
            rmax = 20.0*wunits.pc, # 20 pc box
            n0 = 1000.0, # 1000 H atoms / cm^-3
            T0 = 10.0, # 10 K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    # Now let's make the first timestep
    integrator.Step()
    # Now the current time should be something like 100 kyr
    # This is because of the Courant condition: 
    # https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
    # The box is a uniform gas at 10 K, and a cell is 20 pc / 256,
    #  so the crossing time of a cell is ~ 100 kyr 
    print("Current time: ", integrator.time, "seconds")
    print("Timestep: ", integrator.dt, "seconds")
    integrator.Step()
    # Now the time is twice that again.
    print("Current time: ", integrator.time, "seconds")
    # That's it for the first example. Next we'll look at a simple test problem.


# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
