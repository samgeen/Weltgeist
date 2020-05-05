"""
Example 5 - Saving
Save the state of the simulation

@author: samgeen
"""

# This piece of code basically adds the parent directory to PYTHONPATH
import os, sys
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)

# Import numpy, matplotlib and weltgeist
import numpy as np
import matplotlib.pyplot as plt
import weltgeist
import weltgeist.units as wunits # make this easier to type

def run_example():
    # This is mostly going to be a straight copy of example 2 to get some data we can use
    # First let's pick some variables
    ncells = 128
    niterations = 100
    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    weltgeist.cooling.cooling_on = True
    # And the setup
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = 100.0, # 1000 H atoms / cm^-3
            T0 = 10.0, # 10 K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro
    ncells = hydro.ncells

    # Make a supernova
    snenergy = 1e51
    hydro.KE[0] += snenergy

    # Run it a bunch
    for i in range(0,niterations):
        # Bootleg timer
        if i % 10 == 0:
            print(".")
        integrator.Step()
    
    # OK now plot something to show it's evolved
    plt.clf()
    plt.plot(hydro.x[0:ncells]/wunits.pc,hydro.nH[0:ncells],label=str(integrator.time/wunits.Myr)+" Myr",color="r")
    plt.xlabel("radius / pc")
    plt.ylabel("$n_{\mathrm{H}}$ / cm$^{-3}$")
    plt.yscale("log")
    plt.legend(frameon=False)
    plt.show()

    # OK, great, now the point of the example
    integrator.Save("mycoolsave")
    # WARNING - this will not save sources (yet) - this is a bit harder to do, but I'm working on it...
    # By which I mean I punted it down to the bottom of my unfeasibly long todo list
    # It'll happen sometime

    # In the next module we look at loading this file

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
