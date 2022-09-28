"""
Example 6 - Loading
Loading a save file

@author: samgeen
"""

# Import numpy, matplotlib and weltgeist
import numpy as np
import matplotlib.pyplot as plt
import weltgeist
import weltgeist.units as wunits # make this easier to type

def run_example():
    # This is just example 5 in reverse, kinda
    integrator = weltgeist.integrator.Integrator()
    # If you haven't run example 5 yet, do it now so you have this file
    # You don't need the ".hdf5", it'll add that for you
    print("Loading...")
    integrator.Load("mycoolsave")
    print("Loaded!")
    # Okay! We loaded a file
    # AGAIN - this doesn't contain any information on sources
    # I swear this will happen sometime
    # Serialising nested objects is time-consuming, is all

    # What's the current simulation time?
    print(str(integrator.time/wunits.Myr)+" Myr")
    # Alright! It's not zero, which means something happened. But is the simulation state ok?
    
    # Let's try plotting the same thing as the end of example 5 when we saved it
    # I made this a function so I can do it again
    # Python does some wild stuff with nested functions and scope
    hydro = integrator.hydro
    ncells = hydro.ncells
    def plotstuff():
        # OK now plot something to show it's evolved
        plt.clf()
        plt.plot(hydro.x[0:ncells]/wunits.pc,hydro.nH[0:ncells],label=str(integrator.time/wunits.Myr)+" Myr",color="r")
        plt.xlabel("radius / pc")
        plt.ylabel("$n_{\mathrm{H}}$ / cm$^{-3}$")
        plt.yscale("log")
        plt.legend(frameon=False)
        plt.show()
    plotstuff()
    
    # Run it a bunch
    niterations = 100
    for i in range(0,niterations):
        # Bootleg timer
        if i % 10 == 0:
            print(".")
        integrator.Step()

    # Plot again
    plotstuff()
    # Okay. We loaded a simulation and kept running it

    # And that's it. That's all the examples I could think of.
    # You can explore the code if you like to see how it all works
    # Gravity is also included, but I need to confirm it works with a
    #  test problem, so that's for another time

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
