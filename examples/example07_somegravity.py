"""
Example 7 Some Gravity
Test out the gravity module

@author: samgeen
"""

# Import numpy and weltgeist
import numpy as np
import weltgeist
import weltgeist.units as wunits # make this easier to type

import weltgeist.graphics # Separate in case we don't want to import it

def run_example():
    # You like dinosaurs, right? 
    # Or maybe you prefer space. That probably makes more sense now I think about it.
    # Anyway
    # This example looks like a dinosaur's tail
    # I think so, anyway. You'll see.

    # This time we're going to explore gravity
    # Because this is 1D, it always points towards the centre of the grid at r=0

    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 512
    # Fiducial very low density background
    n0 = 1e-2 # cm^-3
    T0 = 100000.0 # K (this is quite hot but prevents the gas in the centre doing weird stuff)
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    # Let's set up a simple test problem with a big gas mass in the middle and a collapsing ring around it

    # Add a central mass of 1e5 Msun just to be absurd
    weltgeist.gravity.centralmass = 1e5 * wunits.Msun 
    # Initial position of mass to fall from in grid indices (halfway)
    initialindex = ncells//2
    # Add a ring of mass 1 Msun to that position
    hydro.mass[initialindex] = 1e1 * wunits.Msun
    # And let's get the initial radius (this is important for the analytic solution later)
    startingposition = hydro.x[initialindex]

    # I'm forgetting something...
    # Ah yes, we need to turn gravity on
    weltgeist.gravity.gravity_on = True

    # Let's render what's happening
    # First, let's make a red line to render the simulation to screen
    simulationLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for an analytic fit
    analyticLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # Now we make a rendering object to draw them
    # You've seen this before, of course
    renderer = weltgeist.graphics.Renderer([simulationLine,analyticLine])

    # Set up some arrays to store values in the simulations over time
    # You may notice that we predict time from the peak position
    # This is just the way the analytic solution works 
    times = []
    analytictimes = []
    peakpositions = []

    # Because of the way the rendering module pyglet works, it has to
    #  control the stepping. So let's make a function to give it
    def MyStep(dtRender):
        """
        Function for pyglet to run every timestep

        Parameters
        ----------

        dtRender : float
            wall clock timestep each frame of the rendering
        """
        # This time we'll plot the position of the ring of mass (defined by the peak density) over time
        # The x-axis is in Myr, the y-axis is in parsecs
        times.append(integrator.time / wunits.Myr)
        peakposition = hydro.x[np.argmax(hydro.nH)]
        # Update the sequence of 
        peakpositions.append(peakposition / wunits.pc)
        simulationLine.Update(times,peakpositions)

        # If you want to see what is happening to the gas field, 
        #   replace the line above with this one (and stop the analytic line from updating)
        #simulationLine.Update(hydro.x[:],hydro.nH[:])
        
        # Get rid of any background gas pooling at the centre
        hydro.nH[0:20] = n0
        hydro.T[0:20] = hydro.T[-1]
        # Steadily remove thermal energy to prevent over-pressurisation of the gas, which would slow down freefall
        hydro.T *= 0.90

        # Set up the analytic solution
        # As I said above, the solution predicts time from position
        # Check out the analyticsolutions module for the equation used and the source
        analytictime = weltgeist.analyticsolutions.CollapseSolutionPosition(peakposition,startingposition)
        analytictimes.append(analytictime / wunits.Myr)
        analyticLine.Update(analytictimes,peakpositions)

        # Step the integrator
        integrator.Step()
    
    # Now run the renderer and the simulation will evolve!
    # It might move slowly at first - it takes time for the most dense cell to get moving
    # Over time, the ring falls to the centre
    # Because the grid has a finite resolution, initially the peak follows integer cells, hence the "dinosaur tail"
    # I dunno, I thought it was a fun visual metaphor
    renderer.Start(MyStep)
    # Press Escape or click the cross to exit

    # Some final thoughts:
    # This setup is pretty idealised and kind of unstable.
    # When I set this up, it did all kinds of weird things until I found a "stable" version
    # In reality, gas can be supported by turbulence and/or rotation
    # And the problem of the unstable singularity at the origin (x=0) is less of an issue
    # But still, gravity can be important - compressing shells, stalling expansion
    # Maybe you have ideas for how you want to play with this
    # Or you can just ignore it and turn it off. That's probably fine too.

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
