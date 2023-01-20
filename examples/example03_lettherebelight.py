"""
Example 3 - Let There Be Light
Explore radiation sources and visualisation

@author: samgeen
"""

# Import numpy and weltgeist
import numpy as np
import weltgeist
import weltgeist.units as wunits # make this easier to type

import weltgeist.graphics # Separate in case we don't want to import it

def run_example():
    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 256
    nanalytic = np.zeros((ncells))
    n0 = 1000.0 # cm^-3
    T0 = 10.0 # K
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    # Weltgeist has various in-built ways to simulate stars
    # We can use the source module to put energy onto the grid
    # This can be winds, radiation, supernovae
    QH = 1e49 # photons per second emitted, roughly a 35 solar mass star
    weltgeist.sources.Sources().MakeSimpleRadiation(QH)

    # Weltgeist has a real-time visualisation tool for checking quickly
    #  what the code is doing
    # First, let's make a red line to render the simulation to screen
    simulationLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for an analytic fit
    analyticLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # Now we make a rendering object to draw them
    # ("Rendering" just means drawing a picture to screen)
    renderer = weltgeist.graphics.Renderer([simulationLine,analyticLine])

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
        # Get the hydrogen number density and radius and update the line
        # Note that you can also add points to the line over time
        # This is useful for plotting time-dependent functions like example 2
        x = hydro.x[0:ncells]
        nH = hydro.nH[0:ncells]
        simulationLine.Update(x,nH)
        # Now let's get the analytic function to compare to
        # This is called the "Spitzer" solution after Lyman Spitzer
        # It's in his 1978 book on the Interstellar Medium
        nanalytic[:ncells] = n0
        time = integrator.time
        ri = weltgeist.analyticsolutions.SpitzerSolution(QH,n0,time)
        ni = weltgeist.analyticsolutions.SpitzerDensity(QH,n0,time)
        nanalytic[np.where(x <= ri)] = ni
        # NOTE: this won't predict the shell density
        # Exercise: using the Rankine-Hugoniont conditions, include this
        analyticLine.Update(x,nanalytic)
        # Show the time on-screen
        renderer.Text("{:.2f}".format(time/wunits.Myr)+" Myr")
        # Step the integrator
        integrator.Step()
    
    # Now run the renderer and the simulation will evolve!
    # Press Escape or click the cross to exit
    # Note that this doesn't have axis labels, it's super simple so far
    renderer.Start(MyStep)

    # A few more things to try:
    # 1) Change T0 to 3000 K. 
    # What happens? Why does the underdense region stop?
    # Answer: see Raga
    # 2) Add cooling - uncomment this line and move it above renderer.Start(MyStep)
    # weltgeist.cooling.cooling_on = True
    # This adds radiative cooling to the gas
    # Basically, various processes cause the gas to lose thermal energy
    #  to low-energy radiation
    # What happens? Why is the shell thinner?
    # 3) Explore what happens before the solution is stable
    #    Plot hydro.xhii (hydrogen ionisation fraction)
    #    The above solutions assume photoionisation equilibrium ("D-type")
    #    What happpens before that? ("R-type")

    # In the next example we look at more "realistic" stars and clouds

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
