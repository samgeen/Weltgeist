"""
Test recombination in an example source

@author: samgeen
"""

# Import numpy and weltgeist
import numpy as np
import weltgeist
import weltgeist.units as wunits # make this easier to type

import weltgeist.graphics # Separate in case we don't want to import it

# Is our source turned on?
sourceOn = True

def run_example():
    global sourceOn
    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 512
    nanalytic = np.zeros((ncells))
    n0 = 1000.0 # cm^-3
    T0 = 10.0 # K
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    # Make a simple radiation source
    source = weltgeist.sources.SimpleRadiationSource(1e49)
    weltgeist.sources.Sources().AddSource(source)

    # Turn cooling and radiation on
    weltgeist.cooling.cooling_on = True
    weltgeist.radiation.radiation_on = True

    # Now let's render the temperature in red
    temperatureLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for density
    densityLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # And a black line for ionisation fraction
    ionLine = weltgeist.graphics.Line(weltgeist.graphics.black,width=3.0)
    # And make a rendering object
    renderer = weltgeist.graphics.Renderer([temperatureLine,densityLine,ionLine])


    # Because of the way the rendering module pyglet works, it has to
    #  control the stepping. So let's make a function to give it
    def MyStep(dtRender):
        global sourceOn
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
        T = hydro.T[0:ncells]
        xhii = hydro.xhii[0:ncells]*4
        densityLine.Update(x,np.log10(nH))
        temperatureLine.Update(x,np.log10(T))
        ionLine.Update(x,xhii)
        # Set up a courant limiter to prevent a big jump at t=0
        integrator.CourantLimiter(1e6)
        # Show the time on-screen
        renderText = "{:.2f}".format(integrator.Time()/wunits.year/1e6)+" Myr"
        if not sourceOn:
            renderText = "SOURCE OFF "+renderText
        renderer.Text(renderText)
        # Turn off source after 1 Myr
        if integrator.Time() > 1e6*wunits.year and sourceOn:
            weltgeist.sources.Sources().RemoveSource(source)
            sourceOn = False
        # Step the integrator
        integrator.Step()
    
    # Now run the renderer and the simulation will evolve!
    # Press Escape or click the cross to exit
    # Note that this doesn't have axis labels, it's super simple so far
    renderer.Start(MyStep)

    # A few things to notice
    # 1) Winds are a giant pain. They're super fast (up to a few 
    #  1000 km/s) and they make super hot gas (up to 10^8-10^9 K) 
    # This is an issue for the Courant condition. A 3D simulation that 
    # takes 2 days without winds takes weeks with winds included
    # 2) What structures do you see evolving? How do they compare with
    #  https://ui.adsabs.harvard.edu/abs/1977ApJ...218..377W/abstract
    # Try turning radiation on to see what happens. Why?
    # Try plotting xHII instead of log(T) (or add a new line for it)
    # Now try setting the medium to uniform and see what happens

    # In the next example we explore saving the simulation state to file

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
