"""
Example 4 - Starmaker
Explore stellar sources and initial conditions

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
    ncells = 128
    nanalytic = np.zeros((ncells))
    n0 = 100.0 # cm^-3
    T0 = 10.0 # K
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 200.0*wunits.pc, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    # Real stars don't form in completely uniform environments
    # The cores that stars form in have roughly(!) "isothermal" profiles
    # This means that density is proportional to 1/r^2
    # Let's modify the initial density field to reflect this
    r0 = 1.0 * wunits.pc
    hydro.nH[0:ncells] = n0 * (hydro.x[0:ncells] / r0)**(-2.0)
    # Get rid of the singularity at r=0
    hydro.nH[0] = hydro.nH[1]
    # We also have to set the temperature correctly again
    hydro.T[0:ncells] = T0
    # Note that this setup is a bit unstable - as long as feedback 
    #  responds faster, that's OK...

    # Now let's add a star. 
    # First, and this is very important, set the table location
    # This tells the code where the stellar evolution tracks are
    # This will be different for your computer!
    # If you don't have the tables, email me
    # These tables are for Solar metallicity (Z = 0.014)
    # There are also tables for sub-Solar (Z = 0.002)
    # Depends if you want to model our Galaxy or the LMC/SMC
    weltgeist.sources.singlestarLocation = \
        "../StellarSources/data/singlestar_z0.014"

    # Second, make a star using these tables
    # This is a 30 solar mass star 
    # By default it has all the feedback modes turnes on
    # You can turn them off in the function below
    # e.g. star = TableSource(30.0,radiation=False,wind=True)
    star = weltgeist.sources.TableSource(119.0,radiation=True,wind=True,supernova=True)
    weltgeist.sources.Sources().AddSource(star)

    # Turn cooling on
    weltgeist.cooling.cooling_on = True

    # Now let's render the temperature in red
    temperatureLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for an analytic fit
    densityLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # And make a rendering object again as in Exercise 3
    renderer = weltgeist.graphics.Renderer([temperatureLine,densityLine])

    # Set up a courant limiter to prevent a big jump at t=0
    integrator.CourantLimiter(1e7)

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
        T = hydro.T[0:ncells]
        densityLine.Update(x,np.log10(nH))
        temperatureLine.Update(x,np.log10(T))
        # Show the time on-screen
        renderer.Text("{:.2f}".format(integrator.Time()/wunits.year/1000)+" kyr")
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
    
