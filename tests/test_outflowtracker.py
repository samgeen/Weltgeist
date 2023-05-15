"""
Test the outflow tracker

@author: samgeen
"""

# Import numpy and weltgeist
import numpy as np
import weltgeist
import weltgeist.units as wunits # make this easier to type

import weltgeist.graphics # Separate in case we don't want to import it

# Is our source turned on?
sourceOn = True

def run_example(radiation=False,winds=False,testLoading=False):
    global sourceOn
    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 8
    nanalytic = np.zeros((ncells))
    n0 = 10.0 # cm^-3
    T0 = 10.0 # K
    rmax = 0.1*wunits.pc
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = rmax, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    hydro = integrator.hydro

    backgroundMass = n0 * wunits.mH / wunits.X * 4.0 / 3.0 * np.pi * rmax**3

    # Make a simple radiation source
    if radiation:
        Qionsource = 1e49
        source = weltgeist.sources.SimpleRadiationSource(Qionsource)
    elif winds:
        windLuminosity = 1e36
        windMassloss = 1e20
        source = weltgeist.sources.WindSource(windLuminosity,windMassloss)
    else:
        print("Pick eiter winds or radiation in run_example")
        raise ValueError
    weltgeist.sources.Sources().AddSource(source)

    # Turn cooling and radiation if we are doing a radiation test
    # Otherwise keep things adiabatic with no radiation
    weltgeist.cooling.cooling_on = radiation
    weltgeist.radiation.radiation_on = radiation
    # Stop the temperature doing weird things
    weltgeist.radiation.forceTion = radiation

    # Now let's render the temperature in red
    temperatureLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for density
    densityLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # And a black line for ionisation fraction
    ionLine = weltgeist.graphics.Line(weltgeist.graphics.black,width=3.0)
    # And make a rendering object
    renderer = weltgeist.graphics.Renderer([temperatureLine,densityLine,ionLine])

    # Pick a time to stop the source and simulation
    stopSourceTime = 1e11
    stopSimulationTime = 10*stopSourceTime

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
        time = integrator.Time()
        x = hydro.x[:]
        nH = hydro.nH[:]
        # Note - We add no thermal energy so this will be a stupidly low value, that's okay
        T = hydro.T[:]
        hydro.nH[hydro.nH < 1e-10] = 1e-10
        xhii = hydro.xhii[0:ncells]*4
        densityLine.Update(x,np.log10(nH))
        temperatureLine.Update(x,np.log10(T))
        ionLine.Update(x,xhii)
        # Set up a courant limiter to prevent a big jump at t=0
        integrator.CourantLimiter(1e6)
        # Show the time on-screen
        renderText = "{:.2f}".format(time/wunits.year/1e6)+" Myr"
        if not sourceOn:
            renderText = "SOURCE OFF "+renderText
        renderer.Text(renderText)
        # Turn off source after 1 Myr
        if time > stopSourceTime and sourceOn:
            weltgeist.sources.Sources().RemoveSource(source)
            sourceOn = False
        # Step the integrator
        integrator.Step()
        # End if we've reached the finish point
        if time > stopSimulationTime:
            raise StopIteration
    
    # Now run the renderer and the simulation will evolve
    # Press Escape or click the cross to exit, or wait until the end time
    try:
        renderer.Start(MyStep)
    except StopIteration:
        # Stopped iterating here
        pass

    # Output log file for times
    if radiation:
        Nemitted = Qionsource * stopSourceTime
        print("Number of photons that were emitted by the source:", Nemitted)
    if winds:
        Eemitted = windLuminosity * stopSourceTime
        Memitted = windMassloss * stopSourceTime
        # Momentum = sqrt(2 * 1/2 m v^2 * m)
        momentum = np.sqrt(2.0 * Eemitted * Memitted)
        print("Mass (+background), momentum and energy emitted by source: ", Memitted+backgroundMass, momentum, Eemitted)


    # Print reporting from the outflow tracker to see how much
    #  mass, momentum, energy and photons were lost from the grid
    print ("Flows lost to the grid - there will be some things lost to absorption in the gas...")
    print(integrator.outflowTracker)

    if testLoading:
        filename = "savefile_outflowtest"
        integrator.Save(filename)
        integrator.Load(filename)
        print ("Testing saving and loading - should be the same as printed above:")
        print(integrator.outflowTracker)

    # Reset integrator
    integrator.Reset()


# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example(winds=True, testLoading=True)
    run_example(radiation=True, testLoading=True)
    
