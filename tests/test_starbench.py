"""
Starbench test
This runs the tests set by the Starbench project
See https://arxiv.org/abs/1507.05621

@author: samgeen
"""

# Import numpy and weltgeist
import numpy as np
import matplotlib.pyplot as plt

import weltgeist
import weltgeist.units as wunits # make this easier to type
import weltgeist.graphics

# For the Starbench test we use a hydrogen-only gas, so fix the units
wunits.X = 1.0
wunits.mp = wunits.mH

# Fix alphaB to be the value used in Starbench
def alphaB(Tion):
    # Value from Starbench
    return 2.7e-13
weltgeist.radiation.alpha_B_HII = alphaB

def run_test(early=True):


    # Set up the background for both tests
    rhoback = 5.21e-21 # g cm^-3
    nback = rhoback / (wunits.mH / wunits.X)
    print(wunits.X)

    # Set the temperature depending on the test
    if early:
        Tbackground = 100.0 # Kelvin
    else:
        Tbackground = 1000.0 # Kelvin

    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    # And the setup
    ncells = 10000
    nanalytic = np.zeros((ncells))
    n0 = nback # cm^-3
    T0 = Tbackground # Kelvin (will be reset later)
    gamma = 1.0001 # Weird Starbench gamma

    # What radius to set?
    if early:
        rmax = 1.5*wunits.pc
    else:
        rmax = 15.0*wunits.pc # Starbench suggests 5pc but this is too small to include the super wide shell

    integrator.Setup(ncells = ncells, # 256 cells
            rmax = rmax, # 20 pc box
            n0 = n0, # atoms / cm^-3
            T0 = T0, # Kelvin
            gamma = gamma) # weird Starbench isothermal gamma
    hydro = integrator.hydro

    # Set up the source for both tests
    QH = 1e49
    Tion = 1e4
    source = weltgeist.sources.Sources().MakeSimpleRadiation(QH,Tion=Tion) # photons/s
    # This adds the source already so no need to do it separately

    # Test doesn't include dust, so turn it off
    weltgeist.radiation.sigmaDust = 0.0

    # No radiation pressure either
    weltgeist.sources.doRadiationPressure = False

    # Force the photoionised gas to be set to Tion
    weltgeist.radiation.forceTion = True

    # Turn cooling on
    # This isn't strictly necessary
    # However, the shell spreading in the late phase is huge if cooling is not included
    #weltgeist.cooling.cooling_on = True

    # Now let's render the temperature in red
    temperatureLine = weltgeist.graphics.Line(weltgeist.graphics.red,width=3.0)
    # And a blue line for an analytic fit
    densityLine = weltgeist.graphics.Line(weltgeist.graphics.blue,width=3.0)
    # And make a rendering object again as in Exercise 3
    renderer = weltgeist.graphics.Renderer([temperatureLine,densityLine])

    # Set up a courant limiter to prevent a big jump at t=0
    #integrator.CourantLimiter(1e7)

    # Set the end time
    if early:
        endtime = 0.141 * wunits.Myr
    else:
        endtime = 3.0 * wunits.Myr

    # Arrays to allow plotting
    times = []
    radii = []

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
        time = integrator.Time()
        # Get the hydrogen number density and radius and update the line
        # Note that you can also add points to the line over time
        # This is useful for plotting time-dependent functions like example 2
        x = hydro.x[0:ncells]
        nH = hydro.nH[0:ncells]
        T = hydro.T[0:ncells]
        densityLine.Update(x,np.log10(nH))
        temperatureLine.Update(x,np.log10(T))
        # Show the time on-screen
        renderer.Text("{:.2f}".format(time/wunits.year/1e6)+" Myr")
        # Step the integrator
        integrator.Step()
        # Set the logged radii and times
        times.append(time)
        try:
            radius = hydro.x[np.where(hydro.xhii > 0.5)[-1]][-1]
            #radius = hydro.x[np.where(hydro.nH == np.max(hydro.nH))[-1]][-1]
        except:
            radius = 0.0
        radii.append(radius)
        # End if we've reached the finish point
        if time > endtime:
            raise StopIteration
    
    # Now run the renderer and the simulation will evolve!
    # Press Escape or click the cross to exit
    # Note that this doesn't have axis labels, it's super simple so far
    try:
        renderer.Start(MyStep)
    except StopIteration:
        # Stopped iterating here
        pass

    # Now plot
    plt.clf()
    times = np.array(times)
    radii = np.array(radii)
    if early:
        rhi = weltgeist.analyticsolutions.HosokawaInutsuka(QH,nback,times, Tion = Tion)
        rsp = weltgeist.analyticsolutions.SpitzerSolution(QH,nback,times, Tion = Tion)
        plt.plot(times/wunits.Myr,rhi/wunits.pc,"k--",label="Hosokawa-Inutsuka (2006)")
        plt.plot(times/wunits.Myr,rsp/wunits.pc,"k:", label="Spitzer (1978)")
        plotname = "starbench_early.pdf"
    else:
        # Note - this solution is simply extracted from the Starbench paper for comparison
        rSB = weltgeist.analyticsolutions.StarbenchExtracted(QH,nback,times, Tion = Tion)
        plt.plot(times/wunits.Myr,rSB/wunits.pc,"k--",label="Starbench Solution (Bisbas, 2015)")
        plotname = "starbench_late.pdf"
    plt.plot(times/wunits.year/1e6,radii/wunits.pc,"k-",label="Weltgeist")
    plt.xlabel("Time / Myr")
    plt.ylabel("Radius / pc")
    plt.legend(frameon=False)
    plt.savefig(plotname)

    # Reset the integrator
    integrator.Reset()

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_test(early=False)
    run_test(early=True)
    
