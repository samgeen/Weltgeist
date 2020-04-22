"""
Example 2 - Leonid Sedov
Run a simple test problem and compare it to an analytic solution

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
    # First let's pick some variables
    ncells = 256
    niterations = 1000
    #  We need the integrator again
    integrator = weltgeist.integrator.Integrator()
    weltgeist.cooling.cooling_on = True
    # And the setup
    integrator.Setup(ncells = ncells, # 256 cells
            rmax = 10.0*wunits.pc, # 20 pc box
            n0 = 1000.0, # 1000 H atoms / cm^-3
            T0 = 10.0, # 10 K
            gamma = 5.0/3.0) # monatomic gas (close enough...)
    # Now let's grab the hydrodynamic variables and take a look
    hydro = integrator.hydro
    # The hydro object contains all the variables on the grid
    # Lets print a quick view of the first 10 cells
    ncells = hydro.ncells
    print("Number of cells: ", ncells)
    print("Hydrogen number density in cm^-3: ", hydro.nH[0:10])
    print("Temperature in K: ", hydro.T[0:10])
    print("Thermal Energy in ergs: ", hydro.TE[0:10])
    # Note that the thermal energy goes up because the cell mass is not
    #  uniform - the cell size dr is constant, so the cell mass is not
    # (mass ~ 4 pi r^2 dr, and r goes up)
    
    # We can also modify the grid - let's test this with a supernova
    # First, you can modify any cell, any time
    # NOTE: for now you need to subscript like this to prevent 
    #  overwiting the variable - fixing this is on my todo list!
    newnH = 100.0
    hydro.nH[0:ncells] = newnH
    # We can also extract derived quantities
    # hydro tries to keep things consistent internally
    newrho = hydro.rho[0]
    # Note that changing the density also changes the temperature,
    #  so let's fix that
    T0 = 10.0 # K
    hydro.T[0:ncells] = T0 

    # Add 1e51 ergs to the first cell's kinetic energy
    # This is about the energy of a supernova
    snenergy = 1e51
    hydro.KE[0] += snenergy
    # In the next example we'll look at the sources module,
    #  which does this properly
    print("Velocity in cm/s: ", hydro.vel[0:10])

    # Now if we let this go, it'll expand
    # Let's step this a bunch and track its progress
    times = []
    radii = []
    for i in range(0,niterations):
        # Bootleg timer
        if i % 100 == 0:
            print(".")
        integrator.Step()
        times.append(integrator.time)
        # Get the position of the furthest cell where a wave is found
        edgecell = hydro.x[np.where(hydro.T[0:ncells] > 2.0*T0)[-1]][-1]
        radii.append(edgecell)
    # Let's turn these lists into numpy arrays
    times = np.array(times)
    radii = np.array(radii)

    # We can compare this to the "Sedov-Taylor" solution:
    # http://www.astronomy.ohio-state.edu/~ryden/ast825/ch5-6.pdf
    # This is named after Leonid Sedov and Geoffrey Taylor
    # Both discovered it independently in the 1950s because strangely,
    #  the USSR and UK were not sharing blastwave research at the time
    # You can get the solution 90% right by doing dimensional analysis
    radii_sedov = weltgeist.analyticsolutions.SedovTaylorSolution(times, snenergy, newrho)

    # Now let's plot all of this
    plt.clf()
    # Fair warning - I will ban you from using red and green on the same plot
    # Accessibility is important!
    # Tip: go to colorbrewer2.org, check "colorblind safe" and
    #  choose a colormap from that list
    plt.plot(times/wunits.Myr,radii/wunits.pc,label="Simulation",color="r")
    plt.plot(times/wunits.Myr,radii_sedov/wunits.pc,label="Sedov-Taylor",color="b")
    plt.xlabel("Time / Myr")
    plt.ylabel("radius / pc")
    plt.legend()
    plt.show()

    # The simulation line isn't exactly on the analytic solution.
    # There are also steps in the curve
    # This is the price of finite numerical resolution and CPU time
    # Try running this exercise again with niterations=10000 and ncells=512
    # Does the solution match better?
    # You can also play setting hydro.TE (thermal energy) instead of KE
    #  and adding some mass to the central cell as ejecta (hydro.mass[0] = 2e33)

    # In the next module we look at sources and visualisation

# This piece of code runs if you start this module versus importing it
if __name__=="__main__":
    run_example()
    
