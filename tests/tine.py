import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colorsx
import matplotlib.cm as cmx
import os

import weltgeist
import weltgeist.units as wunits

DT_OUTPUT = 0.5

''' helper function to create a simulation directory for parameter study '''
def construct_directory_name(simulation_setup):
    star_mass  = simulation_setup['star mass']
    n0         = simulation_setup['central density']
    alpha      = simulation_setup['density profile slope']
    supernova  = simulation_setup['supernova']
    radiation  = simulation_setup['radiation']
    wind       = simulation_setup['winds']
    boxlen     = simulation_setup['boxlen']
    ncells     = simulation_setup['ncells']

    return 'M{}_n{}_a{}_SN{}_R{}_W{}_l{}_n{}/'.format(star_mass, n0, alpha, int(supernova), int(radiation), int(wind), int(boxlen), int(ncells))

''' run the simulation '''
def run_simulation(simulation_setup):

    # integrator
    integrator = weltgeist.integrator.Integrator()
    weltgeist.cooling.cooling_on = True
    weltgeist.cooling.maskContactDiscontinuity = False

    # variable parameters
    star_mass = simulation_setup['star mass']
    n0        = simulation_setup['central density']
    alpha     = simulation_setup['density profile slope']
    supernova = simulation_setup['supernova']
    radiation = simulation_setup['radiation']
    wind      = simulation_setup['winds']
    boxlen     = simulation_setup['boxlen']
    ncells    = simulation_setup['ncells']

    # fixed parameters for simulation suite
    T0     = 100 #K
    gamma = 5.0/3.0 # monatomic gas (close enough...)

    # make output directory
    output_dir = construct_directory_name(simulation_setup)
    if not os.path.exists(os.path.dirname(output_dir)):
        os.makedirs(os.path.dirname(output_dir))

    # setup initial conditions
    integrator.Setup(ncells = ncells,
                     rmax = boxlen*wunits.pc,
                     n0 = n0,
                     T0 = T0,
                     gamma = gamma)
    # modify density to get powerlaw profile
    if alpha != 0:
        r0 = 1.0 * wunits.pc
        integrator.hydro.nH[0:ncells] = n0 * (integrator.hydro.x[0:ncells] / r0)**(-alpha)
        integrator.hydro.nH[0] = integrator.hydro.nH[1]   # Get rid of the singularity at r=0

    # add the central star
    weltgeist.sources.singlestarLocation = "../StellarSources/data/singlestar_z0.014"
    star = weltgeist.sources.TableSource(star_mass, supernova=supernova, radiation=radiation, wind=wind)
    weltgeist.sources.Sources().AddSource(star)

    # output initial conditions
    output_number = 0
    integrator.Save("{}/output_{:05}".format(output_dir, output_number))
    output_number += 1

    # get sn time
    integrator.Step() # hack
    sn_time=star._supernovaTime/wunits.Myr
    end_time = int(sn_time) + 2#8
    print(sn_time, end_time)
    np.savetxt('SN_time_M{}.txt'.format(star_mass), np.array([sn_time]))

    # run 
    time = 0
    times = []
    sn_shock = []
    rad_shock = []
    momentum = []
    outwards_momentum = []
    while time <= end_time:
        # progress some amount of steps
        for s in range(50):
            integrator.Step()
            # output snapshot at specific times
            time = integrator.time/wunits.Myr
            if time >= output_number*DT_OUTPUT:
                print('output', output_number, 'at time', time, 'dt=', integrator.dt/wunits.Myr)
                integrator.Save("{}/output_{:05}".format(output_dir, output_number))
                output_number += 1

        print('time', integrator.time/wunits.Myr, 'dt=', integrator.dt/wunits.Myr)

        # calculate some intermediate results
        times.append(time)
        # determine HII bubble radius
        if simulation_setup['radiation']:
            try:
                radiation_bubble_radius =  integrator.hydro.x[np.where( integrator.hydro.xhii[0:ncells] > 0.95)[-1]][-1]
            except:
                radiation_bubble_radius = 0
            rad_shock.append(radiation_bubble_radius)
        # determine SN bubble radius
        if simulation_setup['supernova']:
            try:
                sn_bubble_radius =  integrator.hydro.x[np.where( integrator.hydro.T[0:ncells] > 1e5)[-1]][-1]
            except:
                sn_bubble_radius = 0
            sn_shock.append(sn_bubble_radius)
        # calculate total momentum
        momentum.append(np.sum(integrator.hydro.mass[0:ncells] * integrator.hydro.vel[0:ncells]))
        outwards_momentum.append(np.sum(integrator.hydro.mass[0:ncells][integrator.hydro.vel[0:ncells] > 0] * integrator.hydro.vel[0:ncells][integrator.hydro.vel[0:ncells] > 0]))
        
        # write diagnostics
        np.savetxt(output_dir+'/time.txt', np.array(times))
        np.savetxt(output_dir+'/sn_radius.txt', np.array(sn_shock) / wunits.pc)
        np.savetxt(output_dir+'/HII_radius.txt', np.array(rad_shock) / wunits.pc)
        np.savetxt(output_dir+'/momentum.txt', np.array(momentum))
        np.savetxt(output_dir+'/outwards_momentum.txt', np.array(outwards_momentum))

    # write final diagnostics
    np.savetxt(output_dir+'/time.txt', np.array(times))
    np.savetxt(output_dir+'/sn_radius.txt', np.array(sn_shock) / wunits.pc)
    np.savetxt(output_dir+'/HII_radius.txt', np.array(rad_shock) / wunits.pc)
    np.savetxt(output_dir+'/momentum.txt', np.array(momentum))
    np.savetxt(output_dir+'/outwards_momentum.txt', np.array(outwards_momentum))


''' helper function for mapping colors to time '''
def define_colors(n_outputs):
    cmap = cmx.get_cmap('inferno')
    colors = []
    cNorm  = colorsx.Normalize(vmin=0, vmax=n_outputs*1.05)
    for val in range(n_outputs):
        colors.append(cmap(cNorm(val)))
    return colors

''' helper function to retreive SN time '''
def load_sn_time(mstar):
    try:
        sn_time = np.loadtxt('SN_time_M{}.txt'.format(mstar))
    except:
        sn_time = 0
    return sn_time

''' Show time evolution of density, temperature and ionisation profile '''
def plot_time_evolution(simulation_setup, what='output'):
    # gather outputs
    output_dir = construct_directory_name(simulation_setup)
    n_outputs = 0
    for item in os.listdir(output_dir):
        if item.startswith(what):
            n_outputs += 1

    integrator = weltgeist.integrator.Integrator()

    colors = define_colors(n_outputs)
    sn_time=load_sn_time(simulation_setup['star mass'])

    fig, axes = plt.subplots(figsize=(12,3), nrows=1, ncols=3)

    for out in range(0,n_outputs):
        integrator.Load("{}/{}_{:05}".format(output_dir, what, out))
        time = integrator.time/wunits.Myr

        hydro = integrator.hydro
        ncells = hydro.ncells
        if time < sn_time:
            ls = '--'
        else:
            ls ='-'

        # density, temperature and ionisation fraction profile
        axes[0].plot(hydro.x[0:ncells]/wunits.pc,hydro.nH[0:ncells],   ls=ls, color=colors[out])
        axes[1].plot(hydro.x[0:ncells]/wunits.pc,hydro.T[0:ncells],    ls=ls, color=colors[out])
        axes[2].plot(hydro.x[0:ncells]/wunits.pc,hydro.xhii[0:ncells], ls=ls, color=colors[out])

    for ax in axes:
        ax.set_xlabel("radius [pc]")
        ax.grid()
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[0].set_ylabel("$n_{\mathrm{H}}$ [cm$^{-3}$]")
    axes[1].set_ylabel("temperature [K]")
    axes[2].set_ylabel('ionisation fraction')

    plt.savefig('{}_time_evolution_{}.png'.format(output_dir[:-1], what), dpi=200, bbox_inches='tight')
    plt.close()

''' Show time evolution of the total momentum '''
def plot_momentum_evolution(simulation_setup):

    output_dir = construct_directory_name(simulation_setup)
    sn_time=load_sn_time(simulation_setup['star mass'])
    times = np.loadtxt(output_dir+'/time.txt')
    momentum = np.loadtxt(output_dir+'/momentum.txt')

    fig, axes = plt.subplots(figsize=(5,5), nrows=1, ncols=1)

    axes.plot(times-sn_time, momentum)

    #axes.set_xlim(0, 0.002)
    axes.set_xlabel("time [Myr]")
    axes.set_ylabel("momentum [g cm s$^{-1}$]")
    axes.grid()
    #axes.set_xscale('symlog', linthresh=1e-3)
    axes.set_yscale('log')
    #axes.set_xlim(min(times-sn_time), max(times-sn_time))

    #if simulation_setup['supernova'] and not simulation_setup['radiation'] and not simulation_setup['winds']:
    #axes.set_xlim(0, max(times-sn_time))

    plt.savefig('{}_momentum.png'.format(output_dir[:-1]), dpi=200, bbox_inches='tight')
    plt.close()

''' Show the time evolution of the size of the bubble '''
def plot_bubble_radius(simulation_setup):

    output_dir = construct_directory_name(simulation_setup)
    sn_time=load_sn_time(simulation_setup['star mass'])
    times = np.loadtxt(output_dir+'/time.txt')

    fig, axes = plt.subplots(figsize=(7,3), nrows=1, ncols=1)

    if simulation_setup['radiation']:
        HII_radius = np.loadtxt(output_dir+'/HII_radius.txt')
        axes.plot(times - sn_time, HII_radius, label='HII bubble')
    if simulation_setup['supernova']:
        sn_radius = np.loadtxt(output_dir+'/sn_radius.txt')
        axes.plot(times - sn_time, sn_radius, label='SN bubble')

    #    axes.scatter(time-sn_time, radiation_bubble_radius, s=5, edgecolor=colors[out], facecolor=colors[out])

    # plot Sedov-Taylor solution
    sedov_times = np.linspace(0, max(times - sn_time), 1000)
    radii_sedov = weltgeist.analyticsolutions.SedovTaylorSolution(sedov_times*wunits.Myr, 1.e51, simulation_setup['central density']*wunits.mH)
    axes.plot(sedov_times, radii_sedov/wunits.pc, label="Sedov-Taylor",color="black")

    #axes.set_xscale('symlog', linthresh=1e-5)
    #axes.set_xscale('log')
    axes.set_xlabel("time [Myr]")
    axes.set_ylabel("shock radius [pc]")
    axes.grid()
    axes.legend()

    if simulation_setup['supernova'] and not (simulation_setup['radiation'] or simulation_setup['winds']):
        axes.set_xlim(0, max(times-sn_time))
    else:
        axes.set_xlim(min(times-sn_time), max(times-sn_time))

    plt.savefig('{}_bubble_radius.png'.format(output_dir[:-1]), dpi=200, bbox_inches='tight')
    plt.close()

def compare_simulations(simulation_setups):
    fig, axes = plt.subplots(figsize=(5,5), nrows=1, ncols=2, sharex=True)
    #colors = [(0.3,0.9,0.6),(0.1,0.8,0.9),(0.1,0.5,1),(0,0.1,0.7)]
    colors = define_colors(len(simulation_setups))

    for simulation_setup, color in zip(simulation_setups, colors):
        output_dir = construct_directory_name(simulation_setup)
        sn_time=load_sn_time(simulation_setup['star mass'])
        times = np.loadtxt(output_dir+'/time.txt')
        # momentum
        momentum = np.loadtxt(output_dir+'/momentum.txt')
        axes[0].plot(times-sn_time, momentum, color=color, label=simulation_setup['ncells'])
        # HII bubble radius
        if simulation_setup['radiation']:
            HII_radius = np.loadtxt(output_dir+'/HII_radius.txt')
            axes[1].plot(times - sn_time, HII_radius, lw=1, color=color)
        # SN bubble
        if simulation_setup['supernova']:
            sn_radius = np.loadtxt(output_dir+'/sn_radius.txt')
            axes[1].plot(times - sn_time, sn_radius, color=color)

    # plot Sedov-Taylor solution
    sedov_times = np.linspace(0, max(times - sn_time), 1000)
    radii_sedov = weltgeist.analyticsolutions.SedovTaylorSolution(sedov_times*wunits.Myr, 1.e51, simulation_setup['central density']*wunits.mH)
    axes[1].plot(sedov_times, radii_sedov/wunits.pc, label="Sedov-Taylor",color="black")

    axes[1].set_xlabel("time [Myr]")
    axes[0].set_ylabel("momentum [g cm s$^{-1}$]")
    axes[1].set_ylabel("shock radius [pc]")
    axes[0].grid()
    axes[1].grid()
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[0].set_ylim(1e40,1e45)
    axes[1].set_ylim(10,500)
    axes[0].legend()

    plt.savefig('comparison.png', dpi=200, bbox_inches='tight')
    plt.close()

def resolution_study():

    setups = []

    for res in [200, 400, 800, 1600,3200]:
        simulation_setup = {'star mass': 40, 'central density': 1, 'density profile slope': 0,
                            'supernova': True, 'radiation': True, 'winds': False,
                            'boxlen':200, 'ncells':res}
        setups.append(simulation_setup)
  
    compare_simulations(setups)


def momentum_profile(simulation_setup, what='output'):
    # gather outputs
    output_dir = construct_directory_name(simulation_setup)
    n_outputs = 0
    for item in os.listdir(output_dir):
        if item.startswith(what):
            n_outputs += 1

    integrator = weltgeist.integrator.Integrator()

    colors = define_colors(n_outputs)
    sn_time=load_sn_time(simulation_setup['star mass'])

    fig, axes = plt.subplots(figsize=(5,5), nrows=1, ncols=1)

    for out in range(0,n_outputs):
        integrator.Load("{}/{}_{:05}".format(output_dir, what, out))
        time = integrator.time/wunits.Myr

        hydro = integrator.hydro
        ncells = hydro.ncells
        if time < sn_time:
            ls = '--'
        else:
            ls ='-'

        momentum = hydro.mass[0:ncells] * integrator.hydro.vel[0:ncells]
        axes.plot(hydro.x[0:ncells]/wunits.pc, momentum, ls=ls, color=colors[out])

    axes.set_xlabel("radius [pc]")
    axes.grid()
    axes.set_yscale('symlog')
    axes.set_ylabel("momentum  g cm s$^{-1}$]")
    #axes.set_xlim(0,100)

    plt.savefig('{}_time_evolution_momentum_{}.png'.format(output_dir[:-1], what), dpi=200, bbox_inches='tight')
    plt.close()

if __name__=="__main__":
    ''' 
    ---- OBJECTIVES ----
    How does early stellar feedback set the stage for supernovae?
    How do the SN properties change when it explodes in a hole carved by HII?
    Does the picture change when we add stellar winds on top of HII?

    ---- OBSERVABLES ----
    * radius of the bubble (shock position)
    * density,temperature, ionisation profiles of the bubble
    * momentum of the SN

    ---- PARAMETERS ----

    KEPT FIXED
    * resolution
    * end time simulation
    * box size
    * gamma
    * SN ejecta mass
    * intial temperature
    '''

    # variable parameters
    simulation_setup = {'star mass': 40,               # Msun                 [10, 20, 40, 80]
                        'central density': 1,          # atoms / cm^-3        [0.01, 0.1, 1, 10, 100]
                        'density profile slope': 0,    # n = n0 * x^-slope    [0] #[0, 1, 1.5, 2]
                        'supernova': True,             #                      [True]
                        'radiation': True,             #                      [True, False]
                        'winds': False,                #                      [True, False]
                        'boxlen':200,                  # size of the box in pc, n0.01 -> 300, n0.1 -> 200, n1 -> 100, n10
                        'ncells':6400                  # fixed to get certain resolution
                        }

    print(simulation_setup)

    run_simulation(simulation_setup)
    plot_time_evolution(simulation_setup, what='output')    
    plot_momentum_evolution(simulation_setup)
    plot_bubble_radius(simulation_setup)
    momentum_profile(simulation_setup, what='output')

    #resolution_study()
