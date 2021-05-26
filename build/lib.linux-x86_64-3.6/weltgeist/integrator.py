"""
Engine for running the integrator
Sam Geen, February 2018
"""

import h5py
import numpy as np

from . import cooling, hydro, gravity, sources, units, vhone

# Instance the integrator, using singleton pattern
_integrator = None
def Integrator():
    """
    Method that returns the singleton _Integrator object
    The reason for doing this is that f2py doesn't like instancing
     multiple versions of modules run with it
    """
    global _integrator
    if _integrator is None:
        _integrator = _Integrator()
    return _integrator

class _Integrator(object):
    """
    This class controls the operation of the underlying simulation code
    """
    def __init__(self):
        """
        Constructor
        """
        self._initialised = False
        # Internal time values
        self._time_code = 0.0
        self._dt_code = 0.0
        # Hydro variables
        self._hydro = None

    def Save(self,filename):
        '''
        Save the current simulation state to file
        TODO: Add sources
        '''
        file = h5py.File(filename+".hdf5", "w")
        # Save a file format version 
        file.attrs['version']="1.0.0"
        # Save setup parameters
        hydro = self.hydro
        ncells = hydro.ncells
        file.create_dataset("ncells", data=(ncells,), dtype=np.int32)
        rmax = vhone.data.xmax * units.distance
        file.create_dataset("rmax", data=(rmax,), dtype=np.float64)
        # (Note: we don't save n0 and T0 because these overwritten by the grid state)
        file.create_dataset("gamma",data=(hydro.gamma,),dtype=np.float64)
        # Save the time variables
        file.create_dataset("time",data=(vhone.data.time,),dtype=np.float64)
        file.create_dataset("dt",data=(self._dt_code,),dtype=np.float64)
        # Save courant limiter
        file.create_dataset("vcourant",data=(vhone.data.vdtext,),dtype=np.float64)
        # Save the basic hydro variables
        file.create_dataset("rho",data=hydro.rho[0:ncells],dtype=np.float64)
        file.create_dataset("P",data=hydro.P[0:ncells],dtype=np.float64)
        file.create_dataset("vel",data=hydro.vel[0:ncells],dtype=np.float64)
        # Save python-only hydro variables
        file.create_dataset("xhii",data=hydro.xhii[0:ncells],dtype=np.float64)
        file.create_dataset("Zsolar",data=hydro.Zsolar[0:ncells],dtype=np.float64)
        file.create_dataset("grav",data=hydro.grav[0:ncells],dtype=np.float64)
        # Save switches in the modules
        file.create_dataset("switches",data=(cooling.cooling_on,gravity.gravity_on),dtype=np.float64)
        # TODO: Save sources (this is the hard one...)
        # We probably have to have serialisation options inside sources
        # ...
        # Done!
        file.close()

    def Load(self,filename):
        '''
        Load the current simulation state from file
        TODO: Add sources
        '''
        # Check for an already initialised VH1
        if self._initialised:
            print("VH-1 already initialised, will try to load anyway...")
        # Open file
        file = h5py.File(filename+".hdf5", "r")
        # Save a file format version 
        #file.attrs['version']="1.0.0"
        # Save setup parameters
        def loaditem(varname):
            data = np.array(file.get(varname))
            if len(data) == 1:
                data = data[0]
            return data
        ncells = loaditem("ncells")
        rmax = loaditem("rmax")
        gamma = loaditem("gamma")
        # Do a check that the loaded values don't clash with the setup values
        if self._initialised:
            hydro = self.hydro
            if hydro.ncells == ncells:
                print ("Error loading: ncells is incorrect (try not calling integrator.Setup)")
                raise ValueError
            if vhone.data.xmax == rmax / units.distance:
                print ("Error loading: rmax is incorrect (try not calling integrator.Setup)")
                raise ValueError
            if gamma != hydro.gamma:
                print ("Error loading: gamma is incorrect (try not calling integrator.Setup)")
                raise ValueError
        else:
            # Note: n, T will change anyway
            self.Setup(ncells = ncells,
                rmax = rmax,
                n0 = 1000.0, # H atoms / cm^-3
                T0 = 10.0, # K
                gamma = gamma)
            hydro = self.hydro
        # Update time
        time = loaditem("time")
        dt = loaditem("dt")
        self._time_code = time
        vhone.data.time = time
        self._dt_code = dt
        vhone.data.time = time
        # Update the courant limiter
        vhone.data.vdtext = loaditem("vcourant")
        # Update the hydro variables
        hydro.rho[0:ncells] = loaditem("rho")
        hydro.P[0:ncells] = loaditem("P")
        hydro.vel[0:ncells] = loaditem("vel")
        # Update the python-only hydro variables
        hydro.xhii[0:ncells] = loaditem("xhii")
        hydro.Zsolar[0:ncells] = loaditem("Zsolar")
        hydro.grav[0:ncells] = loaditem("grav")
        # Save switches in the modules
        switches = loaditem("switches")
        cooling.cooling_on = bool(switches[0])
        gravity.gravity_on = bool(switches[1])
        # TODO: Load sources (this is the hard one...)
        # We probably have to have serialisation options inside sources
        # ...
        # Done!
        file.close()

    def Setup(self,
            ncells = 512,
            rmax = 20.0*units.pc,
            n0 = 1000.0, # H atoms / cm^-3
            T0 = 10.0, # K
            gamma = 5.0/3.0):
        """
        Main initialisation function
        Note that the grid can be altered at any time using the hydro module

        Parameters
        ----------

        ncells: integer
            Number of cells to use
            default: 512
        rmax: float
            Sets the maximum size of the grid
            default: 20 pc
        n0: float
            Sets the initial Hydrogen number density of the gas
            default: 1000 cm^-3
        T0: float
            Sets the initial temperature of the gas
            default: 10 K
        gamma: float
            Sets the adiabatic index of the gas
            default: 5/3 (monatomic)
        """

        # Derived quantities
        rho0 = n0*units.mp # g cm^-3
        P0 = n0*units.kB*T0 # ergs cm^-3

        # Running twice is probably an error...?
        if not self._initialised:
            """
            Reference from VH1:
            ! Set up geometry and boundary conditions of grid
            !
            ! Boundary condition flags : nleft, nright
            !   = 0  :  reflecting boundary condition
            !   = 1  :  inflow/outflow boundary condition (zero gradients)
            !   = 2  :  fixed inflow boundary condition (values set by dinflo, pinflo, etc.)
            !   = 3  :  periodic (nmax+1 = nmin; nmin-1 = nmax)
            !
            ! Geometry flag : ngeom                         |  Cartesian:
            !   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0   (x,y,z)
            !   = 1  :  cylindrical radial                  |  Cylindrical:
            !   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0   (s,phi,z)
            !   = 3  :  cylindrical angle                   |
            !   = 4  :  spherical polar angle (theta)       |  Spherical:
            !   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5   (r,theta,phi)
            """

            # Define the computational grid...
            vhone.data.imax = ncells
            vhone.data.jmax = 1
            vhone.data.kmax = 1
            
            vhone.data.ngeomx = 2
            vhone.data.ngeomy = 4
            vhone.data.ngeomz = 5
            
            vhone.data.nleftx = 0
            vhone.data.nrightx= 0
            vhone.data.nlefty = 0
            vhone.data.nrighty= 0
            vhone.data.nleftz = 0
            vhone.data.nrightz= 0
            
            vhone.data.xmin   = 0.0
            vhone.data.xmax   = rmax/units.distance
            vhone.data.ymin   = 0.0
            vhone.data.ymax   = 1.0
            vhone.data.zmin   = 0.0
            vhone.data.zmax   = 1.0

            #vhone.data.gam    = 1.4 # Diatomic, Value from RAMSES
            vhone.data.gam    = gamma # Monatomic = 1.66

            # Initialise the computational grid
            vhone.data.setup()

            # Now set up the hydro variables for the problem
            nx = vhone.data.imax
            vhone.data.zro[0:nx,0,0] = rho0/units.density
            vhone.data.zpr[0:nx,0,0] = P0/units.pressure

            # Initialise hydro object for accessing variables
            self._hydro = hydro._Hydro()

            # Set initialised to true, to prevent this running again
            # TODO: check whether this actually results in an error? 
            #       might be useful...
            self._initialised = True

    def Step(self):
        """
        Run a single hydrodynamic step
        """
        hydro = self._hydro
        # Check if the grid has been initialised
        if not self._initialised:
            print("Error: grid not initialised! Run integrator.Init()")
            raise RuntimeError
        # Update time
        #self._UpdateTime()
        # Gravity step
        # NOTE: gravity is currently in-testing, use with caution
        if gravity.gravity_on:
            gravity.calculate_gravity()
        else:
            hydro.grav[0:hydro.ncells] = 0.0
        # Cooling step
        if cooling.cooling_on:
            cooling.solve_cooling(self.dt)
        # Inject sources (includes radiation step!)
        sources.Sources().InjectSources()
        # Hydro step
        vhone.data.step()
        # Update time to make sure the code sees the correct time
        self._UpdateTime()
        # Final sanity check
        if cooling.cooling_on:
            cooling.CheckTemperature()


    # ---------------
    # HYDRO FUNCTIONS
    # ---------------
    # Courant limiter function
    def CourantLimiter(self,vin):
        """
        Limits the timestep to prevent flows faster than vin
        VH1 has this already, so in general this isn't really needed
        This is mainly useful for winds where we need to prepare the
         grid to accept fast flows and prevent a too-long first 
         timestep

        Parameters
        ----------
        vin: float
            velocity to limit the timestep to 
        """
        vnew = vin / units.velocity / (self.hydro.dx / units.distance)
        vhone.data.vdtext = max(vhone.data.vdtext,vnew)

    def _UpdateTime(self):
        """
        Update the time values
        """
        # Get dt from last time
        oldtime = self._time_code
        # The + 0.0 is a fun Python trick to ensure copy by value 
        #   and not reference
        self._time_code = vhone.data.time + 0.0
        self._dt_code = self._time_code - oldtime

    @property
    def time(self):
        """
        Get the time on the grid
        NOTE: don't use self._time_code yourself! Use this

        Returns
        -------
        dt : float
            current simulation time in seconds
        """
        return self._time_code*units.time
    
    @property
    def dt(self):
        """
        Get the last timestep length on the grid
        NOTE: don't use self._dt_code yourself! Use this

        Returns
        -------
        dt : float
            timestep in seconds
        """
        return self._dt_code*units.time

    @property
    def hydro(self):
        """
        Return the object containing the grid of hydrodynamic variables

        Returns
        -------
        hydro : _Hydro
            object containing grid of hydrodynamic variables
        """
        return self._hydro
