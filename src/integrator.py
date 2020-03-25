"""
Engine for running the integrator
Sam Geen, February 2018
"""

import hydro
import cooling, gravity, sources, units, vhone

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
        self._UpdateTime()
        # Cooling step
        if cooling.cooling_on:
            cooling.solve_cooling(self.dt)
        # Inject sources (includes radiation step!)
        sources.InjectSources(self.time, self.dt)
        # Gravity step
        # NOTE: gravity is currently in-testing, use with caution
        if gravity.gravity_on:
            gravity.calculate_gravity()
        else:
            hydro.grav[0:hydro.ncells] = 0.0
        # Hydro step
        vhone.data.step()
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
        vnew = vin / units.velocity / (self.variables.dx / units.distance)
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
        """
        return self._time_code*units.time
    
    @property
    def dt(self):
        """
        Get the last timestep length on the grid
        NOTE: don't use self._dt_code yourself! Use this
        """
        return self._dt_code*units.time

    @property
    def variables(self):
        """
        Return the object for accessing the hydrodynamic variables
        """
        return self._hydro
