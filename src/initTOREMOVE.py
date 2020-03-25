"""
Initialises the simulation
Sets up the grid and hydrodynamic quantities on it

Sam Geen, January 2018
"""

import vhone, units

# Is the grid initialised?
# Ensures that this isn't done twice, which would cause issues
initialised = False

def init(rmaxIn = 20.0*units.pc,
        n0 = 1000.0, # H atoms / cm^-3
        T0In = 10.0, # K
        gammaIn = 5.0/3.0):
    """
    Main initialisation function
    Note that the grid can be altered at any time using the hydro module

    Parameters
    ----------

    rmaxIn: float
        Sets the maximum size of the grid
        default: 20 pc
    n0In: float
        Sets the initial Hydrogen number density of the gas
        default: 1000 cm^-3
    T0In: float
        Sets the initial temperature of the gas
        default: 10 K
    gammaIn: float
        Sets the adiabatic index of the gas
        default: 5/3 (monatomic)
    """

    global initialised

    # Derived quantities
    rho0 = n0*units.mp # g cm^-3
    P0 = n0*units.kB*T0 

    if not initialised:
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
        vhone.data.imax = 512
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

        # Set initialised to true, to prevent this running again
        # TODO: check whether this actually results in an error? 
        #       might be useful...
        initialised = True

if __name__=="__main__":
    init()
