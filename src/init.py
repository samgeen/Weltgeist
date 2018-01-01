
import vhone

class Fields(object):
    def __init__(self):
        pass
    
    def Density(self):
        nx = vhone.data.imax
        return vhone.data.zro[:,0,0]
    
    def Pressure(self):
        nx = vhone.data.imax
        return vhone.data.zpr[0:nx,0,0]

hydro = Fields()
initialised = False


# Physical quantities (base units in cgs)
pc = 3.086e+18
g = 1.66e-24
year = 3.154e+7
kB = 1.3806485279e-16 # in cgs

# Units
ux = pc
udens = g/ux**3.0

# Physical values to initialise grid to
X = 0.74
rmax = 100.0 # pc
rho0 = 1000.0/X # g cm^-3
T0 = 10.0 # K
# P=rho*kB*T
P0 = rho0*kB*T0

def init():
    global hydro, initialised
    '''
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
    '''

    if not initialised:
    
        # Define the computational grid...
        vhone.data.imax = 256
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
        vhone.data.xmax   = 1.0*rmax
        vhone.data.ymin   = 0.0
        vhone.data.ymax   = 1.0
        vhone.data.zmin   = 0.0
        vhone.data.zmax   = 1.0

        vhone.data.gam    = 1.4 # Value from RAMSES

        # Initialise the computational grid
        vhone.data.setup()

        # Now set up the hydro variables for the problem
        #vhone.data.zpr[0,0,0] = 1e7 # simple 1D planar Sedov test

        nx = vhone.data.imax
        vhone.data.zro[0:nx,0,0] = rho0
        vhone.data.zpr[0:nx,0,0] = P0

if __name__=="__main__":
    init()
