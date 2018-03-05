
import vhone, units

initialised = False

# Physical values to initialise grid to
rmax = 10*units.pc # 10 pc
n0 = 1000.0 # H atoms / cm^-3
rho0 = n0*units.g/units.X # g cm^-3
T0 = 10.0 # K
# P=rho*kB*T
P0 = n0*units.kB*T0 

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
        vhone.data.imax = 2048
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

        vhone.data.gam    = 1.4 # Value from RAMSES

        # Initialise the computational grid
        vhone.data.setup()

        # Now set up the hydro variables for the problem
        #vhone.data.zpr[0,0,0] = 1e7 # simple 1D planar Sedov test

        nx = vhone.data.imax
        vhone.data.zro[0:nx,0,0] = rho0/units.density
        vhone.data.zpr[0:nx,0,0] = P0/units.pressure

        initialised = True

if __name__=="__main__":
    init()
