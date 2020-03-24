'''
Hydro solver and grid variables
Sam Geen, January 2018
'''

import vhone, init, units
import numpy as np

class _Field(object):
    def __init__(self,getter,setter,ArrayFunc):
        # Set get and set methods
        setattr(self,"_getitem",getter)
        setattr(self,"_setitem",setter)
        # Function that returns the raw array
        self._ArrayFunc = ArrayFunc

    def __getitem__(self,item):
        #print "GET"
        return self._getitem(item)

    def __setitem__(self,item,val):
        #print "SET", val
        self._setitem(item,val)

    def __getattr__(self, name):
        # Do other operation
        #print "GETATTR"
        try:
            return getattr(self._ArrayFunc(), name)
        except AttributeError:
            raise AttributeError(
                    "'Array' object has no attribute {}".format(name))

    def __str__(self):
        #print "STR", str(self._ArrayFunc())
        return str(self._ArrayFunc())
        

class _Hydro(object):
    '''
    Stores and manages the hydro variables for use by python
    Singleton class (since VH1/F2PY is a singleton)
    '''
    # ---------------
    # HYDRO FUNCTIONS
    # ---------------
    # Courant limiter function
    def CourantLimiter(self,vin):
        vnew = vin / units.velocity / (self.dx / units.distance)
        vhone.data.vdtext = max(vhone.data.vdtext,vnew)
    # ---------------
    # HYDRO FIELDS
    # ---------------
    # TODO: MAKE EVERY FIELD A PYTHON @PROPERTY
    # WARNING: This assumes a constant spherical grid
    def __init__(self):
        init.init()
        self.ncells = vhone.data.imax
        # Views to hydro variables in VH-1
        # POSITION
        def _xget(slicer):
            return vhone.data.zxa[slicer]*units.distance
        def _xset(slicer,val):
            print("Can't set the grid once the simulation has started!")
            raise ValueError
            #vhone.data.x[slicer,0,0] = val/units.distance
        def _xarr():
            return vhone.data.zxa[0:self.ncells,0,0]*units.density
        self.x = _Field(_xget,_xset,_xarr)
        x = self.x[0:self.ncells]
        # Assume that cells are evenly spaced in radius
        dx = self.x[1]
        self.dx = dx
        if self.dx != (self.x[2] - self.x[1]):
            print("Grid not evenly spaced!")
            raise ValueError
        # RHO
        def _rhoget(slicer):
            return vhone.data.zro[slicer,0,0]*units.density
        def _rhoset(slicer,val):
            vhone.data.zro[slicer,0,0] = val/units.density
        def _rhoarr():
            return vhone.data.zro[0:self.ncells,0,0]*units.density
        self.rho = _Field(_rhoget,_rhoset,_rhoarr)
        # HYDROGEN NUMBER DENSITY
        def _nHget(slicer):
            return self.rho[slicer]/units.g*units.X
        def _nHset(slicer,val):
            vhone.data.zro[slicer,0,0] = val*units.g/units.X/units.density
        def _nHarr():
            return self.rho[0:self.ncells]/units.g*units.X
        self.nH = _Field(_nHget,_nHset,_nHarr)
        # PRESSURE
        def _Pget(slicer):
            return vhone.data.zpr[slicer,0,0]*units.pressure
        def _Pset(slicer,val):
            vhone.data.zpr[slicer,0,0] = val/units.pressure
        def _Parr():
            return vhone.data.zpr[0:self.ncells,0,0]*units.pressure
        self.P = _Field(_Pget,_Pset,_Parr)
        # VELOCITY
        def _velget(slicer):
            return vhone.data.zux[slicer,0,0]*units.velocity
        def _velset(slicer,val):
            vhone.data.zux[slicer,0,0] = val/units.velocity
        def _velarr():
            return vhone.data.zux[0:self.ncells,0,0]*units.velocity
        self.vel = _Field(_velget,_velset,_velarr)
        # Derived variables
        # TEMPERATURE
        def _Tget(slicer):
            return self.P[slicer]/self.nH[slicer]/units.kB
        def _Tset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val*self.nH[slicer]*units.kB
            vhone.data.zpr[slicer,0,0] = newP/units.pressure
        def _Tarr():
            return self.P[0:self.ncells]/self.nH[0:self.ncells]/units.kB
        self.T = _Field(_Tget,_Tset,_Tarr)
        # SOUND SPEED
        def _Csget(slicer):
            return np.sqrt(init.gamma*self.P[slicer]/self.rho[slicer])
        def _Csset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val**2.0 * self.rho[slicer] / init.gamma
            vhone.data.zpr[slicer,0,0] = newP/units.pressure
        def _Csarr():
            return np.sqrt(init.gamma*self.P[0:self.ncells]/self.rho[0:self.ncells])
        self.cs = _Field(_Csget,_Csset,_Csarr)
        # KINETIC ENERGY
        # NOTE: In VH1, zpr is SOLELY thermal pressure!!! (unlike RAMSES)
        def _KEget(slicer):
            # 0.5 * mass * v^2
            #print "KEGET"
            return 0.5*self.mass[slicer]*self.vel[slicer]**2.0
        def _KEset(slicer,val):
            # dv = sqrt(2.0 * dKE / mass)
            #print "KESET"
            vhone.data.zux[slicer,0,0] = np.sqrt(2.0*val/self.mass[slicer])/units.velocity
        def _KEarr():
            # 0.5 * mass * v^2
            #print "KEARR"
            return 0.5*self.mass[0:self.ncells]*self.vel[0:self.ncells]**2.0
        self.KE = _Field(_KEget,_KEset,_KEarr)
        # THERMAL ENERGY
        def _TEget(slicer):
            # 0.5 * mass * v^2
            #print "KEGET"
            return 1.5 * self.P[slicer] * self.vol[slicer]
        def _TEset(slicer,val):
            # dv = sqrt(2.0 * dKE / mass)
            #print "KESET"
            vhone.data.zpr[slicer,0,0] = val/(1.5*self.vol[slicer]*units.pressure)
        def _TEarr():
            # 0.5 * mass * v^2
            #print "KEARR"
            return 1.5 * self.P[0:self.ncells] * self.vol[0:self.ncells]
        self.TE = _Field(_TEget,_TEset,_TEarr)
        # GRAVITATIONAL POTENTIAL ENERGY
        def _GPEget(slicer):
            # G*M*m/r
            Minside = np.cumsum(self.mass[0:self.ncells])
            invx = Minside*0.0
            invx[1:] = 1.0/self.x[1:]
            return units.G * self.mass[slicer] * Minside[slicer] * invx[slicer]
        def _GPEset(slicer,val):
            # Not implemented yet!
            print("GPE setting not implemented! (What would you set, though?)")
            raise NotImplementedError
        def _GPEarr():
            # G*M*m/r
            return _GPEget(slice(0,self.ncells))
        self.GPE = _Field(_GPEget,_GPEset,_GPEarr)
        # Variables contained only in this module
        # XHII (Hydrogen ionisation fraction)
        self.xhii = np.zeros(self.ncells)
        # Metallicity in solar units
        self.Zsolar = np.zeros(self.ncells)+1.0
        # GRAVITY
        def _Gget(slicer):
            return vhone.data.zgr[slicer,0,0]*units.gravity
        def _Gset(slicer,val):
            vhone.data.zgr[slicer,0,0] = val/units.gravity
        def _Garr():
            return vhone.data.zgr[0:self.ncells,0,0]*units.gravity
        self.grav = _Field(_Gget,_Gset,_Garr)
        # HARD-CODED EXTERNAL PRESSURE FIELD (TODO)
        self.Pext = np.zeros(self.ncells)
        # VOLUME
        #self.vol = dx*(x*(x+dx)+dx*dx/3.0) # from volume.f90
        # x is the *inside* radius, so vol = 4/3*pi*[(r+dr)**3 - r**3]
        self.vol = 4*np.pi*(x**2*dx + x*dx**2 + dx**3/3.0)
        #self.vol = 4.0/3.0 * np.pi * ((x+0.5*dx)**3.0 -(x-0.5*dx)**3.0)
        #self.vol[0] = 4.0/3.0 * np.pi * (0.5*dx)**3.0
        # MASS
        def _Mget(slicer):
            return self.vol[slicer]*self.rho[slicer]
        def _Mset(slicer,val):
            # Make sure mass injection is elastic!
            oldke = self.KE[slicer]
            self.rho[slicer] = val/self.vol[slicer]
            self.KE[slicer] = oldke
        def _Marr():
            return self.vol*self.rho[0:ncells]
        self.mass = _Field(_Mget,_Mset,_Marr)


hydro = _Hydro()

# Test functions
if __name__=="__main__":
    print(str(hydro.mass))