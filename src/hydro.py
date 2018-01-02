'''
Hydro solver and grid variables
Sam Geen, January 2018
'''

import vhone, init
import numpy as np

class _Field(object):
    def __init__(self,getter,setter,ArrayFunc):
        # Set get and set methods
        setattr(self,"_getitem",getter)
        setattr(self,"_setitem",setter)
        # Function that returns the raw array
        self._ArrayFunc = ArrayFunc

    def __getitem__(self,item):
        return self._getitem(item)

    def __setitem__(self,item,val):
        self._setitem(item,val)

    def __getattr__(self, name):
        # Do other operation
        try:
            return getattr(self._ArrayFunc(), name)
        except AttributeError:
            raise AttributeError(
                    "'Array' object has no attribute {}".format(name))

    def __str__(self):
        return str(self._ArrayFunc())
        

class _Hydro(object):
    '''
    Stores and manages the hydro variables for use by python
    Singleton class (since VH1/F2PY is a singleton)
    '''
    # TODO: MAKE EVERY FIELD A PYTHON @PROPERTY
    # WARNING: This assumes a constant spherical grid
    def __init__(self):
        init.init()
        self._time = vhone.data.time
        self._dt = 0.0
        self.ncells = vhone.data.imax
        # Views to hydro variables in VH-1
        # POSITION
        x = vhone.data.zxa
        self.x = x
        # Assume that cells are evenly spaced in radius
        dx = self.x[1]
        self.dx = dx
        if self.dx != (self.x[2] - self.x[1]):
            print "Grid not evenly spaced!"
            raise ValueError
        # RHO
        def _rhoget(slicer):
            return vhone.data.zro[slicer,0,0]
        def _rhoset(slicer,val):
            vhone.data.zro[slicer,0,0] = val
        def _rhoarr():
            return vhone.data.zro[0:self.ncells,0,0]
        self.rho = _Field(_rhoget,_rhoset,_rhoarr)
        # PRESSURE
        def _Pget(slicer):
            return vhone.data.zpr[slicer,0,0]
        def _Pset(slicer,val):
            vhone.data.zpr[slicer,0,0] = val
        def _Parr():
            return vhone.data.zpr[0:self.ncells,0,0]
        self.P = _Field(_Pget,_Pset,_Parr)
        # VELOCITY
        def _velget(slicer):
            return vhone.data.zux[slicer,0,0]
        def _velset(slicer,val):
            vhone.data.zux[slicer,0,0] = val
        def _velarr():
            return vhone.data.zux[0:self.ncells,0,0]
        self.vel = _Field(_velget,_velset,_velarr)
        # Derived variables
        # TEMPERATURE
        def _Tget(slicer):
            return vhone.data.zpr[slicer,0,0]/vhone.data.zro[slicer,0,0]/init.kB
        def _Tset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val*vhone.data.zro[slicer,0,0]*init.kB
            vhone.data.zpr[slicer,0,0] = newP
        def _Tarr():
            return vhone.data.zpr[0:self.ncells,0,0]/ \
                vhone.data.zro[0:self.ncells,0,0]/init.kB
        self.T = _Field(_Tget,_Tset,_Tarr)
        # Variables contained only in this module
        # XHII (Hydrogen ionisation fraction)
        self.xhii = np.zeros(self.ncells)
        # GRAV (TODO)
        self.grav = np.zeros(self.ncells)
        # HARD-CODED EXTERNAL PRESSURE FIELD (TODO)
        self.Pext = np.zeros(self.ncells)
        # VOLUME
        self.vol = dx*(x*(x+dx)+dx*dx/3.0) # from volume.f90
        # MASS
        def _Mget(slicer):
            return self.vol[slicer]*vhone.data.zro[slicer,0,0]
        def _Mset(slicer,val):
            vhone.data.zro[slicer,0,0] = val/self.vol[slicer]
        def _Marr():
            return self.vol*vhone.data.zro[0:self.ncells,0,0]
        self.mass = _Field(_Mget,_Mset,_Marr)

    def Step(self):
        oldtime = self._time
        vhone.data.step()
        self._time = vhone.data.time
        self._dt = self._time - oldtime

    @property
    def time(self):
        return self._time
    
    @property
    def dt(self):
        return self._dt


hydro = _Hydro()

# Test functions
if __name__=="__main__":
    print hydro.mass