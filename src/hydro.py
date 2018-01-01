'''
Hydro solver and grid variables
Sam Geen, January 2018
'''

import vhone, init
import numpy as np

class _Field(object):
    def __init__(self,getter,setter):
        # Set get and set methods
        setattr(self,"_getitem",getter)
        setattr(self,"_setitem",setter)

    def __getitem__(self,item):
        return self._getitem(item)

    def __setitem__(self,item,val):
        self._setitem(item,val)
        

class _Hydro(object):
    # init init init init
    def __init__(self):
        init.init()
        self.ncells = vhone.data.imax
        # Views to hydro variables in VH-1
        # RHO
        def _rhoget(slicer):
            return vhone.data.zro[slicer,0,0]
        def _rhoset(slicer,val):
            vhone.data.zro[slicer,0,0] = val
        self.rho = _Field(_rhoget,_rhoset)
        # PRESSURE
        def _Pget(slicer):
            return vhone.data.zpr[slicer,0,0]
        def _Pset(slicer,val):
            vhone.data.zpr[slicer,0,0] = val
        self.P = _Field(_Pget,_Pset)
        # VELOCITY (TODO: Check it's actually zuz)
        def _velget(slicer):
            return vhone.data.zuz[slicer,0,0]
        def _velset(slicer,val):
            vhone.data.zuz[slicer,0,0] = val
        self.vel = _Field(_velget,_velset)
        # Derived variables
        # TEMPERATURE
        def _Tget(slicer):
            return vhone.data.zpr[slicer,0,0]/vhone.data.zro[slicer,0,0]/init.kB
        def _Tset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val*vhone.data.zro[slicer,0,0]*init.kB
            vhone.data.zpr[slicer,0,0] = newP
        self.T = _Field(_Tget,_Tset)
        # Variables contained only in this module
        # XHII (Hydrogen ionisation fraction)
        self.xhii = np.zeros(ncells)
        # GRAV (TODO)
        self.grav = np.zeros(ncells)
        # HARD-CODED EXTERNAL PRESSURE FIELD (TODO)
        self.Pext = np.zeros(ncells)

hydro = _Hydro()



# Test functions
if __name__=="__main__":
    hydro.rho[2] = 4.0
    hydro.rho[3:5] = 16.0
    print hydro.rho[0:20]