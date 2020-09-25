"""
Hydro solver and grid variables
Sam Geen, January 2018
"""

import numpy as np

from . import vhone, units

class _Field(object):
    """
    Object that allows the user to access the arrays in the underlying
     hydro code. This allows the user to treat fields like a normal
     array without messing up the underlying data structures
    """
    def __init__(self,getter,setter,ArrayFunc):
        """
        Constructor

        Parameters
        ----------

        getter, setter, ArrayFunc: function
            Functions to call to get and set values in the array
            ArrayFunc is used when getting values from the whole array
        """
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
        """
        Returns a string for printing this
        TODO: maybe make more useful?
        """
        return str(self._ArrayFunc())
        

class _Hydro(object):
    """
    Stores and manages the hydrodynamic variables for use by python
    Singleton class (since VH1/F2PY can't handle multiple instances)
    """
    # ---------------
    # HYDRO FIELDS
    # ---------------
    # TODO: MAKE EVERY FIELD A PYTHON @PROPERTY
    # WARNING: This assumes a constant spherical grid
    def __init__(self):
        self.ncells = vhone.data.imax
        # Views to hydro variables in VH-1
        """ 
        BASIC VARIABLES---------------------------------
        These variables represent the basic fluid values
        """
        """ 
        POSITION
        This just returns the grid positions
        Grid positions are (so far) not changeable
        """
        def _xget(slicer):
            return vhone.data.zxa[slicer]*units.distance
        def _xset(slicer,val):
            print("Can't set the grid once the simulation has started!")
            raise ValueError
            # TODO: perhaps explore if this is indeed possible?
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
        
        """ 
        CELL VOLUME
        Calculated from the cell size and position
        self.vol = dx*(x*(x+dx)+dx*dx/3.0) # from volume.f90
        NOTE: x is the *inside* radius, so vol = 4/3*pi*[(r+dr)**3 - r**3]
        """
        self.vol = 4*np.pi*(x**2*dx + x*dx**2 + dx**3/3.0)
        #self.vol = 4.0/3.0 * np.pi * ((x+0.5*dx)**3.0 -(x-0.5*dx)**3.0)
        #self.vol[0] = 4.0/3.0 * np.pi * (0.5*dx)**3.0
        """ 
        DENSITY (rho)
        This is the total mass density in each grid cell
        """
        def _rhoget(slicer):
            return vhone.data.zro[slicer,0,0]*units.density
        def _rhoset(slicer,val):
            vhone.data.zro[slicer,0,0] = val/units.density
        def _rhoarr():
            return vhone.data.zro[0:self.ncells,0,0]*units.density
        self.rho = _Field(_rhoget,_rhoset,_rhoarr)
        """ 
        HYDROGEN NUMBER DENSITY
        This the density rho converted into a Hydrogen number density
        """
        def _nHget(slicer):
            return self.rho[slicer]/units.mH*units.X
        def _nHset(slicer,val):
            vhone.data.zro[slicer,0,0] = val*units.mH/units.X/units.density
        def _nHarr():
            return self.rho[0:self.ncells]/units.mH*units.X
        self.nH = _Field(_nHget,_nHset,_nHarr)
        """ 
        PRESSURE
        This is the *thermal* pressure in each cell
        Note: in some codes this is the total energy in the cell
         (kinetic + thermal + magnetic + etc etc)
         VH1 does not do this, so this is just thermal energy density
        """
        def _Pget(slicer):
            return vhone.data.zpr[slicer,0,0]*units.pressure
        def _Pset(slicer,val):
            vhone.data.zpr[slicer,0,0] = val/units.pressure
        def _Parr():
            return vhone.data.zpr[0:self.ncells,0,0]*units.pressure
        self.P = _Field(_Pget,_Pset,_Parr)
        """ 
        VELOCITY
        This is the velocity of the gas flow in each cell
        """
        def _velget(slicer):
            return vhone.data.zux[slicer,0,0]*units.velocity
        def _velset(slicer,val):
            vhone.data.zux[slicer,0,0] = val/units.velocity
        def _velarr():
            return vhone.data.zux[0:self.ncells,0,0]*units.velocity
        self.vel = _Field(_velget,_velset,_velarr)
        """
        CELL MASS
        Mass in each grid cell
        NOTE: Updating the mass modifies the velocity to preserve
               kinetic energy
              This is done because it is correct for e.g. wind injection
        """
        def _Mget(slicer):
            return self.vol[slicer]*self.rho[slicer]
        def _Mset(slicer,val):
            # Make sure mass injection is elastic!
            oldke = self.KE[slicer]
            self.rho[slicer] = val/self.vol[slicer]
            self.KE[slicer] = oldke
        def _Marr():
            return self.vol*self.rho[0:self.ncells]
        self.mass = _Field(_Mget,_Mset,_Marr)

        """ 
        DERIVED VARIABLES--------------------------------------
        These variables are combinations of the variables above
        Be a bit careful - changing them usually does something
         sensible, but make sure it does what you expect!
        """
        """ 
        TEMPERATURE
        Calculated using the ideal gas equation
        Changing it will alter the pressure & keep density the same
        """
        def _Tget(slicer):
            return self.P[slicer]/self.nH[slicer]/units.kB
        def _Tset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val*self.nH[slicer]*units.kB
            vhone.data.zpr[slicer,0,0] = newP/units.pressure
        def _Tarr():
            return self.P[0:self.ncells]/self.nH[0:self.ncells]/units.kB
        self.T = _Field(_Tget,_Tset,_Tarr)
        """ 
        SOUND SPEED
        Related to temperature
        Changing it is similar to changing the temperature
        """
        def _Csget(slicer):
            return np.sqrt(self.gamma*self.P[slicer]/self.rho[slicer])
        def _Csset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val**2.0 * self.rho[slicer] / self.gamma
            vhone.data.zpr[slicer,0,0] = newP/units.pressure
        def _Csarr():
            return np.sqrt(self.gamma*self.P[0:self.ncells]/self.rho[0:self.ncells])
        self.cs = _Field(_Csget,_Csset,_Csarr)
        """ 
        KINETIC ENERGY
        1/2 m v^2
        Changing it changes the velocity, not the density
        """
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
        """ 
        THERMAL ENERGY
        Uses ideal gas assuming atomic composition (3/2 P V)
        Changing it changes the pressure
        """
        def _TEget(slicer):
            # 3/2 P V
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
        """ 
        GRAVITATIONAL POTENTIAL ENERGY
        Calculates the GPE up to a cell radius from all the mass inside
        Changing it changes the pressure
        """
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
        """ 
        PYTHON-ONLY VARIABLES---------------------------------------
        These variables do not affect the underlying simulation code
        They are used in the Python-only modules (e.g. radiation)
        """
        """ 
        HYDROGEN IONISATION FRACTION xHII
        What fraction of hydrogen is ionised?
        """
        self.xhii = np.zeros(self.ncells)
        """ 
        METALLICITY Z
        What is the mass fraction of metals?
        Given as a fraction of solar (Z=0.014)
        """
        self.Zsolar = np.zeros(self.ncells)+1.0
        """ 
        IONISING PHOTON RATE
        Rate of ionising photons through each cell
        """
        self.Qion = np.zeros(self.ncells)
        """ 
        DUST ABSORPTION CROSS SECTION
        Absorption cross-section of dust in each radial unit
        """
        self.sigmaDust = np.zeros(self.ncells)
        """ 
        GRAVITY
        The acceleration from gravity towards the centre
        NOTE: This is not fully debugged yet
              We need a good test problem to make sure the units are ok
        """
        def _Gget(slicer):
            return vhone.data.zgr[slicer,0,0]*units.gravity
        def _Gset(slicer,val):
            vhone.data.zgr[slicer,0,0] = val/units.gravity
        def _Garr():
            return vhone.data.zgr[0:self.ncells,0,0]*units.gravity
        self.grav = _Field(_Gget,_Gset,_Garr)

    @property
    def gamma(self):   
        """ 
        Adiabatic index of the gas
        
        Returns
        -------

        gamma: float
            adiabatic index
        """
        return vhone.data.gam