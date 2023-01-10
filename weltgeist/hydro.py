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
    Note that this class does not store the data itself, which can be in a Fortran module
    Instead, it is a set of utility functions designed to facilitate safe, "Pythonic" access
    """
    def __init__(self,docstring=None,getitem=None,setitem=None):
        """
        Constructor

        Contains option to also set empty field and set properly once other data has been gathered

        Parameters
        ----------

        docstring: string
            Explanation of what the field does, units, etc
        getter, setter: function
            Functions to call to get and set values in the array
        """
        # Docstring for field
        self._docstring = docstring
        # Set get and set methods
        self._getitem = getitem
        self._setitem = setitem
        # Set docstring of this class to the docstring explaining the underlying variable
        self.__doc__ = docstring

    @property
    def docstring(self):
        """
        Get the docstring of the field
        """
        return self._docstring

    def __getitem__(self,item):
        return self._getitem(item)

    def __setitem__(self,item,val):
        self._setitem(item,val)

    def _assigngetset(self,docstring,getitem,setitem):
        """
        Assign the getter and setter here
        This is a quasi-private function used by hydro to set up the field
        This is an option to set up the field once other data has been set

        Parameters
        ----------

        docstring: string
            Explanation of what the field does, units, etc
        getter, setter: function
            Functions to call to get and set values in the array
        """
        # Docstring for field
        self._docstring = docstring
        # Set get and set methods
        self._getitem = getitem
        self._setitem = setitem
        # Set docstring of this class to the docstring explaining the underlying variable
        self.__doc__ = docstring

    def all(self):
        """
        Return all cells
        """
        return self[:]

    def __str__(self):
        """
        Returns a string for printing this
        TODO: maybe make more useful?
        """
        return str(self.all())

    """
    OPERATORS
    Allow Field objects to be operated on like numpy arrays
    """
    # Logical
    def __lt__(self,b):
        return self[:] < b
    def __le__(self,b):
        return self[:] <= b
    def __eq__(self,b):
        return self[:] == b
    def __ne__(self,b):
        return self[:] != b
    def __ge__(self,b):
        return self[:] <= b
    def __gt__(self,b):
        return self[:] > b
    # Mathematical/bitwise
    def __abs__(self):
        return abs(self[:])
    def __add__(self,b):
        try:
            return self[:] + b[:]
        except:
            return self[:] + b
    def __radd__(self,b):
        try:
            return self[:] + b[:]
        except:
            return self[:] + b
    def __floordiv__(self,b):
        return self[:] // b
    def __mod__(self,b):
        return self[:] % b
    def __mul__(self,b):
        try:
            return self[:] * b[:]
        except:
            return self[:] * b
    def __rmul__(self,b):
        try:
            return self[:] * b[:]
        except:
            return self[:] * b
    def __neg__(self):
        return -self[:]
    def __pos__(self):
        return +self[:]
    def __pow__(self,b):
        return self[:]**b
    def __sub__(self,b):
        try:
            return self[:] - b[:]
        except:
            return self[:] - b
    def __rsub__(self,b):
        try:
            return b[:] - self[:]
        except:
            return b - self[:]
    def __truediv__(self,b):
        return self[:] / b
    def __concat__(self,b):
        return self[:] + b[:]
    def __contains__(self,b):
        return b in self[:]
    # Allow numpy to process the data in here
    def __array__(self) -> np.ndarray:
        return self[:]
    # Allow len(field)
    def __len__(self):
        return len(self[:])
    # In-place operations are handled already by the getter and setter in hydro

    



_internalvariables = {}
_fieldvariables = {}

class _Hydro(object):
    """
    Stores and manages the hydrodynamic variables for use by python
    Singleton class (since VH1/F2PY can't handle multiple instances)

    Contains arrays that describe the physical state of each cell in the grid
    Each array describes a given physical quantity, e.g. gas density, pressure, velocity, temperature
    Quantities can be derived and should be set up to provide internal consistency via the base units, 
      e.g. setting temperature changes the gas pressure
    Each array will be in cgs units for its public-facing workings
    Arrays can refer to Fortran arrays in the integrator or be entirely pythonic
    
    Each array should implement the following:
    1. Get the whole array
        a = hydro.foo   # sets a to the whole array in foo
    2. Set the whole array
        hydro.foo = a   # sets the whole array in foo to a
    3. Get a slice
        a = hydro.foo[b:c]   # gets the array from indices b to c
    4. Set a slice
        hydro.foo[b:c] = a   # sets slice from indices b to c to value a

    This needs the following underlying logic:
    Firstly, each _Field object contains a getitem and setitem function accepting indexes or slices
    Then, for each case given above:
    1. hydro returns the _Field object foo via a getter
    2. hydro calls the _Field object foo's __setitem__ command with arguments (slice(0:self.ncells) , a)
    3. hydro returns the _Field object foo via a getter, which then gets its __getitem__ command called
    4. hydro returns the _Field object foo via a getter, which then gets its __setitem__ command called
    """
    global _internalvariables
    global _fieldvariables
    # ---------------
    # HYDRO FIELDS
    # ---------------
    # TODO: MAKE EVERY FIELD A PYTHON @PROPERTY
    # WARNING: This assumes a constant spherical grid
    def __init__(self):
        global _internalvariables
        global _fieldvariables
        self.ncells = vhone.data.imax

        # Set up variables that need to be made once the data is initialised
        # Cells should be evenly spaced in radius, where self.x[0] is 0.0, if not throw value error
        # Grid spacing
        x = self.x[:]
        self._dx = x[1]
        if self._dx != (x[2] - x[1]):
            print("Grid not evenly spaced!")
            raise ValueError

        # Cell volume
        # vol = dx*(x*(x+dx)+dx*dx/3.0) # from volume.f90
        # NOTE: x is the *inside* radius, so vol = 4/3*pi*[(r+dr)**3 - r**3]
        self._vol = 4*np.pi*(x**2*self.dx + x*self.dx**2 + self.dx**3/3.0)

        # Set up variables stored here (needs to be done once ncells is set, hence this indirect approach)
        for key in _internalvariables.keys():
            _internalvariables[key] = np.zeros(self.ncells)
        _internalvariables["Zsolar"] += 1

        # Set up derived variables / variables using other Fields
        """ 
        THERMAL PRESSURE
        Remove magnetic pressure from the gas and return thermal pressure
        """
        def _PThget(slicer):
            return self.P[slicer] - self.PMagnetic[slicer]
        def _PThset(slicer,val):
            val2 = val + self.PMagnetic[slicer]
            self.P[slicer] = val2
        _PThermalstring = "Gas thermal pressure in erg/cm^3 (note: does not include magnetic or kinetic energy)"
        self._PThermal._assigngetset(_PThermalstring,_PThget,_PThset)

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
        _massstring = "Gas mass in each cell in g"
        self._mass._assigngetset(_massstring,_Mget,_Mset)

        """ 
        TEMPERATURE
        Calculated using the ideal gas equation
        Changing it will alter the pressure & keep density the same
        """
        def _Tget(slicer):
            return self.PThermal[slicer]/self.nH[slicer]/units.kB
        def _Tset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val*self.nH[slicer]*units.kB
            self.PThermal[slicer] = newP
        _Tstring = "Gas temperature in K"
        self._T._assigngetset(_Tstring,_Tget,_Tset)

        """ 
        SOUND SPEED
        Related to temperature
        Changing it is similar to changing the temperature
        """
        def _Csget(slicer):
            return np.sqrt(self.gamma*self.PThermal[slicer]/self.rho[slicer])
        def _Csset(slicer,val):
            # Set the pressure from the ideal gas equation
            newP = val**2.0 * self.rho[slicer] / self.gamma
            self.PThermal[slicer] = newP
        _csstring = "Gas temperature in K"
        self._cs._assigngetset(_csstring,_Csget,_Csset)

        """ 
        KINETIC ENERGY
        1/2 m v^2
        Changing it changes the velocity, not the density
        """
        def _KEget(slicer):
            return 0.5*self.mass[slicer]*self.vel[slicer]**2.0
        def _KEset(slicer,val):
            vhone.data.zux[slicer,0,0] = np.sqrt(2.0*val/self.mass[slicer])/units.velocity
        _KEstring = "Gas kinetic energy in erg"
        self._KE._assigngetset(_KEstring,_KEget,_KEset)

        """ 
        THERMAL ENERGY
        Uses ideal gas assuming atomic composition (3/2 P V)
        Changing it changes the pressure
        """
        def _TEget(slicer):
            # 3/2 P V
            return 1.5 * self.PThermal[slicer] * self.vol[slicer]
        def _TEset(slicer,val):
            self.PThermal[slicer] = val/(1.5*self.vol[slicer])
        _TEstring = "Gas thermal energy in erg"
        self._TE._assigngetset(_TEstring,_TEget,_TEset)

        """ 
        GRAVITATIONAL POTENTIAL ENERGY
        Calculates the GPE up to a cell radius from all the mass inside
        """
        def _GPEget(slicer):
            # G*M*m/r
            Minside = np.cumsum(self.mass[0:self.ncells])
            invx = Minside*0.0
            invx[1:] = 1.0/self.x[1:]
            return units.G * self.mass[slicer] * Minside[slicer] * invx[slicer]
        def _GPEset(slicer,val):
            # Not implemented yet!
            print("GPE setting not implemented! (What would you set, though? Mass?)")
            raise NotImplementedError
        # TODO: Implement source mass as a parameter
        _GPEstring = "Gas gravitational potential energy in erg (NOTE! DOES NOT IMPLEMENT ANY SOURCE MASS YET)"
        self._GPE._assigngetset(_GPEstring,_GPEget,_GPEset)

        """ 
        MAGNETIC PRESSURE
        Magnetic pressure stored below
        """
        # Magnetic pressure field objects
        def _PMagneticget(slicer):
            return _internalvariables["PMagnetic"][slicer]
        def _PMagneticset(slicer,val):
            # Modify both the global pressure and the magnetic pressure tracker
            PBdiff = val - _internalvariables["PMagnetic"][slicer]
            _internalvariables["PMagnetic"][slicer] += PBdiff
            self.P[slicer] += PBdiff
        _PMagneticstring = "Gas magnetic pressure in erg/cm^3"
        # NOTE: The magnetic pressure is stored as a numpy array below
        self._PMagneticField._assigngetset(_PMagneticstring,_PMagneticget,_PMagneticset)

        """ 
        MAGNETIC FIELD
        Magnetic field
        """
        EIGHTPI = 8.0*np.pi
        def _Bget(slicer):
            return np.sqrt(EIGHTPI * self.PMagnetic[slicer])
        def _Bset(slicer,val):
            # Modify both the global pressure and the magnetic pressure tracker
            Pmag = val**2 / EIGHTPI
            self.PMagnetic[slicer] = Pmag
        _Bstring = "Gas magnetic field in Gauss - note: this is calculated semi-analytically per cell"
        self._Bfield._assigngetset(_Bstring,_Bget,_Bset)


    # Views to hydro variables in VH-1
    """
    Utility functions for setting up fields as class properties
    """
    def propertyargs(field):
        def getfield(field):
            return lambda self: field
        def setfield(field):
            def setfunc(self,value):
                field[:] = value
            return setfunc
        return getfield(field),setfield(field),None,field.docstring

    def PythonField(fieldname,docstring):
        """
        Basic field that uses a Python array as its field
        """
        _internalvariables[fieldname] = None
        def _get(slicer):
            return _internalvariables[fieldname][slicer]
        def _set(slicer,val):
            # Modify both the global pressure and the magnetic pressure tracker
            _internalvariables[fieldname][slicer] = val
        return _Field(docstring,_get,_set)

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
        print("Error: You can't set the grid coordinates by hand - instead run integrator.Reset and integrator.Setup again")
        raise ValueError
    _xstring = "Position of each cell in cm - should be uniformly spaced with 0th position at 0.0 cm"
    _x = _Field(_xstring,_xget,_xset)
    x = property(*propertyargs(_x))

    """ 
    GRID SPACING
    This returns a scalar for the grid cells
    Doesn't need a whole field variable to itself, can't be set
    """
    _dx = None # Define in __init__
    @property
    def dx(self):
        """
        Grid spacing in cm - grid is uniform so this is the same for all cells
        """
        return self._dx
    
    """ 
    CELL VOLUME
    Calculated from the cell size and position
    Doesn't need a whole field variable to itself, can't be set
    """    
    @property
    def vol(self):
        """
        Cell volume in cm^3, calculated via dx*(x*(x+dx)+dx*dx/3.0)
        """
        return self._vol

    """ 
    DENSITY (rho)
    This is the total mass density in each grid cell
    """
    def _rhoget(slicer):
        return vhone.data.zro[slicer,0,0]*units.density
    def _rhoset(slicer,val):
        vhone.data.zro[slicer,0,0] = val/units.density
    _rhostring = "Mass density of the gas in g/cm^3"
    _rho = _Field(_rhostring,_rhoget,_rhoset)
    rho = property(*propertyargs(_rho))
    """ 
    HYDROGEN NUMBER DENSITY
    This the density rho converted into a Hydrogen number density
    """
    def _nHget(slicer):
        #return self.rho[slicer]/units.mH*units.X
        return vhone.data.zro[slicer,0,0]/units.mH*units.X*units.density
    def _nHset(slicer,val):
        vhone.data.zro[slicer,0,0] = val*units.mH/units.X/units.density
    _nHstring = "Hydrogen number density of the gas in cm^{-3}"
    _nH = _Field(_nHstring,_nHget,_nHset)
    nH = property(*propertyargs(_nH))
    """ 
    TOTAL PRESSURE
    This is the thermal + magnetic pressure in each cell
    Note: in some codes this is the total energy in the cell
        (kinetic + thermal + magnetic + etc etc)
        VH1 does not include kinetic energy in this
    """
    def _Pget(slicer):
        return vhone.data.zpr[slicer,0,0]*units.pressure
    def _Pset(slicer,val):
        vhone.data.zpr[slicer,0,0] = val/units.pressure
    _Pstring = "Gas thermal + magnetic pressure in erg/cm^3 (note: does not include kinetic energy)"
    _P = _Field(_Pstring,_Pget,_Pset)
    P = property(*propertyargs(_P))
    """ 
    THERMAL PRESSURE
    Remove magnetic pressure from the gas and return thermal pressure
    Placeholder due to being semi-derived
    """
    _PThermal = _Field()
    PThermal = property(*propertyargs(_PThermal))
    """ 
    VELOCITY
    This is the velocity of the gas flow in each cell
    """
    def _velget(slicer):
        return vhone.data.zux[slicer,0,0]*units.velocity
    def _velset(slicer,val):
        vhone.data.zux[slicer,0,0] = val/units.velocity
    _velstring = "Gas velocity in cm/s (+ve away from centre)"
    _vel = _Field(_velstring,_velget,_velset)
    vel = property(*propertyargs(_vel))
    """
    CELL MASS
    Mass in each grid cell
    NOTE: Updating the mass modifies the velocity to preserve
            kinetic energy
            This is done because it is correct for e.g. wind injection
    """
    _mass = _Field()
    mass = property(*propertyargs(_mass))

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
    _T = _Field()
    T = property(*propertyargs(_T))
    """ 
    SOUND SPEED
    Related to temperature
    Changing it is similar to changing the temperature
    """
    _cs = _Field()
    cs = property(*propertyargs(_cs))
    """ 
    KINETIC ENERGY
    1/2 m v^2
    Changing it changes the velocity, not the density
    """
    _KE = _Field()
    KE = property(*propertyargs(_KE))
    """ 
    THERMAL ENERGY
    Uses ideal gas assuming atomic composition (3/2 P V)
    Changing it changes the pressure
    """
    _TE = _Field()
    TE = property(*propertyargs(_TE))
    """ 
    GRAVITATIONAL POTENTIAL ENERGY
    Calculates the GPE up to a cell radius from all the mass inside
    """
    _GPE = _Field()
    GPE = property(*propertyargs(_GPE))
    """ 
    MAGNETIC PRESSURE
    Magnetic pressure stored below
    """
    # Memory space for magnetic pressure
    _internalvariables["PMagnetic"] = None
    _PMagneticField = _Field()
    PMagnetic = property(*propertyargs(_PMagneticField))
    """ 
    MAGNETIC FIELD
    Magnetic field
    """
    _Bfield = _Field()
    Bfield = property(*propertyargs(_Bfield))
    """ 
    PYTHON-ONLY VARIABLES---------------------------------------
    These variables do not affect the underlying simulation code
    They are used in the Python-only modules (e.g. radiation)
    """
    """ 
    HYDROGEN IONISATION FRACTION xHII
    What fraction of hydrogen is ionised?
    """
    #_xhii = np.zeros(self.ncells)
    _xhiiField = PythonField("xii","Hydrogen ionisation fraction (from 0 to 1)")
    xhii = property(*propertyargs(_xhiiField))
    """ 
    METALLICITY Z
    What is the mass fraction of metals?
    Given as a fraction of solar (Z=0.014)
    """
    _ZsolarField = PythonField("Zsolar","Gas metallicity (in solar units, where Zsolar=0.014)")
    Zsolar = property(*propertyargs(_ZsolarField))
    """ 
    IONISING PHOTON RATE
    Rate of ionising photons through each cell
    """
    _QionField = PythonField("Qion","Ionising photon rate passing through each cell in photons/second")
    Qion = property(*propertyargs(_QionField))
    """ 
    DUST ABSORPTION CROSS SECTION
    Absorption cross-section of dust in each radial unit
    """
    _sigmaDustField = PythonField("sigmaDust","Dust cross section in cm^-2")
    sigmaDust = property(*propertyargs(_sigmaDustField))
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
    _gravstring = "Gravitational acceleration in cm/s^2"
    _grav = _Field(_gravstring,_Gget,_Gset)
    grav = property(*propertyargs(_grav))

    @property
    def gamma(self):   
        """ 
        Adiabatic index of the gas
        Can't be set directly by user after setup to a gettable-only property
        To change, run integrator.Reset and integrator.Setup
        
        Returns
        -------

        gamma: float
            adiabatic index
        """
        return vhone.data.gam

