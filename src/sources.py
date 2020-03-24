'''
Manages sources of energy/mass 
Sam Geen, January 2018
'''

import abc

import numpy as np

import singlestar, units, radiation
from hydro import hydro

# This should be a singleton object really
# OH WELL
_sources = []
_totalke = 0.0
_totalmass = 0.0
_totalphotons = 0.0
_tablesetup = False
_tableloc = "/home/samgeen/Programming/Astro/StellarSources/Compressed/singlestar_z0.014"
# Temp variables to prevent re-definitions
ngroups = 5 # Hard-coded for benefit of Fortran module
_uv = np.zeros(ngroups,dtype="float64")

def TableSetup():
    global _tablesetup
    if not _tablesetup:
        singlestar.star_setup(_tableloc)
        _tablesetup = True

def AddSource(source):
    '''
    Add a source object to the list of sources
    '''
    global _sources
    _sources.append(source)

def RemoveSource(source):
    '''
    Remove a source object
    '''
    global _sources
    _sources.remove(source)
    print(len(_sources))

def MakeStar(mass,tbirth=0.0,rad=True,wind=True):
    '''
    Inject a single star of a given mass (TODO: ADD RADIATION AND SNe)
    '''
    source = TableSource(mass,tbirth,rad,wind)
    AddSource(source)

def MakeSupernova(energy,mass,time=0.0):
    '''
    Inject a supernova
    '''
    source = SupernovaSource(energy,mass,time)
    AddSource(source)

def MakeWind(lum,massloss,time=0.0):
    '''
    Inject a wind source
    '''
    source = WindSource(lum,massloss)
    AddSource(source)

def MakeRadiation(Sphotons):
    '''
    Inject a radiation source
    '''
    source = RadiationSource(Sphotons)
    AddSource(source)

def InjectSources(t,dt):
    global _totalke
    global _totalmass
    global _totalphotons
    _totalke = 0.0
    _totalmass = 0.0
    _totalphotons = 0.0
    # Add to the input arrays
    for source in _sources:
        source.Inject(t,dt)
    # Dump input values
    if _totalmass > 0:
        hydro.mass[0] += _totalmass
    if _totalke > 0:
        hydro.KE[0] += _totalke
    if _totalphotons > 0:
        radiation.trace_radiation(_totalphotons)

# TODO: ABSTRACT THIS INTO WINDS, SN, RADIATION
class AbstractSource(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self):
        pass

    @abc.abstractmethod
    def Inject(self,t,dt):
        pass

class SupernovaSource(AbstractSource):
    def __init__(self,energy,mass,time=0.0):
        self._energy = energy
        self._mass = mass
        self._time = time
        self._exploded = False

    def Inject(self,t,dt):
        global _sources
        global _totalmass
        global _totalke
        # Should the SN happen?
        # TODO: Shorten dt so that the SN happens exactly on time
        if t >= self._time and not self._exploded:
            self._exploded = True
            _totalmass += self._mass
            _totalke += self._energy
            #hydro.TE[0] += self._energy
            #print "KE", hydro.KE[0:10], units.energy
            #print "P", hydro.P[0:10], units.pressure
            #print "T", hydro.T[0:10]
            #print "v", hydro.vel[0:10], units.velocity
            #print "rho", hydro.rho[0:10], units.density
            # Remove the SN so it doesn't happen again next timestep
            # (Unnecessary since self._exploded is set, but makes things a bit quicker)
            RemoveSource(self)

            
class WindSource(AbstractSource):
    def __init__(self,lum,massloss):
        # Luminosity in ergs/s
        self._lum = lum
        # Mass loss in g/s
        self._massloss = massloss
        self._vcourant = np.sqrt(2.0*lum/massloss)

    def Inject(self,t,dt):
        global _sources
        global _totalmass
        global _totalke
        # Inject constant source of wind mass & energy (as pure KE)
        # Convert dt to seconds from internal units
        _totalmass += self._massloss*dt
        _totalke += self._lum*dt
        hydro.CourantLimiter(self._vcourant)

class RadiationSource(AbstractSource):
    def __init__(self,Sphotons):
        self._Sphotons = Sphotons

    def Inject(self,t,dt):
        global _totalphotons
        _totalphotons += self._Sphotons


class TableSource(AbstractSource):
    '''
    Source of energy & photons based on a lookup table
    '''
    global _ewind, _mlwind
    def __init__(self,mass,tbirth=0.0,rad=True,wind=True):
        self._tbirth = tbirth
        self._mass = mass
        self._rad = rad
        self._wind = wind

    def Inject(self,t,dt):
        global _totalmass
        global _totalke
        global _totalphotons
        uv = np.a
        age = t-self._tbirth
        if age > 0.0:
            if self._wind:
                # Get energy and mass to inject
                energy, massloss = singlestar.star_winds(self._mass,age,dt)
                # Add mass FIRST since KE needs to be added elastically
                _totalmass += massloss
                # Add energy to grid as kinetic energy
                _totalke += energy
            if self._rad:
                #DO RAD
                nphotons = singlestar.star_radiation(self._mass,age,dt)
                # NOTE !!! ASSUMES 5 RADIATION BINS!!!
                nuv = np.sum(nphotons[3:5])
                _totalphotons += nuv/dt # Solver needs S* (photons / second)