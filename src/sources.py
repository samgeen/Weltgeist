"""
Manages sources of energy/mass 
This can be constant sources, or using stellar evolution tracks
Sam Geen, January 2018
"""

import abc

import numpy as np

import singlestar, units, radiation, integrator

# This should be a singleton object really
# OH WELL
_sources = []
_totalke = 0.0
_totalmass = 0.0
_totalphotons = 0.0


def Sources():
    """
    Returns the sources object; should only be one of them
    """
    return _sources

class _Injector(object):
    """
    This object controls injection onto the grid
    """
    def __init__(self):
        """
        Constructor
        """
        self._totalke = 0.0
        self._totalmass = 0.0
        self._totalphotons = 0.0

    def InjectSources(self, sources):       
        """
        Runs through all of the sources and injects them onto the grid
        """
        hydro = integrator.Integrator().hydro
        # Clear values first
        self._totalke = 0.0
        self._totalmass = 0.0
        self._totalphotons = 0.0
        # Add to the input arrays
        for source in sources:
            source.Inject(self)
        # Dump input values onto the grid
        if self._totalmass > 0:
            hydro.mass[0] += self._totalmass
        if self._totalke > 0:
            hydro.KE[0] += self._totalke
        if self._totalphotons > 0:
            radiation.trace_radiation(self._totalphotons)

    def AddKE(self, ke):
        """
        Add kinetic energy to grid

        Parameters
        ----------
        ke : float
            kinetic energy to add in erg
        """
        self._totalke += ke
    
    def AddMass(self, mass):
        """
        Add mass to grid

        Parameters
        ----------
        mass : float
            mass to add in g
        """
        self._totalmass += mass

    def AddPhotons(self, photons):
        """
        Add photons to grid

        Parameters
        ----------
        photons : float
            photons to add (number)
        """
        self._totalphotons += photons


class _Sources(object):
    """
    Controls the injection of feedback from different sources
    """
    def __init__(self):
        """
        Constructor
        """
        # List of sources
        self._sources = []
        # Injector object that handles injecting values onto grid
        self._injector = _Injector()

    def AddSource(self, source):
        """
        Add a source object to the list of sources

        Parameters
        ----------
        source : Source
            Source object to add to active sources
        """
        self._sources.append(source)

    def RemoveSource(self, source):
        """
        Remove a source object
        
        Parameters
        ----------
        source : Source
            Source object to remove from active sources
        """
        self._sources.remove(source)

    def AddTableSource(self, tablesource):
        """
        Inject a single star using the singlestar tables
        (TODO: ADD SNe)

        Parameters
        ----------
        tablesource : TableSource
            A table source object
        """
        self.AddSource(tablesource)

    def MakeSupernova(self, energy, mass, time=0.0):
        """
        Inject a supernova

        Parameters
        ----------

        energy : float
            Energy to inject in erg
        mass : float
            Mass to inject in g
        time : float
            Time of supernova in s
        """
        source = SupernovaSource(energy,mass,time)
        self.AddSource(source)

    def MakeWind(self, lum, massloss):
        """
        Inject a wind source
        
        Parameters
        ----------

        lum : float
            Luminosity (energy per unit time) to inject in erg / s
        massloss : float
            Mass to inject per unit time in g/s
        """
        source = WindSource(lum,massloss)
        self.AddSource(source)

    def MakeRadiation(self, Sphotons):
        """
        Inject a radiation source
        """
        source = RadiationSource(Sphotons)
        self.AddSource(source)

    def InjectSources(self):
        """
        Gathers all of the sources and injects them onto the grid
        """
        # This role is passed onto the injector
        self._injector.InjectSources(self._sources)

# Singleton sources object - access it via Sources()
_sources = _Sources()

# TODO: ABSTRACT THIS INTO WINDS, SN, RADIATION
class AbstractSource(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self):
        pass

    @abc.abstractmethod
    def Inject(self,injector):
        """
        Template function for injecting feedback
        """
        pass

class SupernovaSource(AbstractSource):
    def __init__(self,energy,mass,time=0.0):
        self._energy = energy
        self._mass = mass
        self._time = time
        self._exploded = False

    def Inject(self,injector):
        """
        Injects the supernova
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """
        # Should the SN happen?
        # TODO: Shorten dt so that the SN happens exactly on time
        if integrator.Integrator().time >= self._time and not self._exploded:
            self._exploded = True
            injector.AddMass(self._mass)
            injector.AddKE(self._energy)
            # Remove the SN so it doesn't happen again next timestep
            # (Unnecessary since self._exploded is set, but makes things a bit quicker)
            RemoveSource(self)

            
class WindSource(AbstractSource):
    def __init__(self,lum,massloss):
        # Luminosity in ergs/s
        self._lum = lum
        # Mass loss in g/s
        self._massloss = massloss
        # Velocity of wind, used to control the timestep
        self._vcourant = np.sqrt(2.0*lum/massloss)

    def Inject(self,injector):
        """
        Injects wind over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """
        # Inject constant source of wind mass & energy (as pure KE)
        dt = integrator.Integrator().dt
        injector.AddMass(self._massloss*dt)
        injector.AddKE(self._lum*dt)
        integrator.Integrator().CourantLimiter(self._vcourant)

class RadiationSource(AbstractSource):
    """
    A source of radiation
    """
    def __init__(self,QH):
        """
        Constructor

        QH : float
            ionising photon emission rate (in s^-1)
        """
        self._QH = QH

    def Inject(self,injector):
        """
        Injects photons over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """
        dt = integrator.Integrator().dt
        injector.AddPhotons(self._QH*dt)

class TableSource(AbstractSource):
    """
    Source of energy & photons based on a lookup table
    """
    def __init__(self,mass,tbirth=0.0,rad=True,wind=True):
        self._tbirth = tbirth
        self._mass = mass
        self._rad = rad
        self._wind = wind

    def Inject(self,injector):
        """
        Injects feedback from single stars over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """

        # Check that the table is set up
        # NOTE: Make sure singlestarLocation is set before you get here
        self._TableSetup()

        # Calculate the star's current age
        t = integrator.Integrator().time
        dt = integrator.Integrator().dt
        age = t-self._tbirth
        if age > 0.0:
            if self._wind:
                # Get energy and mass to inject
                energy, massloss = singlestar.star_winds(self._mass,age,dt)
                # Add mass FIRST since KE needs to be added elastically
                injector.AddMass(massloss)
                # Add energy to grid as kinetic energy
                injector.AddKE(energy)
            if self._rad:
                #DO RAD
                nphotons = singlestar.star_radiation(self._mass,age,dt)
                # NOTE !!! ASSUMES 5 RADIATION BINS!!!
                # 1 = IR, 2 = optical+FUV, 3 = H-ionising, 
                # 4 = HeI->HeII, 5 = HeII->HeIII and up
                nuv = np.sum(nphotons[3:5])
                injector.AddPhotons(nuv)

    def _TableSetup(self):
        """
        Set up the single star table if not done already
        Note that this relies on globals because f2py can't instance
         multiple versions of a Fortran library
        """
        global _tablesetup
        global singlestarLocation
        if not _tablesetup:
            singlestar.star_setup(self._tableloc)
            _tablesetup = True

# Location of single star tables
# NOTE: this needs to be set correctly before the single star module is used
singlestarLocation = "/home/samgeen/Programming/Astro/StellarSources/Compressed/singlestar_z0.014"
# Are the tables set up?
_tablesetup = False