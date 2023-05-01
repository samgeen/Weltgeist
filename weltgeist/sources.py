"""
Manages sources of energy/mass 
This can be constant sources, or using stellar evolution tracks
Sam Geen, January 2018
"""

import abc

import numpy as np

from . import singlestar, units, radiation, integrator

doRadiationPressure = True
courantWarningMade = False

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
        self._totalte = 0.0
        self._totalke = 0.0
        self._totalmass = 0.0
        self._totalLionising = 0.0
        self._totalLnonionising = 0.0
        self._Eionising = 0.0
        self._Tion = 0.0

    def InjectSources(self, sources):       
        """
        Runs through all of the sources and injects them onto the grid
        """
        global doRadiationPressure
        global alwaysCalculateRadiation
        hydro = integrator.Integrator().hydro
        # Clear values first
        self._totalte = 0.0
        self._totalke = 0.0
        self._totalmass = 0.0
        self._totalNionising = 0.0
        self._totalLionising = 0.0
        self._totalLnonionising = 0.0
        self._Eionising = 0.0
        self._Tion = 0.0

        # Add to the input arrays
        for source in sources.sources:
            source.Inject(self)

        # Dump input values onto the grid
        if self._totalmass > 0:
            hydro.mass[0] += self._totalmass
        if self._totalte > 0:
            hydro.TE[0] += self._totalte
        if self._totalke > 0:
            hydro.KE[0] += self._totalke

        # Turn radiation on if trying to inject sources
        if self._totalLionising > 0 or self._totalLnonionising > 0:
            radiation.radiation_on = True
        
        # If radiation turned on, inject photons
        if radiation.radiation_on:
            timer = integrator.Integrator().ProcessTimer()
            timer.Begin("radiation")
            radiation.trace_radiation(self._totalLionising, self._totalLnonionising, 
                                               self._Eionising, self._Tion, doRadiationPressure)
            timer.End("radiation")


    def AddTE(self, te):
        """
        Add thermal energy to grid

        Parameters
        ----------
        te : float
            thermal energy to add in erg
        """
        self._totalte += te

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

    def AddPhotons(self, Lionising, Lnonionising, Eionising, Tion, sigmaDust=None):
        """
        Add photons to grid

        Parameters
        ----------
        Lionising : float
            Ionising photon luminosity (erg/s)
        Lnonionising : float
            Non-ionising photon luminosity (erg/s)
        Eionising : float
            Energy of ionising photons (ergs)
        Tion : float
            Temperature of photoionised gas (K)
        sigmaDust : float
            If not None, set the sigmaDust of the gas
        """
        # Set photons to inject
        self._totalLionising += Lionising
        self._totalLnonionising += Lnonionising
        self._Eionising = max(self._Eionising,Eionising)
        self._Tion = max(self._Tion,Tion)


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

    def Reset(self):
        """
        Reset the sources object by removing all of the sources
        """
        self._sources = []

    @property
    def sources(self):
        """
        List of sources

        Returns
        -------
        sources : list
            list of sources contained in this class
        """
        return self._sources

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

    def MakeSimpleRadiation(self, QH, Tion=1e4):
        """
        Inject a radiation source

        Parameters
        ----------

        QH : float
            Emission rate of ionising photons in number / s

        Tion : float
            Temperature of the photoionised gas in K (Default: 1e4 K)
        """
        source = SimpleRadiationSource(QH,Tion=Tion)
        self.AddSource(source)

    def InjectSources(self):
        """
        Gathers all of the sources and injects them onto the grid
        """
        # This role is passed onto the injector
        self._injector.InjectSources(self)

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
        integrator.Integrator().ForceTimeTarget(self._time)
        if integrator.Integrator().time >= self._time and not self._exploded:
            self._exploded = True
            injector.AddMass(self._mass)
            injector.AddKE(self._energy)
            # Remove the SN so it doesn't happen again next timestep
            # (Unnecessary since self._exploded is set, but makes things a bit quicker)
            injector.RemoveSource(self)

            
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

class SimpleRadiationSource(AbstractSource):
    """
    A simple source of ionising radiation
    """
    def __init__(self,QH,Tion=1e4):
        """
        Constructor

        QH : float
            ionising photon emission rate (in s^-1)
        """
        self._QH = QH
         # Set average photon energy to 13.6 eV 
        self._ephoton = 13.6 * units.eV
        # Set average temperature of photoionised gas to 10000 K
        self._Tion = Tion

    def Inject(self,injector):
        """
        Injects photons over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """
        dt = integrator.Integrator().dt
        injector.AddPhotons(self._QH*self._ephoton,
                            0.0,
                            self._ephoton,
                            self._Tion,
                            0.0)

class TableSource(AbstractSource):
    """
    Source of energy & photons based on a lookup table
    """
    def __init__(self,mass,tbirth=0.0,radiation=True,wind=True,supernova=False):
        """
        Constructor
    
        Parameters
        ----------

        mass : float
            Mass of star in solar masses
        tbirth : float
            Birth time of the star in seconds
        radiation : bool
            Turn radiation on?
        wind : bool
            Turn winds on?
        """
        self._tbirth = tbirth
        self._mass = mass
        self._radiation = radiation
        self._wind = wind
        self._supernova = supernova
        self._exploded = False
        # Check that the table is set up
        # NOTE: Make sure singlestarLocation is set before you get here
        self._TableSetup()
        # Set time the supernova should go off in the simulation
        self._supernovaTime = self._tbirth + singlestar.star_lifetime(self._mass)

    def Inject(self,injector):
        """
        Injects feedback from single stars over the timestep
        
        Parameters
        ----------

        injector : _Injector object
            object that accepts values to inject to grid
        """
        global starmetal
        global courantWarningMade

        # Calculate the star's current age
        t = integrator.Integrator().time
        dt = integrator.Integrator().dt
        age = t-self._tbirth
        Teff = None
        # Check whether the star is "alive" or not
        if age > 0.0 and not self._exploded:
            # Do supernova
            # Check first whether the star should explode before putting in pre-supernova feedback
            if self._supernova:
                # TODO: Fix timestep to make SN at exact time of supernova
                integrator.Integrator().ForceTimeTarget(self._supernovaTime)
                if t >= self._supernovaTime:
                    self._exploded = True
                    snEnergy, snMassLoss, snYield = singlestar.star_supernovae(self._mass)
                    # Inject 80% of the star's initial mass and 1e51 ergs kinetic energy
                    print("TableSource: Injecting supernova with energy, mass, at time", snEnergy, snMassLoss, self._supernovaTime)
                    injector.AddMass(snMassLoss)
                    injector.AddKE(snEnergy)
            if not self._exploded:
                # Do stellar winds
                if self._wind:
                    # Get energy and mass to inject
                    energy, massloss = singlestar.star_winds(self._mass,age,dt)
                    # Add mass FIRST since KE needs to be added elastically
                    injector.AddMass(massloss)
                    # Add energy to grid as kinetic energy
                    injector.AddKE(energy)
                    # Add some thermal energy to account for star's temperature
                    Teff = singlestar.star_effectivetemperature(self._mass,age) # Kelvin
                    TE = 1.5 * units.kB * massloss/(units.mH/units.X)*Teff
                    injector.AddTE(TE)
                    # Set the Courant condition
                    ecourant, mcourant = singlestar.star_winds(self._mass,age,1.0)
                    # Check whether the courant limiter is trying to put in zero mass
                    try:
                        vwind = np.sqrt(2.0*ecourant/mcourant)
                    except:
                        if not courantWarningMade:
                            print("Incorrect value of either ecourant, mcourant, ignoring courant limiter:",ecourant,mcourant)
                            print("(This warning only made once to prevent warning spam)")
                            courantWarningMade = True
                    else:
                        integrator.Integrator().CourantLimiter(vwind)
                        
                # Do stellar radiation
                if self._radiation:
                    # Read single star tables:
                    # 1. Number of photons emitted per s
                    Qphotons = singlestar.star_radiation(self._mass,age,1.0)
                    # 2. Photon luminosities
                    photonbands = singlestar.star_bandenergies(self._mass,age,1.0)
                    # 3. Ionised gas temperature
                    if Teff is None:
                        Teff = singlestar.star_effectivetemperature(self._mass,age)
                    Tion = radiation.IonisedGasTemperature(Teff, starmetal)
                    # Get the ionising photon band
                    # Assumes Lbolometric (erg/s) in position 2 and 
                    #  Lionising (erg/s) in position 3
                    Lionising = photonbands[3]
                    # Get the non-ionising photon band
                    Lbol = photonbands[2]
                    Lnonionising = Lbol - Lionising
                    # Assumes:
                    # 1 = IR, 2 = optical+FUV, 3 = H-ionising, 
                    # 4 = HeI->HeII, 5 = HeII->HeIII and up
                    QH = np.sum(Qphotons[2:5])
                    # Get energy of an ionising photon
                    Eionising = Lionising / QH
                    injector.AddPhotons(Lionising, Lnonionising, Eionising, Tion)

    def _TableSetup(self):
        """
        Set up the single star table if not done already
        Note that this relies on globals because f2py can't instance
         multiple versions of a Fortran library
        """
        global _tablesetup
        global singlestarLocation
        if not _tablesetup:
            singlestar.star_setup(singlestarLocation)
            _tablesetup = True

# Location of single star tables
# NOTE: this needs to be set correctly before the single star module is used
starmetal = 0.014
singlestarLocation = "/home/samgeen/Programming/Astro/StellarSources/Compressed/singlestar_z"+str(starmetal)
# Are the tables set up?
_tablesetup = False