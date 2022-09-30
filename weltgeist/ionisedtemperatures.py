"""
Calculate the temperature of photoionised gas
Sam Geen, April 2019
"""

# Apologies, this is garbage code because I wrote it quickly for another project

import numpy as np
import scipy.interpolate

import weltgeist

import os
from pathlib import Path

# Empty objects to fill in _readtempertaures
_temperatureFunc = None
_Teffs = None
_nHs = None
_Zs = None

# Get the current source location to read 
source_path = Path(__file__).resolve()
source_dir = source_path.parent

# Get the default gas temperature table location (will be copied in setup.py with the code)
defaulttableloc = str(source_dir)+os.sep+"Tgas.csv"

# What Zsolar does Cloudy assume? TODO: CHECK!!!
Zsolar = 0.014

def _readtemperatures(filename=defaulttableloc):
    global _temperatureFunc, _Teffs, _nHs, _Zs
    # NOTE: we're assuming an ionisation parameter U=-2
    if _temperatureFunc is None:
        # First rip the raw tables
        raw = np.loadtxt(filename,delimiter=",")
        #headers = "Teff,Styp,nH,U,dZ,T_H0,T_H+,T_H2".split(",")
        rawTeff = raw[:,0]
        rawnH = raw[:,2]
        rawU = raw[:,3]
        rawZ = raw[:,4]
        rawTi = raw[:,6]
        # Second, make a grid of values
        _Teffs = np.unique(rawTeff)
        _nHs = np.unique(rawnH)
        _Zs = np.unique(rawZ)
        Tis = np.zeros([len(_Teffs),len(_nHs),len(_Zs)])
        for i,Teff in enumerate(_Teffs):
            for j,nH in enumerate(_nHs):
                for k,Z in enumerate(_Zs):
                    Tis[i,j,k] = 10.0**rawTi[(rawTeff == Teff)*(rawnH == nH)*(rawZ == Z)*(rawU == -2)]
        # Make the interpolation function
        # NOTE: here we are explicitly using the strategy of extrapolating values to account for lower temperatures
        #       THIS IS A HUGE RISK - get Eric to make new tables???
        _temperatureFunc = scipy.interpolate.RegularGridInterpolator((_Teffs,_nHs,_Zs),Tis,bounds_error=False,fill_value=None)

def FindTemperature(Teffin,metalin):
    global _temperatureFunc
    # Make sure the tables are read
    dir = weltgeist.__file__
    dir = dir[:len(dir)-len("__init__.py")]
    _readtemperatures(dir+"Tgas.csv")
    # HACK - Hard code nH to 10^1.75 for simplicity
    lognH = 1.75
    # Clamp values to limits of arrays
    #print 10.0**_logtemperatureFunc(np.array([_Teffs.min(),1.75,0.0]))
    #print Teffin, np.log10(metalin/Zsolar), _Teffs, _Zs
    Teff = Teffin
    metal = metalin
    # Return the interpolation (nH, Z are in log units)
    Tout = _temperatureFunc(np.array([Teff,lognH,np.log10(metal/Zsolar)]))[0]
    return Tout

if __name__=="__main__":
    _readtemperatures()
    Ts = []
    import matplotlib.pyplot as plt
    Teffs = np.linspace(10000,50000,1000)
    for T in Teffs:
        Ts.append(FindTemperature(T,0.014))
    plt.plot(Teffs,Ts)
    plt.savefig("testionisedtemperatures.pdf")
