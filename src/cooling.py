'''
Controls gas cooling
Sam Geen, March 2018
'''

from hydro import hydro
import cooling_module, init

def solve_cooling(dt):
    '''
    Solve the cooling step
    dt - timestep in seconds
    '''
    ncell = hydro.ncells
    nH = hydro.nH[0:ncell]
    T2 = hydro.T[0:ncell]
    xhii = hydro.xhii
    zsolar = hydro.Zsolar[0:ncell]
    gamma = init.gamma
    dT2 = cooling_module.solve_cooling_frig(nH,T2,zsolar,dt,gamma,ncell)
    T2 += dT2
    hydro.T[xhii < 0.1] = T2[xhii < 0.1]
    CheckTemperature()

def CheckTemperature():
    '''
    Check for and fix super low tempertaures
    '''
    ncell = hydro.ncells
    T2 = hydro.T[0:ncell]
    # Extra check for low temperatures
    T2[T2 < 1.0] = 1.0
    hydro.T[0:ncell] = T2
