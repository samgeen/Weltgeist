"""
Tracks flows leaving the grid
Sam Geen, May 2023
"""

import numpy as np

class OutflowTracker():
    '''
    Tracks mass, momentum, enery and photons leaving the grid
    '''
    def __init__(self, hydro,cellindex=-1):
        self._hydro = hydro
        self._cellindex = cellindex
        self.Reset()

    def Reset(self):
        self._mass = 0.0
        self._momentum = 0.0
        self._energy = 0.0
        self._photons = 0.0

    @property
    def mass(self):
        '''
        Return the mass lost so far
        '''
        return self._mass
    
    @property
    def momentum(self):
        '''
        Return the momentum lost so far
        '''
        return self._momentum
    
    @property
    def energy(self):
        '''
        Return the energy lost so far
        '''
        return self._energy
    
    @property
    def photons(self):
        '''
        Return the number of photons lost so far
        '''
        return self._photons

    def TrackForStep(self, dt):
        '''
        Track for one step
        dt (float) : timestep in seconds
        '''
        h = self._hydro
        # Get the flow values in the last cell
        last = self._cellindex
        rho = h.rho[last]
        vel = h.vel[last]
        # If it's the last cell, don't include any likely spurious inflows
        if last == -1:
            vel = max(0.0,vel)
        vol = h.vol[last]
        # Get area of outer edge of cell
        area = 4.0 * np.pi * (h.x[last] + 0.5*h.dx)**2
        # Volume per timestep leaving the grid
        flowVolume = vel * area * dt
        # Calculate flows through the surface of the last cell
        self._mass += rho * flowVolume
        self._momentum += rho*vel * flowVolume
        self._energy += (h.KE[last]+h.TE[last])/vol * flowVolume
        # For photons, we already have a flow rate in Qion, 
        #   so just roughly integrate over the timstep
        self._photons += h.Qion[last] * dt

    def Save(self, fileh5py, version):
        '''
        Save results to file
        fileh5py: HDF5 file to save to
        version: version number in case format changes
        '''
        if version >= 1.04:
            outflows = np.array([self._mass,self._momentum,self._energy,self._photons])
            fileh5py.create_dataset("outflows",data=outflows,dtype=np.float64)

    def Load(self, fileh5py, version):
        '''
        Load results from file
        fileh5py: HDF5 file to load from
        version: version number in case format changes
        '''
        def loaditem(varname):
            data = np.array(fileh5py.get(varname))
            if len(data) == 1:
                data = data[0]
            return data
        if version >= 1.04:
            self._mass,self._momentum,self._energy,self._photons = loaditem("outflows")

    def __str__(self):
        '''
        Return a string with the reported outflow totals
        ''' 
        return "Outflows: (mass: "+str(self._mass)+" g, momentum: "+str(self._momentum)+" g cm/s, " + \
                        "energy: "+str(self._energy)+" erg, photons: "+str(self._photons)+")"
