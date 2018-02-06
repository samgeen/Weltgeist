'''
Engine for running the integrator
Sam Geen, February 2018
'''

from hydro import hydro
import sources, init, vhone,units

class Integrator(object):
    def __init__(self):
        init.init() # Make sure VH1 is initialised
        self._time = vhone.data.time+0.0
        self._dt = 0.0

    def Step(self):
        # Get dt from last time
        oldtime = self._time
        self._time = vhone.data.time+0.0
        self._dt = self._time - oldtime
        # Inject sources
        sources.InjectSources(self._time, self._dt)
        vhone.data.step()

    @property
    def time(self):
        return self._time*units.time
    
    @property
    def dt(self):
        return self._dt*units.time

integrator = Integrator()