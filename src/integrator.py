'''
Engine for running the integrator
Sam Geen, February 2018
'''

from hydro import hydro
import cooling, gravity, init, sources, units, vhone

gravity_on = True
cooling_on = False

class Integrator(object):
    def __init__(self):
        init.init() # Make sure VH1 is initialised
        self._time_code = vhone.data.time+0.0
        self._dt_code = 0.0

    def Step(self):
        # Update time
        self._UpdateTime()
        # Cooling step
        if cooling_on:
            cooling.solve_cooling(self.dt)
        # Inject sources (includes radiation step!)
        sources.InjectSources(self.time, self.dt)
        # Gravity step
        if gravity_on:
            gravity.calculate_gravity()
        # Hydro step
        vhone.data.step()
        # Final sanity check
        if cooling_on:
            cooling.CheckTemperature()

    def _UpdateTime(self):
        # Get dt from last time
        oldtime = self._time_code
        self._time_code = vhone.data.time+0.0
        self._dt_code = self._time_code - oldtime

    @property
    def time(self):
        return self._time_code*units.time
    
    @property
    def dt(self):
        return self._dt_code*units.time

integrator = Integrator()