'''
Created on 3 Jun 2014

@author: samgeen
'''

import numpy as np

from hydro import hydro
from integrator import integrator
import sources, gravity, init, units, windsolutions, graphics


class Tester(object):
    def __init__(self):
        self._step = True
        self._first = True
        self._itick = 0
        self._rvtline = graphics.Line([0.0,0.0,1.0,1.0],width=3.0)
        self._anline1 = graphics.Line([0.0,1.0,0.0,1.0])
        self._anline2 = graphics.Line([0.0,1.0,0.0,1.0])
        self._lines = [self._anline1,self._anline2,self._rvtline]
        self._renderer = graphics.Renderer(self._lines)
        
    def Setup(self):
            
        # Initialise sim data
        init.init()
        #sources.MakeSupernova(1e51,2e33)
        self.windlum = 3e34 # 2e38
        self.windml = 1e19 # 2e22
        #sources.MakeWind(self.windlum,self.windml)
        self.Sphotons = 1e49
        sources.MakeRadiation(self.Sphotons)
        # Display and run
        self._renderer.Start(self.Step)
        
    def Step(self, dt):
        global sn, testfile
        nx = hydro.ncells
        # Some initialisation
        n0 = 1000.0
        rho0 = n0*units.mp
        rhoout = rho0
        r0 = 0.1*init.rmax
        T0 = 10.0
        if self._step:
            if self._first:
                hydro.nH[0:nx] = n0 * (hydro.x[0:nx]/r0)**(-2.0)
                hydro.nH[0] = hydro.nH[1]
                hydro.P[0:nx] = hydro.nH[0:nx] * (units.kB * T0)
                #sources.MakeRadiation(1e49)
                #testfile = open("testdata.txt","w")
                # Do something here if needed
                #hydro.T[0:nx] = 2.0/3.0 * units.kB / units.mp * units.G * np.cumsum(hydro.mass[0:nx]) / hydro.x[0:nx]
                #hydro.T[0] = hydro.T[1]
                #gravity.centralmass = 2e38
                #hydro.nH[:512] = 1e-7
                #hydro.nH[513:] = 1e-6
                #hydro.mass[0] = 2e34
                #cs = 1e5
                #rhoc = 1000.0*units.mp
                #rc = cs / np.sqrt(4.0*np.pi*units.G*rhoc)
                #print rc
                #for i in range(0,100):
                #    r = hydro.x[0:nx]
                #    hydro.cs[0:nx] = cs
                #    hydro.rho[0:nx] = rhoc * (1.0+(r/2.88/rc)**2.0)**(-1.47)
                #self.rhoBE = hydro.rho[0:nx]
                #self._sedov.Update(hydro.x[0:nx],np.log10(self.rhoisoT))
                #T = hydro.T[0:nx]
                #Tmax = T.max()
                #Tmin = T.min()
                #print Tmax, Tmin
                #rho = hydro.rho[0:nx]  
                #self.x0 = 0.0 # hydro.x[np.where(hydro.nH[0:nx] > 1e-2)][0]
                #print "NHMAXFIRST", hydro.nH.max()
                self._first = False
            self._itick += 1
            #hydro.nH[:30] = 1e-6
            #hydro.CourantLimiter(1e5)
            rho = hydro.rho[0:nx]
            rhomax = rho.max()
            nH = hydro.nH[0:nx]
            vel = hydro.vel[0:nx]
            T = hydro.T[0:nx]
            Tmax = T.max()
            Tmin = T.min()
            #print Tmax, Tmin
            integrator.Step()
            #hydro.T[0:nx] = 10.0
            #hydro.rho[0:5] = self.rhoisoT[0:5]
            rs = hydro.nH[0:nx].max()# - init.n0
            t = integrator.time
            #Tback = T[-1]
            #print Tback
            try:
                rT = hydro.x[np.where(T > 1e3)[0][-1]]
            except:
                rT = 0.0
            #rhoback = rho[-1]
            #rs = hydro.x[np.where(vel > 1e1)[0][-1]]
            try:
                rrho = hydro.x[np.where(rho > 1.0001*rho[-1])[0][-1]]
            except:
                rrho = 0.0
            rs = rrho
            rs = hydro.x[np.where(rho == rho.max())[0][0]]
            rhoout = rho0 * (rs/r0)**(-2.0)
            #rs = np.max(hydro.nH[0:nx])# - init.n0
            #GPE = hydro.GPE[0:nx]
            GPE = 0.0
            TE = hydro.TE[0:nx]
            KE = hydro.KE[0:nx]
            allenergy = np.sum(GPE+TE+KE)
            allge = np.sum(GPE)
            allkte = np.sum(TE+KE)
            #rpeak = hydro.x[np.where(nH > 100)][0]
            #tstart = units.time*0.1
            #if self.x0 == 0.0 or t < tstart:
            #    self.x0 = rpeak
            #rs = hydro.TE[10] / hydro.KE[10]
            #hydro.vel[0:nx] -= np.sqrt(2.0*0.9*hydro.TE[0:nx]/hydro.mass[0:nx])
            #hydro.TE[0:nx] *= 0.1
            #hydro.T[0:nx] = 10.0
            # beta = 0.968 for gamma = 1.4
            beta = 0.968
            beta = 1.1
            # Sedov blast
            #rsedov = beta*(1e51*(t**2.0) / init.rho0)**0.2
            # Winds
            #rcooled = windsolutions.AdiabaticWind(self.windlum,self.windml,init.n0,integrator.time,model="Cooled")
            rcooled = windsolutions.AdiabaticWind(self.windlum,self.windml,init.n0,integrator.time,model="SlowCool")
            rsedov = windsolutions.AdiabaticWind(self.windlum,self.windml,init.n0,integrator.time,model="WeaverIntermediate")
            #rsedov = windsolutions.AdiabaticWind(self.windlum,self.windml,init.n0,integrator.time,model="CooledNoAcc")
            # Radiation
            #rsedov = windsolutions.SpitzerSolution(self.Sphotons,init.n0,t)
            # Gravity
            #tgrav = windsolutions.CollapseSolution(rs*units.mp,init.rho0)
            #tgrav = windsolutions.CollapseSolutionPosition(rpeak,self.x0+0.0)
            #self._rvtline.MinMax(0.0,max(rsedov,rs))
            # LOG LOG SPACE
            #self._rvtline.Update(np.log10(hydro.x[1:nx]/hydro.x[1]),np.log10(hydro.rho[1:nx]/hydro.rho[1]))
            # LINEAR SPACE
            self._rvtline.Update(hydro.x[0:nx],np.log10(hydro.rho[0:nx]))
            #self._rvtline.Update(hydro.x[0:nx],np.log10(hydro.rho[0:nx]/self.rhoBE))
            #self._sedov.Update(hydro.x[0:nx],hydro.x[0:nx]*0.0)
            #self._rvtline.Update(hydro.x[0:nx],hydro.T[0:nx])
            #tff = np.sqrt(3.0*np.pi/(32.0*units.G*init.rho0))
            #rsedov = init.n0/(1.0-t/tff)**1.5# - init.n0
            #print rpeak
    
            if (rs > 0.0):
            #if (t > tstart):
                #testfile = open("testdata.txt","w")
                #testfile.write(str(t/tff)+" "+str(rs/init.n0)+"\n")
                #testfile.flush()
                #self._rvtline.Append(t,rs)
                #self._anline1.Append(t,rsedov)
                #self._anline2.Append(t,rcooled)
                if self._itick == 1:
                    print "t, Rsim, rhomax/rho0", t/units.Myr, rs/units.pc, rhomax / rhoout
                    #print "t, Rsim, Rsedov, ratio, Tmin/max, evir", t/units.Myr, rs/units.pc, rsedov/units.pc, \
                    #                     rs / rsedov, allkte-0.5*allge
                    #print "t, tgrav, t/grav, n, n0", t/units.time, tgrav/units.time, t/tgrav, rs, init.n0
        if self._itick > 10:
            #import pdb; pdb.set_trace()
            self._itick = 0
        


if __name__=="__main__":
    test = Tester()
    test.Setup()
    
