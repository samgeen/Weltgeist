'''
Created on 3 Jun 2014

@author: samgeen
'''

import numpy as np

import pyglet
from pyglet.gl import *

from hydro import hydro
from integrator import integrator
import sources, gravity, init, units, windsolutions

class Line(object):
    def __init__(self,colour,width=1.0):
        self._colour = colour
        self._currlen = 0
        self._x = np.zeros(100)
        self._y = np.zeros(100)
        self._pts = np.zeros((100,2))
        self._ymin = None
        self._ymax = None
        self._width = width

    def Update(self,x,y):
        '''
        Replace arrays
        '''
        newlen = len(x)
        self._currlen = newlen
        if newlen > self._x.shape[0]:
            self._x = np.zeros(2*newlen)
            self._y = np.zeros(2*newlen)
            self._pts = np.zeros((2*newlen,2))
        self._x[0:newlen] = x
        self._y[0:newlen] = y

    def Append(self,x,y):
        '''
        Add a single point
        '''
        if self._currlen >= len(self._x):
            xnew = np.zeros(self._currlen*2)
            ynew = np.zeros(self._currlen*2)
            self._pts = np.zeros((2*self._currlen,2))
            xnew[0:self._currlen] = self._x
            ynew[0:self._currlen] = self._y
            self._x = xnew
            self._y = ynew
        self._x[self._currlen] = x
        self._y[self._currlen] = y
        self._currlen += 1

    def MinMax(self,yminin=None,ymaxin=None):
        '''
        Return min/max in arrays
        '''
        if yminin is not None:
            self._ymin = yminin
        if ymaxin is not None:
            self._ymax = ymaxin
        x = self._x[0:self._currlen]
        y = self._y[0:self._currlen]
        if len(x) == 0:
            return 0.0,0.0,0.0,0.0
        else:
            if self._ymin is not None:
                ymin = self._ymin
            else:
                ymin = y.min()
            if self._ymax is not None:
                ymax = self._ymax
            else:
                ymax = y.max()
            return x.min(), x.max(), ymin, ymax

    def Draw(self,xmin,xmax,ymin,ymax):
        '''
        Draw line to screen
        '''
        if self._currlen > 0:
            glLineWidth(self._width)
            self._pts[:,0] = (self._x - xmin) / (xmax-xmin)
            self._pts[:,1] = (self._y - ymin) / (ymax-ymin)
            pts = self._pts[0:self._currlen,0:2].flatten()
            ptslen = len(pts)
            glColor4f(*self._colour)
            glPoints = (GLdouble * ptslen)(*pts)
            glVertexPointer(2, GL_DOUBLE, 0, glPoints)
            glDrawArrays(GL_LINE_STRIP,0,len(pts)//2)


class Lines(object):
    def __init__(self,lines,width=1.0):
        self._lines = lines
        self._xmin = 1e30
        self._xmax = -1e30
        self._ymin = 1e30
        self._ymax = -1e30

    def _FindMinMax(self):
        first = True
        self._xmin = 1e30
        self._xmax = -1e30
        self._ymin = 1e30
        self._ymax = -1e30
        for line in self._lines:
            if first:
                self._xmin, self._xmax, self._ymin, self._ymax = line.MinMax()
                first = False
            else:
                xl,xh,yl,yh = line.MinMax()
                self._xmin = min(self._xmin,xl)
                self._xmax = max(self._xmax,xh)
                self._ymin = min(self._ymin,yl)
                self._ymax = max(self._ymax,yh)

    def Draw(self):
        self._FindMinMax()
        for line in self._lines:
            line.Draw(self._xmin,self._xmax,self._ymin,self._ymax)

class Tester(object):
    def __init__(self):
        self._glPoints = None
        self._window = pyglet.window.Window(512,512)
        self._window._integrator = self
        self._mx = 0.5
        self._my = 0.5
        self._step = True
        self._itick = 0
        self._rvtline = Line([0.0,0.0,1.0,1.0],width=3.0)
        self._sedov = Line([0.0,1.0,0.0,1.0])
        self._lines = Lines([self._sedov,self._rvtline])
        self._first = True
        
    def Setup(self):

        @self._window.event
        def on_draw():
            glClearColor(1,1,1,1)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glDisable(GL_DEPTH_TEST)
            glDisable(GL_TEXTURE_2D)
            glDisable(GL_LIGHTING)
            glViewport(0, 0, self._window.width, self._window.height)
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            glOrtho(-1e-2, 1.0+1e-2,
                    -1e-2, 1.0+1e-2,
                    -2, 2)
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            glEnableClientState(GL_VERTEX_ARRAY)
            self._lines.Draw()
            
        @self._window.event
        def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
            w, h = self._window.get_size()
            x = x/float(w)
            y = y/float(h)
            mx, my = self._mx, self._my
            if x >= 1.0:
                return
            if x <= 0.0:
                return
            if y >= 1.0:
                return
            if y <= 0.0:
                return
            self._mx = x
            self._my = y
            #self._step = buttons == 1
            
        @self._window.event
        def on_mouse_release(x, y, button, modifiers):
            pass
            #self._step = False
            
        # Initialise sim data
        init.init()
        #sources.MakeSupernova(1e51,2e33)
        self.windlum = 3e34 # 2e38
        self.windml = 1e19 # 2e22
        #sources.MakeWind(self.windlum,self.windml)
        self.Sphotons = 1e48
        #sources.MakeRadiation(self.Sphotons)
        # Set up rendering
        s = 0.02
        glEnable(GL_BLEND)
        glPointSize(2)
        # Set up pyglet to run
        pyglet.clock.set_fps_limit(60)
        pyglet.clock.schedule_interval(self.Step,1.0/60.0)
        pyglet.app.run()
        
    def Step(self, dt):
        global sn, testfile
        nx = hydro.ncells
        if self._step:
            if self._first:
                #testfile = open("testdata.txt","w")
                # Do something here if needed
                #hydro.T[0:nx] = 2.0/3.0 * units.kB / units.mp * units.G * np.cumsum(hydro.mass[0:nx]) / hydro.x[0:nx]
                #hydro.T[0] = hydro.T[1]
                gravity.centralmass = 2e38
                hydro.nH[:512] = 1e-7
                #hydro.nH[513:] = 1e-6
                #hydro.mass[0] = 2e34
                rho = hydro.rho[0:nx]  
                self.x0 = 0.0 # hydro.x[np.where(hydro.nH[0:nx] > 1e-2)][0]
                hydro.T[0:nx] = 10
                print "NHMAXFIRST", hydro.nH.max()
                self._first = False
            self._itick += 1
            #hydro.nH[:30] = 1e-6
            hydro.CourantLimiter(1e6)
            integrator.Step()
            rs = hydro.nH[0:nx].max()# - init.n0
            t = integrator.time
            rho = hydro.rho[0:nx]
            nH = hydro.nH[0:nx]
            vel = hydro.vel[0:nx]
            T = hydro.T[0:nx]
            Tmax = T.max()
            Tmin = T.min()
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
            rs = rT
            rs = np.max(hydro.nH[0:nx])# - init.n0
            rpeak = hydro.x[np.where(nH > 100)][0]
            tstart = units.time*0.1
            if self.x0 == 0.0 or t < tstart:
                self.x0 = rpeak
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
            #rsedov = windsolutions.AdiabaticWind(self.windlum,init.n0,integrator.time)#,model="Avedisova")
            # Radiation
            #rsedov = windsolutions.SpitzerSolution(self.Sphotons,init.n0,t)
            # Gravity
            #tgrav = windsolutions.CollapseSolution(rs*units.mp,init.rho0)
            tgrav = windsolutions.CollapseSolutionPosition(rpeak,self.x0+0.0)
            #self._rvtline.MinMax(0.0,max(rsedov,rs))
            self._rvtline.Update(hydro.x[0:nx],np.log10(hydro.rho[0:nx]))
            #self._rvtline.Update(hydro.x[0:nx],hydro.T[0:nx])
            #tff = np.sqrt(3.0*np.pi/(32.0*units.G*init.rho0))
            #rsedov = init.n0/(1.0-t/tff)**1.5# - init.n0
            #print rpeak
            if (t > tstart):
                #testfile = open("testdata.txt","w")
                #testfile.write(str(t/tff)+" "+str(rs/init.n0)+"\n")
                #testfile.flush()
                #self._rvtline.Append(t,rpeak)
                #self._sedov.Append(t,rsedov)
                #self._sedov.Append(tgrav,rpeak)
                if self._itick == 1:
                    #print "t, Rsim, Rsedov, ratio, Tmin/max", t/units.time, rs, rsedov, rs / rsedov, Tmin, Tmax
                    print "t, tgrav, t/grav, rpeak, r0", t/units.time, tgrav/units.time, t/tgrav, rpeak/units.pc, self.x0 / units.pc
        if self._itick > 10:
            #import pdb; pdb.set_trace()
            self._itick = 0
        


if __name__=="__main__":
    test = Tester()
    test.Setup()
    
