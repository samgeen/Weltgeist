'''
Created on 3 Jun 2014

@author: samgeen
'''

import numpy as np

import pyglet
from pyglet.gl import *

from hydro import hydro
from integrator import integrator
import sources, init, units, windsolutions

sn = False

class Line(object):
    def __init__(self,colour):
        self._colour = colour
        self._currlen = 0
        self._x = np.zeros(10000)
        self._y = np.zeros(10000)
        self._pts = np.zeros((10000,2))

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
            self._x = np.zeros(self._currlen*2)
            self._y = np.zeros(self._currlen*2)
            self._pts = np.zeros(2*self._currlen,2)
        self._x[self._currlen] = x
        self._y[self._currlen] = y
        self._currlen += 1

    def MinMax(self):
        '''
        Return min/max in arrays
        '''
        x = self._x[0:self._currlen]
        y = self._y[0:self._currlen]
        if len(x) == 0:
            return 0.0,0.0,0.0,0.0
        else:
            return x.min(), x.max(), y.min(), y.max()

    def Draw(self,xmin,xmax,ymin,ymax):
        '''
        Draw line to screen
        '''
        if self._currlen > 0:
            self._pts[:,0] = (self._x - xmin) / (xmax-xmin)
            self._pts[:,1] = (self._y - ymin) / (ymax-ymin)
            pts = self._pts[0:self._currlen,0:2].flatten()
            ptslen = len(pts)
            glColor4f(*self._colour)
            glPoints = (GLdouble * ptslen)(*pts)
            glVertexPointer(2, GL_DOUBLE, 0, glPoints)
            glDrawArrays(GL_LINE_STRIP,0,len(pts)//2)


class Lines(object):
    def __init__(self,lines):
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
        self._rvtline = Line([0.0,0.0,1.0,1.0])
        self._sedov = Line([0.0,1.0,0.0,1.0])
        self._lines = Lines([self._sedov,self._rvtline])
        
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
        sources.MakeSupernova(1e51,2e33)
        self.windlum = 2e38
        self.windml = 2e22
        #sources.MakeWind(self.windlum,self.windml)
        # Set up rendering
        s = 0.02
        glEnable(GL_BLEND)
        glPointSize(2)
        # Set up pyglet to run
        pyglet.clock.set_fps_limit(60)
        pyglet.clock.schedule_interval(self.Step,1.0/60.0)
        pyglet.app.run()
        
    def Step(self, dt):
        global sn
        nx = hydro.ncells
        if self._step:
            self._itick += 1
            integrator.Step()
            t = integrator.time
            rho = hydro.rho[0:nx]
            vel = hydro.vel[0:nx]
            T = hydro.T[0:nx]
            #Tback = T[-1]
            #print Tback
            #try:
            #    rs = hydro.x[np.where(T > 1e5)[0][-1]]
            #except:
            #    rs = 0.0
            #rhoback = rho[-1]
            #rs = hydro.x[np.where(vel > 1e1)[0][-1]]
            try:
                rs = hydro.x[np.where(rho > 1.0001*rho[-1])[0][-1]]
            except:
                rs = 1e-7
            # beta = 0.968 for gamma = 1.4
            beta = 0.968
            beta = 1.1
            # Sedov blast
            rsedov = beta*(1e51*(t**2.0) / init.rho0)**0.2
            # Winds
            #rsedov = windsolutions.AdiabaticWind(self.windlum,init.n0,integrator.time)#,model="Avedisova")
            if (rs > 0):
                #self._rvtline.Update(hydro.x[0:nx],hydro.rho[0:nx])
                self._rvtline.Append(t,rs)
                self._sedov.Append(t,rsedov)
                if self._itick == 1:
                    print "t, Rsim, Rsedov, ratio", t, rs, rsedov, rs / rsedov
        if self._itick > 10:
            #import pdb; pdb.set_trace()
            self._itick = 0
        


if __name__=="__main__":
    test = Tester()
    test.Setup()
    
