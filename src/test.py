'''
Created on 3 Jun 2014

@author: samgeen
'''

import init
import numpy as np

import pyglet
from pyglet.gl import *

from hydro import hydro

sn = False

class Integrator(object):
    def __init__(self):
        self._glPoints = None
        self._window = pyglet.window.Window(512,512)
        self._window._integrator = self
        self._mx = 0.5
        self._my = 0.5
        self._step = True
        self._itick = 0
        self._savex = []
        self._savey = []
        
    def Setup(self):

        @self._window.event
        def on_draw():
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            
            glDisable(GL_DEPTH_TEST)
            glLoadIdentity()
            glViewport(0, 0, self._window.width, self._window.height)
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            glOrtho(0,1,0,1,-1,2)
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            
            glColor4f(1.0,1.0,1.0,0.5)
            #nx = hydro.ncells
            nx = 100000 # len(self._savex)
            sx = np.array(self._savex)
            sy = np.array(self._savey)
            try:
                check = sx.min() != sx.max()
            except:
                check = False
            if check:
                pts = np.zeros((nx,2),dtype="float64")
                #pts[:,1] = hydro.vel[0:nx]*1e6
                #pts[:,0] = hydro.x
                #pts[:,1] = np.log10(hydro.rho)
                #print len(self._savex)
                lx = len(sx)
                pts[0:lx,0] = sx
                pts[0:lx,1] = sy
                pts[lx:,0] = pts[lx-1,0]
                pts[lx:,1] = pts[lx-1,1]
                #print pts[:,0].max(), pts[:,1].max()
                #print pts[:,0]

                glMatrixMode(GL_PROJECTION)
                glLoadIdentity()
                #float(vhone.data.xmin),float(vhone.data.xmax
                # No idea what's wrong here if I don't do this... gah
                xmin = pts[:,0].min()-1e-7
                xmax = pts[:,0].max()+1e-7
                ymin = pts[:,1].min()-1e-7
                ymax = pts[:,1].max()+1e-7
                glOrtho(xmin, xmax,
                        ymin, ymax,
                        -2, 2)
                glMatrixMode(GL_MODELVIEW)

                #print "rho MAX/MIN:", vhone.data.zro.min(),vhone.data.zro.max()
                #print "P MAX/MIN:", vhone.data.zpr.min(),vhone.data.zpr.max()
                #print vhone.data.dt
                #pts = rampy.data.xp
                pts = pts.flatten()
                self._glPoints = (GLdouble * len(pts))(*pts)
                
                glEnableClientState(GL_VERTEX_ARRAY)
                glVertexPointer(2, GL_DOUBLE, 0, self._glPoints)
                glDrawArrays(GL_LINE_STRIP,0,len(pts)//2)
            
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
        hydro.P[0] *= 1e8
        # Set up rendering
        s = 0.02
        #glOrtho(-s,+s,-s,+s,-100*s,+100*s)
        glEnable(GL_BLEND)
        glPointSize(2)
        # Set up pyglet to run
        pyglet.clock.set_fps_limit(60)
        pyglet.clock.schedule_interval(self.Step,1.0/60.0)
        pyglet.app.run()
        
    def Step(self, dt):
        global sn
        nx = hydro.ncells
        if self._itick == 0:
            T = hydro.T
            v = hydro.vel
            print "VMAX/MIN:",v.min(),v.max(), hydro.time/init.year
        if self._step:
            self._itick += 1
            # Try a simple wind test
            #if not sn:
            #    vhone.data.zux[0,0,0] += 1e6*vhone.data.dt
            #    vhone.data.zro[0,0,0] += 1e2*vhone.data.dt
            hydro.Step()
            rho = hydro.rho[0:nx]
            rs = hydro.x[np.where(rho == rho.max())]
            if (rs > 0):
                self._savex.append(hydro.time+0.0)
                self._savey.append(rs[0])
            #print hydro.time, rs[0]
            #if not sn and hydro.time > 0.00:
            #    print "BANG"
            #    hydro.P[0] *= 1e8
            #    sn = True
            #    print hydro.time, hydro.dt
        if self._itick > 10:
            self._itick = 0
        


if __name__=="__main__":
    intg = Integrator()
    intg.Setup()
    
