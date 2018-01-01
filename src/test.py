'''
Created on 3 Jun 2014

@author: samgeen
'''

import vhone
import init
import numpy as np

import pyglet
from pyglet.gl import *

from init import hydro

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
            nx = vhone.data.imax
            pts = np.zeros((nx,2),dtype="float64")
            #import pdb; pdb.set_trace()
            pts[:,0] = vhone.data.zxa
            #pts[:,1] = np.log10(vhone.data.zpr[0:nx,0,0] \
            #                    /vhone.data.zro[0:nx,0,0])
            T = vhone.data.zpr[0:nx,0,0]/vhone.data.zro[0:nx,0,0]/init.kB
            pts[:,1] = np.log10(vhone.data.zpr[0:nx,0,0])
            #pts[:,1] = np.log10(T)


            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            #float(vhone.data.xmin),float(vhone.data.xmax
            # No idea what's wrong here if I don't do this... gah
            xmin = pts[:,0].min()
            xmax = pts[:,0].max()
            ymin = pts[:,1].min()-0.5
            ymax = pts[:,1].max()+0.5
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
        if self._itick == 0:
            nx = vhone.data.imax
            T = vhone.data.zpr[0:nx,0,0]/vhone.data.zro[0:nx,0,0]/init.kB
            print "DMAX/MIN:",T.min(),T.max(), vhone.data.time/init.year
        if self._step:
            self._itick += 1
            # Try a simple wind test
            #if not sn:
            #    vhone.data.zux[0,0,0] += 1e6*vhone.data.dt
            #    vhone.data.zro[0,0,0] += 1e2*vhone.data.dt
            vhone.data.step()
            if not sn and vhone.data.time > 0.00:
                print "BANG"
                vhone.data.zpr[0,0,0] *= 1e8
                sn = True
                print vhone.data.time, vhone.data.dt
        if self._itick > 10:
            self._itick = 0
        


if __name__=="__main__":
    intg = Integrator()
    intg.Setup()
    
