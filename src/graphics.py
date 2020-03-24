'''
Created on 8 Feb 2019

@author: samgeen
'''

import numpy as np

import pyglet
from pyglet.gl import *


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

class Renderer(object):
    '''
    Controls the basic line interface
    '''
    def __init__(self, lines):
        self._glPoints = None
        self._window = pyglet.window.Window(512,512)
        self._window._integrator = self
        self._mx = 0.5
        self._my = 0.5
        self._lines = Lines(lines)

    def Start(self, stepfunc):
        '''
        This will set up and start the renderer
        stepfunc: function to run every graphical step
        '''
        # Graphical setup stuff
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

        # Set up rendering
        s = 0.02
        glEnable(GL_BLEND)
        glPointSize(2)
        # Set up pyglet to run
        #pyglet.clock.set_fps_limit(60)
        pyglet.clock.schedule_interval(stepfunc,1.0/60.0)
        pyglet.app.run()