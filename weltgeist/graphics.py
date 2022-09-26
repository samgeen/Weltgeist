# -*- coding: utf-8 -*-
"""
Graphics module

This module contains code for visualising the output of the simulation
 in real-time. This is useful for diagnostics and to get a feel for how
 the setup works. Use matplotlib/pyplot to produce "official" figures

Created on 8 Feb 2019

@author: samgeen
"""

import numpy as np

import pyglet
from pyglet.gl import *

red = [1.0,0.0,0.0,1.0]
blue = [0.0,0.0,1.0,1.0]
black = [0.0,0.0,0.0,1.0]

WINSIZEX = 512 
WINSIZEY = 512

def _newtext(x,y):
    """
    Just make a basic text object
    """
    return  pyglet.text.Label('',
                font_name='Helvetica',
                font_size=10,
                x=x*WINSIZEX,y=y*WINSIZEY,
                anchor_x='left', anchor_y='bottom',color=(0,0,0,255))

class Line(object):
    """
    A line object, used to draw the simulation state to screen in
     real time
    """
    def __init__(self,colour,width=1.0):
        """
        Constructor

        Parameters
        ----------

        colour : list 
            The colour - 4 elements - red, green, blue, alpha
            e.g. [1.0, 0.0, 0.0, 0.5] = red, half transparent
        width : float
            width of the line on-screen in pixels
        """
        self._colour = colour
        self._currlen = 0
        self._x = np.zeros(100)
        self._y = np.zeros(100)
        self._pts = np.zeros((100,2))
        self._ymin = None
        self._ymax = None
        self._width = width

    def Update(self,x,y):
        """
        Replace arrays in the line completely

        Parameters
        ----------

        x, y : numpy array
            numpy arrays containing x and y data for the line
            NOTE: x and y must be the same length
        """
        newlen = len(x)
        self._currlen = newlen
        if newlen > self._x.shape[0]:
            self._x = np.zeros(2*newlen)
            self._y = np.zeros(2*newlen)
            self._pts = np.zeros((2*newlen,2))
        self._x[0:newlen] = x
        self._y[0:newlen] = y

    def Append(self,x,y):
        """
        Add a single point

        Parameters
        ----------

        x, y : float
            points to add to the line
        """
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
        """
        Return min/max in arrays

        Parameters
        ----------

        yminin, ymaxin : float
            (optional) fixes y values for plotting
            Useful if you want to stop the line moving up and down
              over time

        Returns
        -------

        (xmin, xmax, ymin, ymax): tuple of floats
            min/max values in x and y
        """
        if yminin is not None:
            self._ymin = yminin
        if ymaxin is not None:
            self._ymax = ymaxin
        # Get non-empty x and y values
        x = self._x[0:self._currlen]
        y = self._y[0:self._currlen]
        # If no data, return zeros
        if len(x) == 0:
            return 0.0,0.0,0.0,0.0
        else:
            # Set ymin, ymax to chosen values if given
            if self._ymin is not None:
                ymin = self._ymin
            else:
                ymin = y.min()
            if self._ymax is not None:
                ymax = self._ymax
            else:
                ymax = y.max()
            # Return min and max values in x and y
            return x.min(), x.max(), ymin, ymax

    def Draw(self,xmin,xmax,ymin,ymax):
        """
        Draw line to screen
        Called by the Lines object Draw function
        
        Parameters
        ----------

        xmin, xmax, ymin, ymax : float
            Gives limits of the screen in x and y
            Lines scaled to these values
        """
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
    """
    A container for lines you want to draw to screen
    """
    def __init__(self,lines,width=1.0):
        """
        Constructor

        Parameters
        ----------

        lines : list (or other iterable) of Line objects
            The lines you want to draw 
            (external function should set the data points)
        width : float
            width of the lines on-screen in pixels
        """
        self._lines = lines
        self._xmin = 1e30
        self._xmax = -1e30
        self._ymin = 1e30
        self._ymax = -1e30
        # Text on screen for max/min labels
        self._xmintext = _newtext(0.05,0.01)
        self._xmaxtext = _newtext(0.97,0.01)
        self._ymintext = _newtext(0.01,0.05)
        self._ymaxtext = _newtext(0.01,0.97)
        self._xmaxtext.anchor_x='right'

    def _FindMinMax(self):
        """
        Finds the min/max of the lines for clamping to the screen
        (private method, don't call outside)
        """
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
        """
        Draws the lines, called by the renderer
        """
        self._FindMinMax()
        for line in self._lines:
            line.Draw(self._xmin,self._xmax,self._ymin,self._ymax)
    
    def DrawAxes(self):
        # Draw axis markers
        self._xmintext.text = "%.2g" % self._xmin
        self._xmaxtext.text = "%.2g" % self._xmax
        self._ymintext.text = "%.2g" % self._ymin
        self._ymaxtext.text = "%.2g" % self._ymax
        self._xmintext.draw()
        self._xmaxtext.draw()
        self._ymintext.draw()
        self._ymaxtext.draw()

class Renderer(object):
    """
    Controls the basic line drawing interface.
    This uses pyglet, which is a module for drawing to the screen with
     OpenGL
    WARNING! Due to the way pyglet works, you need to feed it the
     step function from the integrator, and pyglet will control the
     stepping at a maximum of 120 frames per second
     See example scripts for an example of this
    """
    def __init__(self, lines):
        """
        Constructor

        Parameters
        ----------

        lines : list (or other iterable) of Line objects
            The lines you want to draw 
            (external function should set the data points)
        """
        self._glPoints = None
        # Pyglet window. This controls most of the pyglet drawing
        self._window = pyglet.window.Window(WINSIZEX,WINSIZEY)
        # Mouse position - this isn't really used, but is left in
        #  in case you want to use it
        self._mx = 0.5
        self._my = 0.5
        # Lines object
        self._lines = Lines(lines)
        # Text on screen
        self._text = pyglet.text.Label('',
                          font_name='Helvetica',
                          font_size=16,
                          x=self._window.width - 10, y=self._window.height - 10,
                          anchor_x='right', anchor_y='top',color=(0,0,0,255))
        # Pause the simulation?
        self._pause = False

    def Text(self,newtext):
        self._text.text = newtext

    def Start(self, stepfunc):
        """
        This will set up and start the renderer

        Parameters
        ----------

        stepfunc : function
            function to run every graphical step
        """
        # Graphical setup stuff
        # Mostly just simplifying OpenGL to display in basic 2D
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
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            glOrtho(0, self._window.width, 
                    0, self._window.height,
                    -2, 2)
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            self._lines.DrawAxes()
            self._text.draw()
        
        stepInterval = 1.0/120.0
            
        # Mouse movement controls
        # This isn't really used, but kept in in case it's useful later
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
            
        # Mouse button press
        # Also not used here, but could be in future
        @self._window.event
        def on_mouse_release(x, y, button, modifiers):
            pass

        @self._window.event
        def on_key_press(symbol, modifiers):
            if symbol == pyglet.window.key.SPACE:
                if self._pause:
                    pyglet.clock.schedule_interval(stepfunc,stepInterval)
                    self._pause = False
                else:
                    pyglet.clock.unschedule(stepfunc)
                    self._pause = True

        @self._window.event
        def on_move(x,y):
            # Screen goes black on moving under WSL
            # Might just be a WSL bug?
            print ("MOVING", x, y)
            #self._window.set_location(x, y)
            #self._window.clear()

        @self._window.event
        def on_resize(w,h):
            print ("RESIZING", w, h)

        @self._window.event
        def on_context_lost():
            print ("CONTEXT LOST!!!")

        @self._window.event
        def on_context_state_lost():
            print ("CONTEXT STATE LOST!!!")

        # Set up rendering
        glEnable(GL_BLEND)
        glPointSize(2)

        # Set up pyglet to run at up to 120 frames every second
        pyglet.clock.schedule_interval(stepfunc,stepInterval)
        pyglet.app.run()