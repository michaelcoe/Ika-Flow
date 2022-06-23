# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 12:47:13 2020

@author: MC
"""
import numpy as np

class Freestream:
    """
    Freestream conditions.
    """
    def __init__(self, u_inf=1.0, alpha=0.0):
        """
        Sets the freestream speed and angle (in degrees).
        
        Input Parameters
        ----------
        u_inf: float, optional
            Freestream speed;
            default: 1.0.
        alpha: float, optional
            Angle of attack in degrees;
            default 0.0.
        """
        self.u_inf = u_inf
        self.alpha = np.radians(alpha)  # degrees to radians

class Panel:
    """
    Contains information related to a panel.
    """
    def __init__(self, xa, ya, xb, yb):
        """
        Initializes the panel.
        
        Sets the end-points and calculates the center-point, length,
        and angle (with the x-axis) of the panel.
        Defines if the panel is located on the upper or lower surface of the geometry.
        Initializes the source-strength, tangential velocity, and pressure coefficient
        of the panel to zero.
        
        Input Parameters
        ---------_
        xa: float
            x-coordinate of the first end-point.
        ya: float
            y-coordinate of the first end-point.
        xb: float
            x-coordinate of the second end-point.
        yb: float
            y-coordinate of the second end-point.
        """
        self.xa, self.ya = xa, ya  # panel starting-point
        self.xb, self.yb = xb, yb  # panel ending-point
        
        self.x_bound, self.y_bound = (xa + xb) / 4, (ya + yb) / 4 # panel bound vortex
        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2 # panel center
        self.xcp, self.ycp = (3/2)*(xa + xb), (3/2)*(ya + yb) # panel control point
        self.length = np.sqrt((xb - xa)**2 + (yb - ya)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb - xa <= 0.0:
            self.beta = np.arccos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = np.pi + np.arccos(-(yb - ya) / self.length)
        
        # panel location
        if self.beta <= np.pi:
            self.loc = 'upper'  # upper surface
        else:
            self.loc = 'lower'  # lower surface
        
        self.gamma = 0.0  # vortex strength
        self.vt = 0.0  # tangential velocity
        self.cp = 0.0  # pressure coefficient
        
    def update_position(self, x, y):
        """
        Allows for the update of the x and y position of the panel
        
        Sets the end-points and calculates the center-point, length,
        and angle (with the x-axis) of the panel.
        
        Input Parameters
        ---------
        x: 1D numpy array of x-coordinates.
            New x-coordinates of the panel.
        y: 1D numpy array of y-coordinates.
            New y-coordinates of the panel.
        """
        # define the new values to work with
        xa_new, xb_new = self.xa + x[0], self.xb + x[1]
        ya_new, yb_new = self.ya + y[0], self.yb + y[1]
        self.xa, self.ya = xa_new, ya_new
        self.xb, self.yb = xb_new, yb_new
        
        self.x_bound, self.y_bound = (xa_new + xb_new) / 4, (ya_new + yb_new) / 4 # panel bound vortex
        self.xcp, self.ycp = (3/2)*(xa_new + xb_new), (3/2)*(ya_new + yb_new) # panel control point
        self.xc, self.yc = (xa_new + xb_new) / 2, (ya_new + yb_new) / 2  # panel center
        self.length = np.sqrt((xb_new - xa_new)**2 + (yb_new - ya_new)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb_new - xa_new <= 0.0:
            self.beta = np.arccos((yb_new - ya_new) / self.length)
        elif xb_new - xa_new > 0.0:
            self.beta = np.pi + np.arccos(-(yb_new - ya_new) / self.length)
            
class Wake_panel:
    """
    Contains information related to the wake.
    """
    def __init__(self, xa, ya, xb, yb):
        """
        Initializes a wake panel.
            
        Input Parameters
        ---------_
        gamma: float
            vortex strength of a panel (default = 0).
        U: float
            velocity in the x-direction.
        dt: float
            time step.
        """
        self.xa, self.ya = xa, ya        # set the start of the panel
        self.xb, self.yb = xb, yb # set the the other end of the panel
        
        self.x_bound, self.y_bound = (self.xa + self.xb) / 4, (self.ya + self.yb) /4
        self.xcp, self.ycp = (xa + xb) / 2, (ya + yb) / 2 # panel center
        self.xc, self.yc = (xa + xb) / 2, (ya + yb) / 2  # panel center
        self.length = np.sqrt((xb - xa)**2 + (yb - ya)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb - xa <= 0.0:
            self.beta = np.arccos((yb - ya) / self.length)
        elif xb - xa > 0.0:
            self.beta = np.pi + np.arccos(-(yb - ya) / self.length)
        
        self.gamma = 0.0  # vortex strength
        self.vt = 0.0  # tangential velocity
        self.cp = 0.0  # pressure coefficient
        
    def update_position(self, x, y):
        """
        Allows for the update of the x and y position of the panel
        
        Sets the end-points and calculates the center-point, length,
        and angle (with the x-axis) of the panel.
            
        Input Parameters
        ---------
        x: 1D numpy array of x-coordinates.
            New x-coordinates of the panel.
        y: 1D numpy array of y-coordinates.
            New y-coordinates of the panel.
        """
        # define the new values to work with
        xa_new, xb_new = x[0], x[1]
        ya_new, yb_new = y[0], y[1]
        self.xa, self.ya = xa_new, ya_new
        self.xb, self.yb = xb_new, yb_new
        
        self.x_bound, self.y_bound = (xa_new + xb_new) / 4, (ya_new + yb_new) / 4 # panel bound vortex
        self.xcp, self.ycp = (3/2)*(xa_new + xb_new), (3/2)*(ya_new + yb_new) # panel control point
        self.xc, self.yc = (xa_new + xb_new) / 2, (ya_new + yb_new) / 2  # panel center
        self.length = np.sqrt((xb_new - xa_new)**2 + (yb_new - ya_new)**2)  # panel length
        
        # orientation of panel (angle between x-axis and panel's normal)
        if xb_new - xa_new <= 0.0:
            self.beta = np.arccos((yb_new - ya_new) / self.length)
        elif xb_new - xa_new > 0.0:
            self.beta = np.pi + np.arccos(-(yb_new - ya_new) / self.length)
    