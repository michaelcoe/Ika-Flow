# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:51:46 2020

@author: MC
"""
import numpy as np

class Carangiform:
    """
    Carangiform Motion Class.
    """
    def __init__(self, w, k, A_n, A_t, x_pivot, length):
        """
        Initiallizes a carangiform motion object.
        
        Input Parameters
        ---------
        x: 1d numpy array.
            x-coordinates of the midline.
        w: float.
            oscillation frequency.
        k: float.
            body wave number.
        A_n: float. 
            amplitude of the nose.
        A_t: float.
            amplitude of the tail
        x_pivot: float.
            position along the body (in percent body_lengths) of the pivot.
        length: float.
            length of the fish / robot.
        """
        
        # solving for the coefficients
        X = np.array([[1, 0, 0],[1, x_pivot, x_pivot**2],
                      [1, length, length**2]])
        A = np.array([[A_n], [0], [A_t]])
        C = np.linalg.solve(X, A)
        
        self.c0 = C[0]
        self.c1 = C[1]
        self.c2 = C[2]
        self.w = w
        self.k = k
        self.x_pivot = x_pivot
        self.length = length
        
    def rigid_motion(self, x, t):       
        """
        Computes the midline curvature as a function of time
        
        Input Parameters
        ----------
        x: 1d numpy array
            x-coordinates of midline.
        t: float.
            the timestep in seconds.
        
        Returns
        -------
        h(x,t): 1d numpy array.
            h-coordinates of the midline.
        """
        # mask out the rigid part from flexible part
        indexes = np.where(x <= self.x_pivot)[0]
        mask = np.ones(x.size, dtype=bool)
        mask[indexes] = False
        x_leading = x[indexes]
        x_trailing = x[mask]
        
        # calculate the rigid leading edge from flexible trialing edge
        h_leading = (-self.c0/self.x_pivot*np.sin(self.w*t)*x_leading + 
                     self.c0*np.sin(self.w*t))
        h_trailing = ((self.c0 + self.c1*x_trailing + self.c2*x_trailing**2) *
                      np.sin(self.w*t-self.k*x_trailing))
        
        return np.append(h_leading, h_trailing)
        
    def single_rigid_motion(self, x, t):       
        """
        Computes the midline curvature as a function of time 
        for a single value
        
        Input Parameters
        ----------
        x: 1D numpy array of panel x-coordinates.
            x-coordinates of midline.
        t: float.
            the timestep in seconds.
        
        Returns
        -------
        h(x,t): 1d numpy array.
            h-coordinates of the midline.
        """
        if x[-1] <= self.x_pivot:
            return (-self.c0/self.x_pivot*np.sin(self.w*t)*x + 
                     self.c0*np.sin(self.w*t))
        elif x[0] <= self.x_pivot and x[-1] > self.x_pivot:
            return np.append(-self.c0/self.x_pivot*np.sin(self.w*t)*x[0] + 
                     self.c0*np.sin(self.w*t), (self.c0 + self.c1*x[-1] + 
                     self.c2*x[-1]**2) * np.sin(self.w*t-self.k*x[-1]))
        else:
            return ((self.c0 + self.c1*x + self.c2*x**2) *
                      np.sin(self.w*t-self.k*x))
        
    