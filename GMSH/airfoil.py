# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 15:36:49 2020

@author: MC
"""
import numpy as np

def camber_line( x, m, p, c ):
    """
    Computes the y-coordinate of the camber line.
    
    Input Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.
    m: float
        maximum camber.
    p: float
        location of maximum camber
    c: float
        chord length
    
    Returns
    -------
    y-coordinates: 1D numpy array of floats.
        y-coordinates of the camber line.
    """
    return np.where((x>=0)&(x<=(c*p)),
                    m * (x / np.power(p,2)) * (2.0 * p - (x / c)),
                    m * ((c - x) / np.power(1-p,2)) * (1.0 + (x / c) - 2.0 * p ))

def dyc_over_dx( x, m, p, c ):
    """
    Computers the slope used to find the perpindicular angle for the cambered line.
    
    Input Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.
    m: float
        maximum camber.
    p: float
        location of maximum camber
    c: float
        chord length
    
    Returns
    -------
    dyc_dx: 1D numpy array of floats.
        slope of the camber line.
    """
    return np.where((x>=0)&(x<=(c*p)),
                    ((2.0 * m) / np.power(p,2)) * (p - x / c),
                    ((2.0 * m ) / np.power(1-p,2)) * (p - x / c ))

def thickness( x, t, c, a):
    """
    Computes the y-coordinates or thickness of one side of a symmetric or
    non-cambered airfoil.
    
    Input Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.
    t: float
        maximum thickness as a fraction of chord.
    c: float
        chord length.
    a: boolean
        Determines if trialing edge of airfoil is closed or open.
        (0.1015) for a non-zero trailing edge thickness and (0.1036) for
        closed trailing edge.
        
    Returns
    -------
    yt: 1D numpy array of floats
        y-coordinates or thickness of airfoil for with no modification.
    """
    if a:
        a4_coeff = 0.1036 # closed trailing edge
    else:
        a4_coeff = 0.1015 # open trailing edge
    
    a0 =  0.2969 * (np.sqrt(x/c))
    a1 = -0.1260 * (x/c)
    a2 = -0.3516 * np.power(x/c,2)
    a3 =  0.2843 * np.power(x/c,3)
    a4 = -a4_coeff * np.power(x/c,4)
    
    return 5 * t * c * (a0 + a1 + a2 + a3 + a4)

def naca4(x, t, m, p, a4, c):
    """
    Computes the x and y coordinates of a naca4 airfoil.
    
    Input Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.    
    t: float
        maximum thickness as a fraction of chord.    
    c: float
        chord length.
    m: float (default 0.0)
        maximum camber.
    p: float (default = 0.3)
        location of maximum camber.
    a4: boolean
        Determines if trailing edge of airfoil is closed or not (default = True).
        True = closed trailing edge, False is open trailing edge.
    
    Returns
    -------
    X, Y: 1D numpy array of floats.
        x and y-coordinates for upper and bottom part of airfoil.
    """
    dyc_dx = dyc_over_dx(x, m, p, c)
    th = np.arctan(dyc_dx)
    yt = thickness(x, t, c, a4)
    yc = camber_line(x, m, p, c)
    
    # calculate the upper and lower x and y values respectively
    xu = x-yt*np.sin(th)
    yu = yc + yt*np.cos(th)
    xl = x + yt*np.sin(th)
    yl = yc - yt*np.cos(th)
       
    return np.asarray([xu, xl]), np.asarray([yu, yl])
    
def naca4_modified(x, t, p, d0, c):
    """
    Computes the x and y coordinates of a modified naca4 airfoil.
    
    Input Parameters
    ----------
    x: 1D array of floats
        x-coordinate of the points defining the geometry.
    m: float
        location of maximum thickness.
    t: float
        maximum thickness as a fraction of chord.
    c: float
        chord length.
    d0: float
        thickness at the trailing edge for an open trailing edge airfoil.
    
    Returns
    -------
    X, Y: 1D numpy array of floats.
        x and y-coordinates for upper and bottom part of airfoil. The format supports
        the contour functions.
    """
    a, d = naca4Coefficients(t, p, d0, c) # solve for a and d coefficients

    xLead = x[x<=p] # x values for leading edge
    xTrail = x[x>p] # x values for trialing edge
    # calculate the thickness
    yLead = a[0]*np.sqrt(xLead/c) + a[1]*(xLead/c) + a[2]*np.power(xLead/c, 2) + a[3]*np.power(xLead/c, 3)
    yTrail = d[0] + d[1]*(1-(xTrail/c)) + d[2]*np.power(1-(xTrail/c), 2) + d[3]*np.power(1-(xTrail/c), 3)
    
    X = np.append(xLead, xTrail)
    Y = np.append(yLead, yTrail)
    
    return np.asarray([X, X]), np.asarray([Y, -Y])

def naca4Coefficients(t, p, d0, c):
    """
    Computes the x and y coordinates of a modified naca4 airfoil.
    
    Input Parameters
    ----------
    m: float
        maximum camber.
    c: float
        chord length.
    d0: float
        thickness at the trailing edge for an open trailing edge airfoil.
    
    Returns
    -------
    a, d: 1d numpy array of floats
        a is coefficients for leading edge of airfoil and d is coefficients of
        trailing edge of airfoil.
    """
    m = p/c
    LEindex = 6
    # Trailing edge angle needed to calculate d1
    TREA = np.poly1d([15.83333333, -2.17857143, -0.24047619, 1.009])
    d1 = TREA(m)*t
    
    # a0 is calculated using sqrt(2*r_t) where r_t is 1.1019*(t*LEindex/6)^2
    a0 = np.sqrt(2*1.1019*(t*LEindex/6)**2)
    
    # system of equations to get d coefficients
    dA = np.array([[(1-m)**2, (1-m)**3], [-2*(1-m), -3*(1-m)**2]])
    dB = np.array([(t/2)-(d1*(1-m))-d0, d1])
    [d2,d3] = np.linalg.solve(dA,dB)
    
    # system of equations to get a coefficients
    aA = np.array([[m, m**2, m**3],[1, 2*m, 3*m**2],[0, 2, 6*m]])
    aB = np.array([(t/2)-a0*m**0.5, -a0/(2*m**0.5), 2*d2+6*d3*(1-m)+(1/4)*a0*m**(-3/2)])
    [a1, a2, a3] = np.linalg.solve(aA, aB)

    return np.array([a0, a1, a2, a3]), np.array([d0, d1, d2, d3])

def convert_to_coordinates(X, Y):
    """
    Computes the x and y coordinates of a modified naca4 airfoil.
    
    Input Parameters
    ----------
    x: 2D numpy array.
        X-coordinates of airfoil.
    y: 2D numpy array.
        Y-coordinates of airfoil.
    
    Returns
    -------
    X, Y: 1D array of floats.
        X and Y are in the format of coordinates for a .dat file.  
        X starts at 1, goes to 0, and then back to one. Y is the corresponding y-coordinates.
    """
    xu = X[0]
    xl = X[1]
    yu = Y[0]
    yl = Y[1]

    # format the output to be inline with a typical .dat file
    x_reverse = np.flipud(xu)
    y_reverse = np.flipud(yu)
    # build array with last coordinate taken out
    X = np.append(x_reverse, xl[1::])
    Y = np.append(y_reverse, yl[1::])
        
    return X, Y