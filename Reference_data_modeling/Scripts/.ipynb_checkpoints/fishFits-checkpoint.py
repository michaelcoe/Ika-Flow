# in airfoil.py
import cv2 as cv2
import numpy as np
import Scripts.airfoil as af
from scipy import integrate

#----- Method List -----------------------
# fit_side_contours(topHull, bottomHull, degree)
# fit_top_contours(cntsArea, p, t, a4)
# fit_top_contours_fixed(cntsArea, p, t, a4)
# fit_top_contours_poly(topHull, bottomHull, degree)
#----------------------------------------

def fit_side_contours(topHull, bottomHull, degree):
    topFit = np.polyfit(topHull[0,:], topHull[1,:], degree)
    bottomFit = np.polyfit(bottomHull[0,:], bottomHull[1,:], degree)
    
    return topFit, bottomFit

def fit_top_contours(cntsArea, m, t, d0):
    # naca profile parameters
    c = 1.0
    
    x = np.linspace(0,1,200)
    X, Y = af.naca4_modified(x, m, t, c, d0)
    
    xu = np.asarray(X[0])
    xl = np.asarray(X[1])
    yu = np.asarray(Y[0])
    yl = np.asarray(Y[1])
    nacaArea = 2*integrate.simps(yu, xu)
    
    convergence = 0.001
    
    if(np.abs(np.round(nacaArea, 3) - np.round(cntsArea,3)) <= convergence):
        return nacaArea, m, t, d0, xu, yu, xl, yl
    elif(np.round(nacaArea, 3) - np.round(cntsArea,3) > convergence):
        return (fit_top_contours(cntsArea, m, t-0.001, d0))
    else:
        return (fit_top_contours(cntsArea, m, t+0.001, d0))
    
def fit_top_contours_fixed(cntsArea, p, t, a4):
    # naca profile parameters
    m = 0.0
    c = 1.0
    
    x = np.linspace(0,1,200)
    X, Y = af.naca4(x, m, p, t, c, 0.1015-a4)
    
    xu = np.asarray(X[0])
    xl = np.asarray(X[1])
    yu = np.asarray(Y[0])
    yl = np.asarray(Y[1])
    nacaArea = 2*integrate.simps(yu, xu)

    return nacaArea, t, 0.1015-a4, xu, yu, xl, yl

def fit_top_contours_poly(topHull, bottomHull, degree):
    topFit = np.polyfit(topHull[0,:], topHull[1,:], degree)
    bottomFit = np.polyfit(bottomHull[0,:], bottomHull[1,:], degree)
    
    topArea = integrate.simps(topHull[1,:], topHull[0,:])
    bottomArea = integrate.simps(bottomHull[1,:], bottomHull[0,:])
    
    polyArea = np.abs(topArea) + np.abs(bottomArea)
    
    return polyArea, topFit, bottomFit