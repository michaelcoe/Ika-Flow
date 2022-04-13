import math
import matplotlib.pyplot as plt
import numpy as np

#https://en.wikipedia.org/wiki/NACA_airfoil#Equation_for_a_cambered_4-digit_NACA_airfoil
def camber_line( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    m * (x / np.power(p,2)) * (2.0 * p - (x / c)),
                    m * ((c - x) / np.power(1-p,2)) * (1.0 + (x / c) - 2.0 * p ))

def dyc_over_dx( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    ((2.0 * m) / np.power(p,2)) * (p - x / c),
                    ((2.0 * m ) / np.power(1-p,2)) * (p - x / c ))

def thickness( x, t, c, a):
    a0 =  0.2969 * (np.sqrt(x/c))
    a1 = -0.1260 * (x/c)
    a2 = -0.3516 * np.power(x/c,2)
    a3 =  0.2843 * np.power(x/c,3)
    a4 = -a * np.power(x/c,4)
    return 5 * t * c * (a0 + a1 + a2 + a3 + a4)

def naca4(x, m, p, t, c, a4):
    dyc_dx = dyc_over_dx(x, m, p, c)
    th = np.arctan(dyc_dx)
    yt = thickness(x, t, c, a4)
    yc = camber_line(x, m, p, c)
    
    # calculate the upper and lower x and y values respectively
    X = np.asarray([x-yt*np.sin(th), x + yt*np.sin(th)])
    Y = np.asarray([yc + yt*np.cos(th), yc - yt*np.cos(th)])
    
    return X, Y

def naca4Coefficients(m, t, c, d0):
    m = m/c
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
    
def naca4_modified(x, m, t, c, d0):
    a, d = naca4Coefficients(m, t, c, d0)
    # x values for leading edge
    xLead = x[x<=m]
    # x values for trialing edge
    xTrail = x[x>m]
    # calculate the thickness
    yLead = a[0]*np.sqrt(xLead/c) + a[1]*(xLead/c) + a[2]*np.power(xLead/c, 2) + a[3]*np.power(xLead/c, 3)
    yTrail = d[0] + d[1]*(1-(xTrail/c)) + d[2]*np.power(1-(xTrail/c), 2) + d[3]*np.power(1-(xTrail/c), 3)
    X = np.append(xLead, xTrail)
    Y = np.append(yLead, yTrail)
    
    return np.asarray([X, X]), np.asarray([Y, -Y])