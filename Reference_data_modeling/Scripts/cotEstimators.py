import numpy as np

def phillips_pp(SR, mass, efficiency, kinVisc, fluidDensity, velocity):
    c = np.power(np.divide(np.multiply(6, np.power(SR,2)), np.multiply(fluidDensity, np.pi)), 1/3)
    d = np.power(np.divide(1,fluidDensity),2/3)*(-0.0122*np.power(SR,2)+0.5196*SR+4.2732)
    k = 1 + (1.5*np.power(SR, -3/2)) + (7*np.power(SR, -3))
    length = np.multiply(c, np.power(mass, 1/3))
    Re = np.divide(np.multiply(length, velocity), kinVisc)
    velLow = velocity[np.where(Re<500000)[0]]
    velHigh = velocity[np.where(Re>1000000)[0]]
    velMid = velocity[np.where(np.logical_and(Re>500000, Re<1000000))]
    
    # calculate the low propulsion power
    alpha = 1.327
    beta = -1/2
    e = fluidDensity/(np.power(kinVisc,beta)*2)*np.multiply(np.multiply(k, c),d)
    propLow = (alpha/efficiency)*np.multiply(np.multiply(e,np.power(mass,(beta+2)/3)),np.power(velLow, 3+beta))
    # calculate the high propulsion power
    alpha = 0.072
    beta = -1/5
    e = fluidDensity/(np.power(kinVisc,beta)*2)*np.multiply(np.multiply(k, c),d)
    propHigh = (alpha/efficiency)*np.multiply(np.multiply(e,np.power(mass,(beta+2)/3)),np.power(velHigh, 3+beta))
    
    ppFit = np.polyfit(np.concatenate((velLow, velHigh), axis=0), np.concatenate((propLow, propHigh), axis=0), 6)
    ppLine = np.poly1d(ppFit)
    propPower = ppLine(velocity)
        
    return propPower

def adapted_pp(SR, As, length, mass, efficiency, kinVisc, fluidDensity, velocity):
    k = 1 + (1.5*np.power(SR, -3/2)) + (7*np.power(SR, -3))
    Re = np.divide(np.multiply(length, velocity), kinVisc)
    indLow = np.where(Re<500000)[0]
    indHigh = np.where(Re>1000000)[0]
    indMid = np.where(np.logical_and(Re>500000, Re<1000000))
    ReLow = Re[indLow]
    ReHigh = Re[indHigh]
    ReMid = Re[indMid]
    velLow = velocity[indLow]
    velHigh = velocity[indHigh]
    velMid = velocity[indMid]
    
    # calculate the low propulsion power
    alpha = 1.327
    beta = -1/2
    CdLow = k*alpha*np.power(ReLow, beta)
    #CdLow = CdLow*extraDrag
    propLow = (fluidDensity/(2*efficiency))*As*np.multiply(CdLow, np.power(velLow,3))
    # calculate the high propulsion power
    alpha = 0.072
    beta = -1/5
    CdHigh = k*alpha*np.power(ReHigh, beta)
    #CdHigh = CdHigh*extraDrag
    propHigh = (fluidDensity/(2*efficiency))*As*np.multiply(CdHigh, np.power(velHigh,3))
    # calculate the transition propulsion power
    CdMid = k*((0.074*np.power(ReMid, -1/5))-(1700*np.power(ReMid, -1)))
    #CdMid = CdMid*extraDrag
    propMid = (fluidDensity/(2*efficiency))*As*np.multiply(CdMid, np.power(velMid, 3))
    
    propPower = np.concatenate((propLow, propMid, propHigh), axis=0)       
    
    return propPower

def ITTC_pp(SR, As, length, efficiency, kinVisc, fluidDensity, velocity):
    k = 1 + (1.5*np.power(SR, -3/2)) + (7*np.power(SR, -3))
    Re = np.divide(np.multiply(length, velocity), kinVisc)
    logRe = np.log10(Re)
    Cf = np.divide(0.075, np.power(logRe-2,2))
    Cd = Cf*k
    #Cd = Cd*extraDrag
    
    propPower = (fluidDensity/(2*efficiency))*As*np.multiply(Cd, np.power(velocity,3)) 
            
    return propPower