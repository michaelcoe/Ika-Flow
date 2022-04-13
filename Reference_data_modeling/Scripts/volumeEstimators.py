import numpy as np
import Scripts.airfoil as af
from scipy import integrate
from scipy.special import ellipe

#-------------Method List----------------------------
# determine_volume(topForm, crossSection, aspectRatio, length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop)
# volume_angulliform(xu, topPolySide, bottomPolySide, heightWidthRatio, length)
# volume_fusiform(topPolySide, bottomPolySide, xu, yu, length)
# volume_skate(topPolySide, bottomPolySide, xu, yu, length)
# volume_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length)
# volume_teardrop(topPolySide, bottomPolySide, xu, yu, length)
# volume_oval(topPolySide, bottomPolySide, xu, yu, length)
# volume_box(topPolySide, bottomPolySide, xu, yu, length)
# volume_triangle(topPolySide, bottomPolySide, xu, yu, length)
# equivalentSpheroid(length, mass, fluidDensity)
#------------------------------------------------------

def determine_volume(topForm, crossSection, aspectRatio, length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop):
    dx = np.linspace(0.0, 1.0, 100)
    # determine if airfoil or not
    if(topForm):
        topPolySide = np.poly1d(topCoeffSide)
        bottomPolySide = np.poly1d(bottomCoeffSide)
        # Naca Airfoil equation (x, camber height, location of max thickness,
        # max thickness, chord length, last coefficient of equation)
        X, Y = af.naca4_modified(dx, topCoeffTop[0], topCoeffTop[1], 1.0, topCoeffTop[2])
        xu = X[0]
        yu = Y[0]  
    else:
        topPolySide = np.poly1d(topCoeffSide)
        bottomPolySide = np.poly1d(bottomCoeffSide) 
        topPolyTop = np.poly1d(topCoeffTop)
        bottomPolyTop = np.poly1d(bottomCoeffTop)
        xu = dx
        yu= topPolyTop(dx)
    # if else statement for the cross section 
    if(crossSection==1):
        return volume_angulliform(xu, topPolySide, bottomPolySide, aspectRatio, length)
    elif(crossSection==2):
        return volume_fusiform(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==3):
        return volume_skate(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==4):
        return volume_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==5):
        return volume_teardrop(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==6):
        return volume_oval(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==7):
        return volume_box(topPolySide, bottomPolySide, xu, yu, length)
    else:
        return volume_triangle(topPolySide, bottomPolySide, xu, yu, length)

def volume_angulliform(xu, topPolySide, bottomPolySide, heightWidthRatio, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    topPointsTop = np.divide(topPointsSide+bottomPointsSide, heightWidthRatio*2)
    # distance between each circumference
    t = (xu[1]-xu[0])
    # calculate the volume based on circumference and thickness
    topVolume = np.sum(np.pi * t/2 * np.multiply(topPointsTop, topPointsSide))
    bottomVolume = np.sum(np.pi * t/2 * np.multiply(topPointsTop, bottomPointsSide))
    # combine the two volumes  
    return (topVolume + bottomVolume)

def volume_fusiform(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    # surface area sums up the N-2 frustrums
    topPointsSide = np.abs(topPolySide(xu[1:]))*length
    bottomPointsSide = np.abs(bottomPolySide(xu[1:]))*length
    xu = xu[1:]*length
    yu = yu[1:]*length
    height = np.max(np.abs(topPointsSide) + np.abs(bottomPointsSide))
    width = np.max(np.abs(2*yu))
    # distance between each circumference
    t = (xu[1]-xu[0])
    # calculate the volume based on circumference and thickness
    topVolume = np.sum(np.pi * t/2 * np.multiply(yu, topPointsSide))
    bottomVolume = np.sum(np.pi * t/2 * np.multiply(yu, bottomPointsSide))
    # combine the two volumes   
    return (topVolume + bottomVolume), height, width

def volume_skate(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    xu = xu*length
    yu = yu*length
    # Thickness of the frustrums
    t = (xu[1]-xu[0])
    # note that skates are the inverse of the fusilform body
    # calculate the volume based on circumference and thickness
    topVolume = np.sum(np.pi * t/2 * np.multiply(yu, topPointsSide))
    bottomVolume = np.sum(np.pi * t/2 * np.multiply(yu, bottomPointsSide))
    # combine the two volumes 
    return (topVolume + bottomVolume)

def volume_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    xu = xu*length
    yu = yu*length
    # distance between each point
    t = (xu[1]-xu[0])
    # split the dx into the first 2/3 and last 1/3
    firstIndex = int(len(xu)*2/3)
    lastIndex = len(xu)-int(len(xu)*1/3)
    #--------------------------------------------------------
    # Calculate the area of the inverted teardrop portion of the fish (first 2/3)
    # Side of the triangle is sqrt(a^2 + b^2)
    # top portion is an ellipse and bottom portion is a triangle
    topEllipseVolume = np.sum(np.pi * t/2 * np.multiply(yu[0:firstIndex], topPointsSide[0:firstIndex]))
    # cross section is flattened at the bottom, assume that around 5 percent of the triangle side will fold over to make flattened bottom
    bottomTriangleVolume = np.sum(np.pi * t/2 * np.multiply(yu[0:firstIndex], topPointsSide[0:firstIndex]))
    #-------------------------------------------------------
    # calculate the area of the ellipse portion of the fish (last 1/3)
     # calculate the volume based on circumference and thickness
    topTailVolume = np.sum(np.pi * t/2 * np.multiply(yu[lastIndex:-1], topPointsSide[lastIndex:-1]))
    bottomTailVolume = np.sum(np.pi * t/2 * np.multiply(yu[lastIndex:-1], bottomPointsSide[lastIndex:-1]))
         
    return (topEllipseVolume + bottomTriangleVolume + topTailVolume + bottomTailVolume)

def volume_teardrop(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    xu = xu*length
    yu = yu*length
    # distance between each point
    t = (xu[1]-xu[0])
    # split the dx into the first 2/3 and last 1/3
    firstIndex = int(len(xu)*2/3)
    lastIndex = len(xu)-int(len(xu)*1/3)
    #--------------------------------------------------------
    # Calculate the area of the inverted teardrop portion of the fish (first 2/3)
    # Side of the triangle is sqrt(a^2 + b^2)
    # bottom portion is an ellipse and top portion is a triangle
    # cross section is flattened at the bottom, assume that around 5 percent of the triangle side will fold over to make flattened bottom
    topTriangleArea = np.sum(np.pi * t/2 * np.multiply(yu[0:firstIndex], topPointsSide[0:firstIndex]))
    # bottom Ellipse
    bottomEllipseArea = np.sum(np.pi * t/2 * np.multiply(yu[0:firstIndex], bottomPointsSide[0:firstIndex]))
    #-------------------------------------------------------
    # calculate 1/4 of the circumference representing one side of the fish
    topTailArea = np.sum(np.pi * t/2 * np.mulitply(yu[lastIndex:-1], topPointsSide[lastIndex:-1]))
    bottomTailArea= np.sum(np.pi * t/2 * np.multiply(yu[lastIndex:-1], bottomPointsSide[lastIndex:-1]))
    # calculate the surface area based on circumference and thickness
         
    return (bottomEllipseArea+topTriangleArea+topTailArea+bottomTailArea)

def volume_oval(dx, topPolySide, bottomPolySide, topPolyTop, bottomPolyTop, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(dx))*length
    bottomPointsSide = np.abs(bottomPolySide(dx))*length
    topPointsTop = np.abs(topPolyTop(dx))*length
    bottomPointsTop = np.abs(bottomPolyTop(dx))*length
    # distance between each circumference
    t = (dx[1]-dx[0])*length
    # calculate the eccentricity of the top and bottom
    eSquaredTopRight = 1.0 - np.divide(np.power(topPointsTop, 2), np.power(topPointsSide, 2))
    eSquaredBottomRight = 1.0 - np.divide(np.power(topPointsTop, 2), np.power(bottomPointsSide, 2))
    eSquaredTopLeft = 1.0 - np.divide(np.power(bottomPointsTop, 2), np.power(topPointsSide, 2))
    eSquaredBottomLeft = 1.0 - np.divide(np.power(bottomPointsTop, 2), np.power(bottomPointsSide, 2))
    # calculate 1/4 of the circumference representing one side of the fish
    # calculate the surface area based on circumference and thickness
    topAreaRight = np.sum(np.multiply(np.multiply(topPointsSide, ellipe(eSquaredTopRight)),t))
    bottomAreaRight = np.sum(np.multiply(np.multiply(bottomPointsSide, ellipe(eSquaredBottomRight)),t))
    topAreaLeft = np.sum(np.multiply(np.multiply(topPointsSide, ellipe(eSquaredTopLeft)),t))
    bottomAreaLeft = np.sum(np.multiply(np.multiply(bottomPointsSide, ellipe(eSquaredBottomLeft)),t))
    # combine to get surface area of one side and multiple by 2 for total surface area    
    return topAreaRight+bottomAreaRight+topAreaLeft+bottomAreaLeft

def volume_box(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    topPointsTop = np.abs(yu)*length
    # distance between each point
    t = (xu[1]-xu[0])
    # split the dx into the first 2/3 and last 1/3
    firstIndex = int(len(xu)*2/3)
    lastIndex = len(xu)-int(len(xu)*1/3)
    #--------------------------------------------------------
    # Calculate the area of the box portion of the fish (first 2/3)
    topBoxArea = np.sum(2 * t * np.multiply(topPointsTop[0:firstIndex], topPointsSide[0:firstIndex]))
    bottomBoxArea = np.sum(2 * t * np.multiply(topPointsTop[0:firstIndex], bottomPointsSide[0:firstIndex]))
    #-------------------------------------------------------
    # calculate the surface area based on circumference and thickness
    topTailArea = np.sum(np.pi * t * np.multiply(topPointsSide[lastIndex:-1], ellipe(eSquaredTop)))
    bottomTailArea= np.sum(np.pi * t * np.multiply(bottomPointsSide[lastIndex:-1], ellipe(eSquaredBottom)))
         
    return topBoxArea+bottomBoxArea+topTailArea+bottomTailArea

def volume_triangle(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    topPointsTop = np.abs(yu)*length
    midPointsTop = topPointsTop*2/3
    # distance between each point
    t = (xu[1]-xu[0])
    # split the dx into the first 2/3 and last 1/3
    firstIndex = int(len(xu)*2/3)
    lastIndex = len(xu)-int(len(xu)*1/3)
    #--------------------------------------------------------
    # Calculate the area of the triangle portion of the fish (first 2/3)
    # Side of the triangle is sqrt(a^2 + b^2)
    topTriangleArea = np.sum(np.multiply(np.sqrt(np.power(topPointsSide[0:firstIndex],2)+np.power(midPointsTop[0:firstIndex],2)),t))
    bottomTriangleArea = np.sum(np.multiply(topPointsTop[0:firstIndex]+np.sqrt(np.power(bottomPointsSide[0:firstIndex],2)
                                                     +np.power(topPointsTop[0:firstIndex]-midPointsTop[0:firstIndex],2)),t))
    #-------------------------------------------------------
    # calculate the area of the ellipse portion of the fish (last 1/3)
    # calculate the eccentricity of the top and bottom
    eSquaredTop = 1.0 - np.divide(np.power(topPointsTop[lastIndex:-1], 2), np.power(topPointsSide[lastIndex:-1], 2))
    eSquaredBottom = 1.0 - np.divide(np.power(topPointsTop[lastIndex:-1], 2), np.power(bottomPointsSide[lastIndex:-1], 2))
    # calculate 1/4 of the circumference representing one side of the fish
    topTailArea = np.sum(np.multiply(np.multiply(topPointsSide[lastIndex:-1], ellipe(eSquaredTop)),t))
    bottomTailArea= np.sum(np.multiply(np.multiply(bottomPointsSide[lastIndex:-1], ellipe(eSquaredBottom)),t))
         
    return 2.0*(topTriangleArea+bottomTriangleArea+topTailArea+bottomTailArea)

def equivalentSpheroid(length, mass, fluidDensity):
    # Equivalent Diameter
    D_s = np.sqrt(np.divide(6*mass,fluidDensity*np.pi*length))
    # volume of spheroid    
    return (4*np.pi/3) * np.power(D_s/2, 2) * length/2

def ellipsoidApproximation(length, width, thickness):
    a = length/2
    b = width/2
    c = thickness/2
                              
    return 4*np.pi/3*np.power(b, 2) * a

def partitionDisc(length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop):
    dx = np.linspace(0.0, 1.0, 100)
    topPolySide = np.poly1d(topCoeffSide)
    bottomPolySide = np.poly1d(bottomCoeffSide) 
    topPolyTop = np.poly1d(topCoeffTop)
    bottomPolyTop = np.poly1d(bottomCoeffTop)
    # scale top and bottom values by length
    # surface area sums up the N-2 frustrums
    topPointsSide = np.abs(topPolySide(dx))*length
    bottomPointsSide = np.abs(bottomPolySide(dx))*length
    topPointsTop = np.abs(topPolyTop(dx))*length
    bottomPointsTop = np.abs(bottomPolyTop(dx))*length
    # distance between each circumference
    t = (dx[1]-dx[0])*length
    width = topPointsTop + bottomPointsTop
    height = topPointsSide + bottomPointsSide
    
    b = height/2
    a = width/2
    
    estVolume = np.zeros(len(dx))
    
    for i in range(len(estVolume)-1):
        estVolume[i] = np.pi * t * a[i] * b[i]
        
    return np.sum(estVolume)