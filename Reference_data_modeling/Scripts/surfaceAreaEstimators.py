import numpy as np
import Scripts.airfoil as af
from scipy import integrate
from scipy.special import ellipe

#-------------Method List----------------------------
# determine_surface_area(topForm, crossSection, aspectRatio, length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop)
# surface_area_angulliform(xu, topPolySide, bottomPolySide, heightWidthRatio, length)
# surface_area_fusiform(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_skate(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_teardrop(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_oval(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_box(topPolySide, bottomPolySide, xu, yu, length)
# surface_area_triangle(topPolySide, bottomPolySide, xu, yu, length)
# equivalentSpheroid(length, mass, fluidDensity)
#------------------------------------------------------

def determine_surface_area(topForm, crossSection, aspectRatio, length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop):
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
        return surface_area_angulliform(xu, topPolySide, bottomPolySide, aspectRatio, length)
    elif(crossSection==2):
        return surface_area_fusiform(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==3):
        return surface_area_skate(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==4):
        return surface_area_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==5):
        return surface_area_teardrop(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==6):
        return surface_area_oval(topPolySide, bottomPolySide, xu, yu, length)
    elif(crossSection==7):
        return surface_area_box(topPolySide, bottomPolySide, xu, yu, length)
    else:
        return surface_area_triangle(topPolySide, bottomPolySide, xu, yu, length)

def surface_area_angulliform(xu, topPolySide, bottomPolySide, heightWidthRatio, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    topPointsTop = np.divide(topPointsSide+bottomPointsSide, heightWidthRatio*2)
    # distance between each circumference
    t = (xu[1]-xu[0])
    # calculate the eccentricity of the top and bottom
    eSquaredTop = 1.0 - np.divide(np.power(topPointsSide, 2), np.power(topPointsTop, 2))
    eSquaredBottom = 1.0 - np.divide(np.power(bottomPointsSide, 2), np.power(topPointsTop, 2))
    # calculate the surface area based on circumference and thickness
    topArea = np.sum(np.multiply(np.multiply(topPointsTop, ellipe(eSquaredTop)),t))
    bottomArea = np.sum(np.multiply(np.multiply(topPointsTop, ellipe(eSquaredBottom)),t))
    # combine to get surface area of one side and multiple by 2 for total surface area    
    return 2.0*(topArea+bottomArea)

def surface_area_fusiform(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    # surface area sums up the N-2 frustrums
    topPointsSide = np.abs(topPolySide(xu[1:]))*length
    bottomPointsSide = np.abs(bottomPolySide(xu[1:]))*length
    xu = xu[1:]*length
    yu = yu[1:]*length
    # distance between each circumference
    t = (xu[1]-xu[0])
    height = np.abs(topPointsSide) + np.abs(bottomPointsSide)
    width = np.abs(2*yu) 
    # calculate the eccentricity of the top and bottom
    #eSquaredTop = 1.0 - np.divide(np.power(yu,2), np.power(topPointsSide,2))
    #eSquaredBottom = 1.0 - np.divide(np.power(yu,2), np.power(bottomPointsSide,2))
    eSquaredTop = 1.0 - np.divide(np.power(topPointsSide,2),np.power(yu,2))
    eSquaredBottom = 1.0 - np.divide(np.power(bottomPointsSide,2), np.power(yu,2))
    # calculate 1/4 of the circumference representing one side of the fish
    #topArea = np.sum(np.multiply(topPointsSide, ellipe(eSquaredTop)))*t
    #bottomArea = np.sum(np.multiply(bottomPointsSide, ellipe(eSquaredBottom)))*t
    topArea = np.sum(np.multiply(yu, ellipe(eSquaredTop)))*t
    bottomArea = np.sum(np.multiply(yu, ellipe(eSquaredBottom)))*t
    # combine to get surface area of one side and multiple by 2 for total surface area    
    return 2.0*(topArea+bottomArea), np.max([height, width])

def surface_area_skate(topPolySide, bottomPolySide, xu, yu, length):
    # scale top and bottom values by length
    topPointsSide = np.abs(topPolySide(xu))*length
    bottomPointsSide = np.abs(bottomPolySide(xu))*length
    xu = xu*length
    yu = yu*length
    # Thickness of the frustrums
    t = (xu[1]-xu[0])
    # find the area of the Side contour which is the top of the fish
    topArea = integrate.simps(topPointsSide, xu)
    bottomArea = integrate.simps(bottomPointsSide, xu)
    # calculate the eccentricity of the top and bottom
    # note that skates are the inverse of the fusilform body
    eSquaredTop = 1.0 - np.divide(np.power(topPointsSide,2), np.power(yu,2),)
    eSquaredBottom = 1.0 - np.divide(np.power(bottomPointsSide,2), np.power(yu,2))
    # calculate 1/4 of the circumference representing one side of the fish
    topArea = np.sum(np.multiply(np.multiply(yu, ellipe(eSquaredTop)),t))
    bottomArea = np.sum(np.multiply(np.multiply(yu, ellipe(eSquaredBottom)),t))
        
    return 2.0*(topArea+bottomArea)

def surface_area_invertedTeardrop(topPolySide, bottomPolySide, xu, yu, length):
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
    eSquaredTop = 1.0 - np.divide(np.power(yu[1:firstIndex], 2), np.power(topPointsSide[1:firstIndex], 2))
    topEllipseArea = np.sum(np.multiply(np.multiply(topPointsSide[1:firstIndex], ellipe(eSquaredTop)),t))
    # cross section is flattened at the bottom, assume that around 5 percent of the triangle side will fold over to make flattened bottom
    bottomTriangleArea = np.sum(np.multiply(np.power(bottomPointsSide[0:firstIndex]+(bottomPointsSide[0:firstIndex]*0.05),2)+
                                            np.power(yu[0:firstIndex],2),t))
    #-------------------------------------------------------
    # calculate the area of the ellipse portion of the fish (last 1/3)
    # calculate the eccentricity of the top and bottom
    eSquaredTop = 1.0 - np.divide(np.power(yu[lastIndex:-1], 2), np.power(topPointsSide[lastIndex:-1], 2))
    eSquaredBottom = 1.0 - np.divide(np.power(yu[lastIndex:-1], 2), np.power(bottomPointsSide[lastIndex:-1], 2))
    # calculate 1/4 of the circumference representing one side of the fish
    topTailArea = np.sum(np.multiply(np.multiply(topPointsSide[lastIndex:-1], ellipe(eSquaredTop)),t))
    bottomTailArea= np.sum(np.multiply(np.multiply(bottomPointsSide[lastIndex:-1], ellipe(eSquaredBottom)),t))
    # calculate the surface area based on circumference and thickness
         
    return 2.0*(topEllipseArea+bottomTriangleArea+topTailArea+bottomTailArea)

def surface_area_teardrop(topPolySide, bottomPolySide, xu, yu, length):
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
    topTriangleArea = np.sum(np.multiply(np.power(yu[0:firstIndex]+(topPointsSide[0:firstIndex]*0.05),2)+np.power(yu[0:firstIndex],2),t))
    # bottom Ellipse
    eSquaredBottom = 1.0 - np.divide(np.power(yu[0:firstIndex], 2), np.power(bottomPointsSide[0:firstIndex], 2))
    bottomEllipseArea = np.sum(np.multiply(np.multiply(bottomPointsSide[0:firstIndex], ellipe(eSquaredBottom)),t))
    #-------------------------------------------------------
    # calculate the area of the ellipse portion of the fish (last 1/3)
    # calculate the eccentricity of the top and bottom
    eSquaredTop = 1.0 - np.divide(np.power(yu[lastIndex:-1], 2), np.power(topPointsSide[lastIndex:-1], 2))
    eSquaredBottom = 1.0 - np.divide(np.power(yu[lastIndex:-1], 2), np.power(bottomPointsSide[lastIndex:-1], 2))
    # calculate 1/4 of the circumference representing one side of the fish
    topTailArea = np.sum(np.multiply(np.multiply(topPointsSide[lastIndex:-1], ellipe(eSquaredTop)),t))
    bottomTailArea= np.sum(np.multiply(np.multiply(bottomPointsSide[lastIndex:-1], ellipe(eSquaredBottom)),t))
    # calculate the surface area based on circumference and thickness
         
    return 2.0*(bottomEllipseArea+topTriangleArea+topTailArea+bottomTailArea)

def surface_area_oval(dx, topPolySide, bottomPolySide, topPolyTop, bottomPolyTop, length):
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

def surface_area_box(topPolySide, bottomPolySide, xu, yu, length):
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
    topBoxArea = np.sum(np.multiply(topPointsTop[0:firstIndex]+topPointsSide[0:firstIndex],t))
    bottomBoxArea = np.sum(np.multiply(topPointsTop[0:firstIndex]+bottomPointsSide[0:firstIndex],t))
    #-------------------------------------------------------
    # calculate the area of the ellipse portion of the fish (last 1/3)
    # calculate the eccentricity of the top and bottom
    eSquaredTop = 1.0 - np.divide(np.power(topPointsTop[lastIndex:-1], 2), np.power(topPointsSide[lastIndex:-1], 2))
    eSquaredBottom = 1.0 - np.divide(np.power(topPointsTop[lastIndex:-1], 2), np.power(bottomPointsSide[lastIndex:-1], 2))
    # calculate 1/4 of the circumference representing one side of the fish
    # calculate the surface area based on circumference and thickness
    topTailArea = np.sum(np.multiply(np.multiply(topPointsSide[lastIndex:-1], ellipe(eSquaredTop)),t))
    bottomTailArea= np.sum(np.multiply(np.multiply(bottomPointsSide[lastIndex:-1], ellipe(eSquaredBottom)),t))
         
    return 2.0*(topBoxArea+bottomBoxArea+topTailArea+bottomTailArea)

def surface_area_triangle(topPolySide, bottomPolySide, xu, yu, length):
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
    # ellipsicity
    ellipsicity = np.sqrt(1-np.divide(np.square(D_s), np.square(length)))
    # Equivalent surface area
    A_s = 2*np.pi*np.power(D_s,2)/4 + 2*np.pi*np.multiply(np.divide(np.multiply(D_s, length),4*ellipsicity), np.arcsin(ellipsicity))
    # Slenderness Ratio
    S_r = np.divide(length, D_s)
    # Geometric constant for length
    # slenderness ratio
    slenderness = np.divide(length, D_s)
    c = np.power(np.divide(np.multiply(6, np.power(S_r,2)), np.multiply(fluidDensity, np.pi)), 1/3)
    calcLength = c*np.power(mass, 1/3)
    newSr = np.divide(calcLength, D_s)
    # Geometric constant for surface area
    d = np.power(np.divide(1,fluidDensity),2/3)*(-0.0122*np.power(slenderness,2)+0.5196*slenderness+4.2732)
    A_sr = np.multiply(d, np.power(mass, 2/3))
    
    return D_s, A_s, slenderness

def ellipsoidApproximation(length, width, thickness):
    a = length/2
    b = width/2
    c = thickness/2

    eccentricity = np.sqrt(1.0-(b**2/a**2))

    return 2*np.pi*b**2 + 2*np.pi*np.divide(a*b,eccentricity)*np.arcsin(eccentricity)

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
    
    estArea = np.zeros(len(dx))
    
    for i in range(len(estArea)-1):
        estArea[i] = (((b[i]+a[i])+(b[i+1]+a[i+1]))/2)*t*np.pi
        
    return np.sum(estArea)