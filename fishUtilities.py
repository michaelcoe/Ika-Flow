import os
import math
import cv2 as cv2
import numpy as np
from scipy import integrate
from scipy.special import ellipe

#----- Method List -----------------------
# list of images = get_image_files(imagePath, fileName)
# contours, hull = get_contours(image, contourNumber)
# (cx, cy), areaTotal, topHull, bottomHull = split_by_centroid_side(imageShape, cnts)
# (cx, cy), tHullUnique, bHullUnique = split_by_centroid_top(imageShape, cnts)
# (cx, cy), aspectRatio, topHull, bottomHull = split_by_centroid_cs(imageShape, cnts)
# scale factor, tophull, bottomhull = scale_data(topHull, bottomHull)
# total area, tophull, bottomhull = scale_data_top(topHull, bottomHull)
# tophull, bottomhull = scale_data_cs(topHull, bottomHull)
# [max location, max], [min location, min] = get_min_max(topHull, bottomHull)
# ellipseTop, ellipseBottom = fit_ellipse(maxPoints, minPoints)
#----------------------------------------

def get_image_files(imagePath, fileName):
    imageFiles = []

    for root, dirs, files in os.walk(imagePath):
        for file in files:
            if file.endswith(fileName):
                imageFiles.append(os.path.join(root, file))
                
    return imageFiles

def get_contours(image, contourNumber):
    # change the image to a gray color
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    #kernal
    kernel = np.ones((9,9), np.uint8)
    # blur the image
    blurImg = cv2.medianBlur(gray, 5)
    # threshold and renormalize
    thrImg = cv2.adaptiveThreshold(blurImg,255,cv2.ADAPTIVE_THRESH_GAUSSIAN_C,\
                    cv2.THRESH_BINARY_INV,7,2)
    # find the contours
    contours, hierarchy = cv2.findContours(thrImg, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    # sort the contours
    cnts = sorted(contours, key=cv2.contourArea, reverse=True)
    cnts = cnts[contourNumber]
    # get the convex hull for the contours
    hull = cv2.convexHull(cnts, False)
    
    return cnts, hull

def split_by_centroid_side(imageShape, cnts):
    # create a blank image
    mask = np.zeros(imageShape, dtype=np.uint8)
    # find the pixel points
    cv2.drawContours(mask, [cnts], -1, 255, 1)
    pixelPoints = np.transpose(np.nonzero(mask))
    hullX = pixelPoints[:,1]
    hullY = pixelPoints[:,0]
    
    # sort the hull components
    arrind = hullX.argsort()
    hullX = hullX[arrind]
    hullY = hullY[arrind]
    # find centroid of contour
    M = cv2.moments(cnts)
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    # find the dividing lines
    leftMost = tuple(cnts[cnts[:,:,0].argmin()][0])
    rightMost = tuple(cnts[cnts[:,:,0].argmax()][0])
    #topMost = tuple(cnts[cnts[:,:,1].argmin()][0])
    #bottomMost = tuple(cnts[cnts[:,:,1].argmax()][0])
    
    # find the dividing line
    rightAngle = np.arctan2(cy-rightMost[1], rightMost[0]-cx)*180/np.pi
    leftAngle = np.arctan2(cy-leftMost[1], cx-leftMost[0])*180/np.pi
    if(np.abs(rightAngle) > 5 or np.abs(leftAngle) > 5):
        dlFit = np.polyfit([leftMost[0], cx, rightMost[0]], [cy, cy, cy], 2)
        divideLine = np.poly1d(dlFit)
    else:
        dlFit = np.polyfit([leftMost[0], cx, rightMost[0]], [leftMost[1], cy, rightMost[1]], 2)
        divideLine = np.poly1d(dlFit)
    
    # divide contours into a top and bottom contour
    topHullDLX = []
    topHullDLY = []
    bottomHullDLX = []
    bottomHullDLY = []
    
    for idx, x in enumerate(hullX):
        Yoffset = cy-divideLine(hullX[idx])
        if(hullY[idx] < divideLine(hullX[idx])):
            topHullDLX.append(hullX[idx])
            topHullDLY.append(hullY[idx]+Yoffset)
        else:
            bottomHullDLX.append(hullX[idx])
            bottomHullDLY.append(hullY[idx]+Yoffset)
    
    topHull = np.array([topHullDLX, topHullDLY])
    bottomHull = np.array([bottomHullDLX, bottomHullDLY])
    
    # shift x over to 0
    bottomHull[0,:] = bottomHull[0,:] - bottomHull[0,0]
    topHull[0,:] = topHull[0,:] - topHull[0,0]
    # shift y to 0
    bottomHull[1,:] = bottomHull[1,:] - cy
    topHull[1,:] = topHull[1,:] - cy
    # flip over the x axis
    bottomHull[1,:] = -bottomHull[1,:]
    topHull[1,:] = -topHull[1,:]
    
    _, tUniqueIndex = np.unique(topHull[0,:], return_index=True)
    _, bUniqueIndex = np.unique(bottomHull[0,:], return_index=True)
    
    tHullUnique = np.array([topHull[0, tUniqueIndex], topHull[1, tUniqueIndex]])
    bHullUnique = np.array([bottomHull[0, bUniqueIndex], bottomHull[1, bUniqueIndex]])
    
    areaTop = integrate.simps(tHullUnique[1,:], tHullUnique[0,:])
    areaBottom = integrate.simps(bHullUnique[1,:], bHullUnique[0,:])
    
    areaTotal = areaTop + np.abs(areaBottom)
    
    return (cx, cy), areaTotal, topHull, bottomHull

def split_by_centroid_top(imageShape, cnts):
    # create a blank image
    mask = np.zeros(imageShape, dtype=np.uint8)
    # find the pixel points
    cv2.drawContours(mask, [cnts], -1, 255, 1)
    pixelPoints = np.transpose(np.nonzero(mask))
    hullX = pixelPoints[:,1]
    hullY = pixelPoints[:,0]
    
    # sort the hull components
    arrind = hullX.argsort()
    hullX = hullX[arrind]
    hullY = hullY[arrind]
    
    # find centroid of contour
    M = cv2.moments(cnts)
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    # find the dividing lines
    leftMost = tuple(cnts[cnts[:,:,0].argmin()][0])
    rightMost = tuple(cnts[cnts[:,:,0].argmax()][0])
    #topMost = tuple(cnts[cnts[:,:,1].argmin()][0])
    #bottomMost = tuple(cnts[cnts[:,:,1].argmax()][0])
    
    # find the dividing line
    rightAngle = np.arctan2(cy-rightMost[1], rightMost[0]-cx)*180/np.pi
    leftAngle = np.arctan2(cy-leftMost[1], cx-leftMost[1])*180/np.pi
    if(np.abs(rightAngle) > 10 or np.abs(leftAngle) > 10):
        dlFit = np.polyfit([leftMost[0], cx, rightMost[0]], [cy, cy, cy], 2)
        divideLine = np.poly1d(dlFit)
    else:
        dlFit = np.polyfit([leftMost[0], cx, rightMost[0]], [leftMost[1], cy, rightMost[1]], 2)
        divideLine = np.poly1d(dlFit)
    
    # divide contours into a top and bottom contour
    topHullDLX = []
    topHullDLY = []
    bottomHullDLX = []
    bottomHullDLY = []
    
    for idx, x in enumerate(hullX):
        Yoffset = cy-divideLine(hullX[idx])
        if(hullY[idx] < divideLine(hullX[idx])):
            topHullDLX.append(hullX[idx])
            topHullDLY.append(hullY[idx]+Yoffset)
        else:
            bottomHullDLX.append(hullX[idx])
            bottomHullDLY.append(hullY[idx]+Yoffset)
    
    topHull = np.array([topHullDLX, topHullDLY])
    bottomHull = np.array([bottomHullDLX, bottomHullDLY])
    
    # shift x over to 0
    bottomHull[0,:] = bottomHull[0,:] - bottomHull[0,0]
    topHull[0,:] = topHull[0,:] - topHull[0,0]
    # shift y to 0
    bottomHull[1,:] = bottomHull[1,:] - cy
    topHull[1,:] = topHull[1,:] - cy
    # flip over the x axis
    bottomHull[1,:] = -bottomHull[1,:]
    topHull[1,:] = -topHull[1,:]
    
    _, tUniqueIndex = np.unique(topHull[0,:], return_index=True)
    _, bUniqueIndex = np.unique(bottomHull[0,:], return_index=True)
    
    tHullUnique = np.array([topHull[0, tUniqueIndex], topHull[1, tUniqueIndex]])
    bHullUnique = np.array([bottomHull[0, bUniqueIndex], bottomHull[1, bUniqueIndex]])
        
    return (cx, cy), tHullUnique, bHullUnique

def split_by_centroid_cs(imageShape, cnts):
    # create a blank image
    mask = np.zeros(imageShape, dtype=np.uint8)
    # find the pixel points
    cv2.drawContours(mask, [cnts], -1, 255, 1)
    pixelPoints = np.transpose(np.nonzero(mask))
    hullX = pixelPoints[:,1]
    hullY = pixelPoints[:,0]
    
    # sort the hull components
    arrind = hullX.argsort()
    hullX = hullX[arrind]
    hullY = hullY[arrind]
    
    # find centroid of contour
    M = cv2.moments(cnts)
    cx = int(M['m10']/M['m00'])
    cy = int(M['m01']/M['m00'])
    
    # find the dividing lines
    leftMost = tuple(cnts[cnts[:,:,0].argmin()][0])
    rightMost = tuple(cnts[cnts[:,:,0].argmax()][0])
    topMost = tuple(cnts[cnts[:,:,1].argmin()][0])
    bottomMost = tuple(cnts[cnts[:,:,1].argmax()][0])
    
    aspectRatio = np.abs((rightMost[0]-leftMost[0])/(topMost[1]-bottomMost[1]))
    dlFit = np.polyfit([leftMost[0], cx, rightMost[0]], [cy, cy, cy], 1)
    divideLine = np.poly1d(dlFit)
    
    # divide contours into a top and bottom contour
    topHullDLX = []
    topHullDLY = []
    bottomHullDLX = []
    bottomHullDLY = []
    
    for idx, x in enumerate(hullX):
        Yoffset = cy-divideLine(hullX[idx])
        if(hullY[idx] < divideLine(hullX[idx])):
            topHullDLX.append(hullX[idx])
            topHullDLY.append(hullY[idx]+Yoffset)
        else:
            bottomHullDLX.append(hullX[idx])
            bottomHullDLY.append(hullY[idx]+Yoffset)
    
    topHull = np.array([topHullDLX, topHullDLY])
    bottomHull = np.array([bottomHullDLX, bottomHullDLY])
    
    # shift x over to 0
    bottomHull[0,:] = bottomHull[0,:] - bottomHull[0,0]
    topHull[0,:] = topHull[0,:] - topHull[0,0]
    # shift y to 0
    bottomHull[1,:] = bottomHull[1,:] - cy
    topHull[1,:] = topHull[1,:] - cy
    # flip over the x axis
    bottomHull[1,:] = -bottomHull[1,:]
    topHull[1,:] = -topHull[1,:]
    
    return (cx, cy), aspectRatio, topHull, bottomHull

def scale_data(topHull, bottomHull):
    # scale the hull
    if(topHull[0,-1] > bottomHull[0,-1]):
        scaleFactor = topHull[0,-1]
        tHull = np.array([[np.divide(topHull[0,:], topHull[0,-1])], [np.divide(topHull[1,:], topHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], topHull[0,-1])], [np.divide(bottomHull[1,:], topHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
    else:
        scaleFactor = bottomHull[0,-1]
        tHull = np.array([[np.divide(topHull[0,:], bottomHull[0,-1])], [np.divide(topHull[1,:], bottomHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], bottomHull[0,-1])], [np.divide(bottomHull[1,:], bottomHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
       
    return scaleFactor, tHull, bHull

def scale_data_top(topHull, bottomHull):
    # scale the hull
    if(topHull[0,-1] > bottomHull[0,-1]):
        tHull = np.array([[np.divide(topHull[0,:], topHull[0,-1])], [np.divide(topHull[1,:], topHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], topHull[0,-1])], [np.divide(bottomHull[1,:], topHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
    else:
        tHull = np.array([[np.divide(topHull[0,:], bottomHull[0,-1])], [np.divide(topHull[1,:], bottomHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], bottomHull[0,-1])], [np.divide(bottomHull[1,:], bottomHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
    
    areaTop = integrate.simps(tHull[1,:], tHull[0,:])
    areaBottom = integrate.simps(bHull[1,:], bHull[0,:])
    
    areaTotal = areaTop + np.abs(areaBottom)
    
    return areaTotal, tHull, bHull

def scale_data_cs(topHull, bottomHull):
    # scale the hull
    if(topHull[0,-1] > bottomHull[0,-1]):
        tHull = np.array([[np.divide(topHull[0,:], topHull[0,-1])], [np.divide(topHull[1,:], topHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], topHull[0,-1])], [np.divide(bottomHull[1,:], topHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
    else:
        tHull = np.array([[np.divide(topHull[0,:], bottomHull[0,-1])], [np.divide(topHull[1,:], bottomHull[0,-1])]])
        bHull = np.array([[np.divide(bottomHull[0,:], bottomHull[0,-1])], [np.divide(bottomHull[1,:], bottomHull[0,-1])]])
        tHull = tHull.reshape(2,-1)
        bHull = bHull.reshape(2,-1)
    
    return tHull, bHull

def get_min_max(topHull, bottomHull):
    minIndex = np.argmin(bottomHull[1,:])
    maxIndex = np.argmax(topHull[1,:])
    return (np.array([topHull[0, maxIndex], topHull[1, maxIndex]]),
    np.array([bottomHull[0, minIndex], bottomHull[1, minIndex]]))

def fit_ellipse(maxPoints, minPoints):
    aEllipse = maxPoints[0]
    bEllipseTop = maxPoints[1]
    bEllipseBottom = np.abs(minPoints[1])
    tEllipseTop = np.linspace(0, np.pi, 100)
    tEllipseBottom = np.linspace(0,-np.pi, 100)
    ellipseTop = np.array([0.5+aEllipse*np.cos(tEllipseTop), bEllipseTop*np.sin(tEllipseTop)])
    ellipseBottom = np.array([0.5+aEllipse*np.cos(tEllipseBottom), bEllipseBottom*np.sin(tEllipseBottom)])
        
    return ellipseTop, ellipseBottom