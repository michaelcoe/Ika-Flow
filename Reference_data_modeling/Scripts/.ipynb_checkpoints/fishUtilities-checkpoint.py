import math
import cv2 as cv2
import numpy as np
from scipy import integrate
from scipy.special import ellipe
from pathlib import Path
#----- Method List -----------------------
# list of images = get_image_files(imagePath, fileName)
# image = kmeans_color_quantization(image, clusters, rounds)
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
    for path in Path(imagePath).rglob('**/'+fileName):
        imageFiles.append(path)
                
    return imageFiles

def kmeans_color_quantization(image, clusters=2, rounds=1):
    Z = image.reshape((-1, 3))
    Z = np.float32(Z)

    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 10, 1.0)
    _, label, center = cv2.kmeans(Z, clusters, None, criteria, rounds, 
    cv2.KMEANS_RANDOM_CENTERS)

    center = np.uint8(center)
    result = center[label.flatten()]

    return result.reshape((image.shape))

def get_contours(image):
    # convert to hsv
    hsvImg = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    # split the channels
    h, s, v = cv2.split(hsvImg)
    # Gaussian filter the s and v component
    s = cv2.GaussianBlur(s, (3,3), 1)
    v = cv2.GaussianBlur(v, (3,3), 1)
    # append the images
    hsvImg_processed = cv2.merge([h, s, v])
    # color quantize the image
    quantized = kmeans_color_quantization(hsvImg_processed)
    # morphological close and open
    kernal = np.ones((21, 21), np.uint8)
    closeImg = cv2.morphologyEx(quantized, cv2.MORPH_CLOSE, kernal)
    openImg = cv2.morphologyEx(closeImg, cv2.MORPH_OPEN, kernal)

    # Canny Edge Detection
    edged = cv2.Canny(openImg,100, 200)
    # find the contours
    contours, _ = cv2.findContours(edged, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    # sort the contours
    cnt = max(contours, key=cv2.contourArea)
    # get the convex hull for the contours
    hull = cv2.convexHull(cnt, False)
    
    return cnt, hull

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
    extLeft = tuple(cnts[cnts[:,:,0].argmin()][0])
    extRight = tuple(cnts[cnts[:,:,0].argmax()][0])
    midPoint_left = (int(extLeft[0] + (cx-extLeft[0])/2), cy)
    midPoint_right = (int(cx + (extRight[0]-cx)/2), cy)
    
    # find the dividing line
    dlFit = np.polyfit([extLeft[0], midPoint_left[0], cx, midPoint_right[0], extRight[0]], 
                       [extLeft[1], midPoint_left[1], cy, midPoint_right[1], extRight[1]], 4)
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