import os
import cv2 as cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import Scripts.surfaceAreaEstimators as sae
import Scripts.fishUtilities as fu
import Scripts.fishFits as ff
from scipy import stats
from scipy.optimize import curve_fit
from scipy.stats.distributions import  t

def get_data(imagePath, specimen):
    imagePathTop = fu.get_image_files(imagePath, specimen+'Top.jpg')
    imagePathSide = fu.get_image_files(imagePath, specimen+'SideNF.jpg')
    imagePathCS = fu.get_image_files(imagePath, specimen+'CS.jpg')
    imagePathSideFull = fu.get_image_files(imagePath, specimen+'Side.jpg')
    
    # import images
    imageTop = cv2.imread(imagePathTop[0])
    imageSide = cv2.imread(imagePathSide[0])
    imageCS = cv2.imread(imagePathCS[0])
    imageSideFull = cv2.imread(imagePathSideFull[0])
    # get contours
    cntsTop, hullTop = fu.get_contours(imageTop, 1)
    cntsSide, hullSide = fu.get_contours(imageSide, 1)
    cntsCS, hullCS = fu.get_contours(imageCS, 1)
    cntsFull, hullFUll = fu.get_contours(imageSideFull, 1)
    # Find the area ratio of fins and no fins
    areaFull = cv2.contourArea(cntsFull)
    areaNF = cv2.contourArea(cntsSide)
    areaRatio = (1-np.abs(areaNF/areaFull))
    # split the hull into top and bottom
    (cxTop, cyTop), topHullTop, bottomHullTop = fu.split_by_centroid_top(imageTop.shape, cntsTop)
    (cxSide, cySide), Area, topHullSide, bottomHullSide = fu.split_by_centroid_side(imageSide.shape, cntsSide)
    (cxCS, cyCS), aspectRatio, topHullCS, bottomHullCS = fu.split_by_centroid_cs(imageCS.shape, cntsCS)
    # get contour data
    areaTop, tHullTop, bHullTop = fu.scale_data_top(topHullTop, bottomHullTop)
    areaSide, tHullSide, bHullSide = fu.scale_data(topHullSide, bottomHullSide)
    tHullCS, bHullCS = fu.scale_data_cs(topHullCS, bottomHullCS)
    # get min and max points
    maxPointsSide, minPointsSide = fu.get_min_max(tHullSide, bHullSide)
    maxPointsTop, minPointsTop = fu.get_min_max(tHullTop, bHullTop)
    maxPointsCS, minPointsCS = fu.get_min_max(tHullCS, bHullCS)
        
    width = np.array([maxPointsTop[1]+np.abs(minPointsTop[1])])
    height = np.array([maxPointsSide[1]+np.abs(minPointsSide[1])])
    # fit top data
    nacaArea, m, t, d0, xu, yu, xl, yl = ff.fit_top_contours(areaTop, np.mean([minPointsTop[0],maxPointsTop[0]]),
                                              width, (tHullTop[1,-5]+np.abs(bHullTop[1,-5]))/2.0)
    nacaData = np.array([m, t, d0])
    polyArea, topCoeffTop, bottomCoeffTop = ff.fit_top_contours_poly(tHullTop, bHullTop, 6)
    # fit side data
    topCoeffSide, bottomCoeffSide = ff.fit_side_contours(tHullSide, bHullSide, 6)
       
    return (width, height, areaRatio, topCoeffSide, bottomCoeffSide,
            topCoeffTop, bottomCoeffTop, nacaData)

def get_surface_area(species, imagePath, testLength, testMass, scale):
    # fluid parameters
    density = 1025
    kinViscosity = 0.00000089266
    
    (width, height, areaRatio, topCoeffSide, bottomCoeffSide,
     topCoeffTop, bottomCoeffTop, nacaData) = get_data(imagePath, species)
    # pre-allocate arrays
    method = np.zeros(testLength.shape)
    ellipsoid = np.zeros(testLength.shape)
    partition = np.zeros(testLength.shape)
    esAs = np.zeros(testLength.shape)
    index = np.arange(0, len(testLength),1)
    
    for idx, length, mass in zip(index, testLength, testMass):
            _, _, _, esAs[idx], _ = sae.equivalentSpheroid(np.divide(length,100), mass/1000, density)
            method[idx] = sae.determine_surface_area(True, 2, 0.15, length, topCoeffSide, bottomCoeffSide,
                                                        nacaData, nacaData)
            ellipsoid[idx] = sae.ellipsoidApproximation(length, width*length, height*length)
            partition[idx] = sae.partitionDisc(length, topCoeffSide, bottomCoeffSide, topCoeffTop, bottomCoeffTop)
    
    print(method)
    method = method + (method*areaRatio)
    print(method)
    #ellipsoid = ellipsoid + (ellipsoid*areaRatio)
    #partition = partition + (partition*areaRatio)
    #esAs = esAs + (esAs*areaRatio)
    
    esAs = np.multiply(esAs, np.power(100,2))
    
    return method*scale, partition*scale, ellipsoid*scale, esAs

def power_law(x, a, b):
    return a*np.power(x,b)

def fit_data(x, y, p0, dx):
    fitParams, fitCov = curve_fit(power_law, x, y,np.asarray(p0))
    
    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)
    
    n = len(y)    # number of data points
    p = len(fitParams) # number of parameters

    dof = max(0, n - p) # number of degrees of freedom

    # student-t value for the dof and confidence level
    tval = t.ppf((1.0-alpha)/2., dof) 
    
    lower = []
    upper = []
    
    for p, var in zip(fitParams, np.diag(fitCov)):
        sigma = var**0.5
        lower.append(p-sigma*tval)
        upper.append(p+sigma*tval)
    
    return power_law(dx, *fitParams), power_law(dx, *upper), power_law(dx, *lower)

def main():
    homeImagePath = r'./Pictures'
    homeFigurePath = r'./'
    homeDatabasePath = r'C:\Users\micha\Dropbox\UUV Project\Databases'
    # UniFilePaths
    uniFigurePath = r'D:\Dropbox\Journal Articles\Surface area and COT\figures'
    uniDatabasePath = r'D:\Dropbox\UUV Project\Databases'
    uniImagePath = r'D:\Dropbox\UUV Project\Pictures'
    
    database = r'osheaData.xlsx'
    dataPath = os.path.join(uniDatabasePath, database)
    
    imagePath = uniImagePath
    figPath = uniFigurePath
    # get path file names
    species = ['atlanticCod', 'scannedSalmon']
        
    # Published Data format Length [cm], Surface Area [cm^2], Mass [g], Surface Area [cm^2]
    dfSalmon = pd.read_excel(dataPath, 'Salmon')
    dfCod = pd.read_excel(dataPath, 'Cod')
    
    OsheaSalmonLength = dfSalmon['Length [cm]'].values
    salmonSA = dfSalmon['SA [cm^2]'].values
    OsheaCodLength = dfCod['Length [cm]'].values
    codSA = dfCod['SA [cm^2]'].values
        
    # Scanned Data
    scannedAs = np.array([np.multiply(0.0674504559,np.power(100,2))])
    scannedLength = np.array(35.599)
    scannedData = np.array([scannedLength, scannedAs])
    salmonLength = np.append(OsheaSalmonLength, scannedLength)
    salmonSA = np.append(salmonSA, scannedAs)
    codLength = OsheaCodLength
    
    # Length-weight relationships (LWR) are given at a*L^b and gotten from fishbase
    LWR = np.array([[0.0069, 3.08],[0.0120, 3.0]])
    salmonLength = np.linspace(salmonLength[0], salmonLength[-1], 20)
    codLength = np.linspace(OsheaCodLength[0], OsheaCodLength[-1]+0.1, 20)
    
    codMass = LWR[0,0]*np.power(codLength,LWR[0,1])
    salmonMass = LWR[1,0]*np.power(salmonLength, LWR[1,1])
    
    methodCod, partitionCod, ellipsoidCod, esAsCod = get_surface_area(species[0], imagePath, codLength, codMass, 1)
    methodSalmon, partitionSalmon, ellipsoidSalmon, esAsSalmon = get_surface_area(species[1], imagePath, salmonLength, salmonMass, 1)
    
    # Atlantic Cod and Atlantic Salmon regression lines
    salmonLRes, salmonLUncUpper, salmonLUncLower = fit_data(OsheaSalmonLength, salmonSA[:-1], [0.72, 1.88], OsheaSalmonLength)
    codLRes, codLUncUpper, codLUncLower = fit_data(OsheaCodLength, codSA, [0.78, 1.89], OsheaCodLength)

    # Regression lines for generated data
    OsheaParams, OsheaConv = curve_fit(power_law, OsheaSalmonLength, salmonSA[:-1],np.asarray([0.72, 1.88]))
    methodParams, methodConv = curve_fit(power_law, salmonLength, methodSalmon,np.asarray([0.72, 1.88]))
    partitionParams, partitionConv = curve_fit(power_law, salmonLength, partitionSalmon,np.asarray([0.72, 1.88]))
    ellipsoidParams, ellipsoidConv = curve_fit(power_law, salmonLength, ellipsoidSalmon,np.asarray([0.72, 1.88]))
    esParams, esConv = curve_fit(power_law, salmonLength, esAsSalmon,np.asarray([0.72, 1.88]))
    salmonParams = np.array([OsheaParams, methodParams, partitionParams, ellipsoidParams, esParams])
    
    # Regression lines for generated data
    OsheaParams, OsheaConv = curve_fit(power_law, OsheaCodLength, codSA,np.asarray([0.72, 1.88]))
    methodParams, methodConv = curve_fit(power_law, codLength, methodCod,np.asarray([0.72, 1.88]))
    partitionParams, partitionConv = curve_fit(power_law, codLength, partitionCod,np.asarray([0.72, 1.88]))
    ellipsoidParams, ellipsoidConv = curve_fit(power_law, codLength, ellipsoidCod,np.asarray([0.72, 1.88]))
    esParams, esConv = curve_fit(power_law, codLength, esAsCod,np.asarray([0.72, 1.88]))
    codParams = np.array([OsheaParams, methodParams, partitionParams, ellipsoidParams, esParams])
    '''
    # calculate the Error
    methodErrorSalmon = np.abs(np.divide(methodSalmon[:-1]-salmonSA[:-1], salmonSA[:-1]))*100
    ellipsoidErrorSalmon = np.abs(np.divide(ellipsoidSalmon[:-1]-salmonSA[:-1], salmonSA[:-1]))*100
    partitionErrorSalmon = np.abs(np.divide(partitionSalmon[:-1]-salmonSA[:-1], salmonSA[:-1]))*100
    ESErrorSalmon = np.abs(np.divide(esAsSalmon[:-1]-salmonSA[:-1], salmonSA[:-1]))*100
    print(np.array([np.min(methodErrorSalmon), np.mean(methodErrorSalmon), np.max(methodErrorSalmon)]),
      np.array([np.min(partitionErrorSalmon), np.mean(partitionErrorSalmon), np.max(partitionErrorSalmon)]),
      np.array([np.min(ellipsoidErrorSalmon), np.mean(ellipsoidErrorSalmon), np.max(ellipsoidErrorSalmon)]),
      np.array([np.min(ESErrorSalmon), np.mean(ESErrorSalmon), np.max(ESErrorSalmon)]))
    
    methodErrorScanned = np.abs(np.divide(methodSalmon[-1]-salmonSA[-1], salmonSA[-1]))*100
    ellipsoidErrorScanned = np.abs(np.divide(ellipsoidSalmon[-1]-salmonSA[-1], salmonSA[-1]))*100
    partitionErrorScanned = np.abs(np.divide(partitionSalmon[-1]-salmonSA[-1], salmonSA[-1]))*100
    ESErrorScanned = np.abs(np.divide(esAsSalmon[-1]-salmonSA[-1], salmonSA[-1]))*100 
    print(np.array([np.min(methodErrorScanned), np.mean(methodErrorScanned), np.max(methodErrorScanned)]),
      np.array([np.min(partitionErrorScanned), np.mean(partitionErrorScanned), np.max(partitionErrorScanned)]),
      np.array([np.min(ellipsoidErrorScanned), np.mean(ellipsoidErrorScanned), np.max(ellipsoidErrorScanned)]),
      np.array([np.min(ESErrorScanned), np.mean(ESErrorScanned), np.max(ESErrorScanned)]))
      
    methodErrorCod = np.abs(np.divide(methodCod-codSA, codSA))*100
    ellipsoidErrorCod = np.abs(np.divide(ellipsoidCod-codSA, codSA))*100
    partitionErrorCod = np.abs(np.divide(partitionCod-codSA, codSA))*100
    ESErrorCod = np.abs(np.divide(esAsCod-codSA, codSA))*100
    print(np.array([np.min(methodErrorCod), np.mean(methodErrorCod), np.max(methodErrorCod)]),
      np.array([np.min(partitionErrorCod), np.mean(partitionErrorCod), np.max(partitionErrorCod)]),
      np.array([np.min(ellipsoidErrorCod), np.mean(ellipsoidErrorCod), np.max(ellipsoidErrorCod)]),
      np.array([np.min(ESErrorCod), np.mean(ESErrorCod), np.max(ESErrorCod)]))
    '''
    
    # plot everything
    fig1, ax1 = plt.subplots()
    ax1.plot(scannedLength, scannedAs, 'ko', label='Scanned Salmon Data')
    ax1.plot(OsheaSalmonLength, salmonSA[:-1], 'bo', label='O\'shea et al. data')
    ax1.plot(OsheaSalmonLength, salmonLRes, 'b', label='O\'shea et al. data regression line\n with 95% confidence interval')
    ax1.fill_between(OsheaSalmonLength, salmonLUncLower, salmonLUncUpper, alpha=0.2)
    ax1.plot(salmonLength, methodSalmon, 'k--', label='Developed Method')
    ax1.plot(salmonLength, partitionSalmon, 'g:', label='Rangtung et al. Partition Disc')
    ax1.plot(salmonLength, ellipsoidSalmon, 'g--', label='Rangtung et al. Prolate Spheroid')
    ax1.plot(salmonLength, esAsSalmon, 'r:', label='Harotounian et al. Prolate Spheroid')
    #ax1.set_title('Salmon')
    ax1.set_xlabel('Lengths [$cm$]')
    ax1.set_ylabel('Surface Area [$cm^2$]')
    ax1.legend()
    #fig1.savefig(os.path.join(figPath, 'salmonWithArea.pdf'), bbox_inches='tight')
    #fig1.savefig(os.path.join(figPath, 'salmonNoArea.pdf'), bbox_inches='tight')
    #fig1.savefig(os.path.join(figPath, 'salmonNoScale.pdf'), bbox_inches='tight')
   
    fig4, ax4 = plt.subplots()
    ax4.plot(OsheaCodLength, codSA, 'bo', label='O\'shea et al. data')
    ax4.plot(OsheaCodLength, codLRes, 'b', label='O\'shea et al. data regression line\n with 95% confidence interval')
    ax4.fill_between(OsheaCodLength, codLUncLower, codLUncUpper, alpha=0.2)
    ax4.plot(codLength, methodCod, 'k--', label='Developed Method')
    ax4.plot(codLength, partitionCod, 'g:', label='Rantung et al. Partition')
    ax4.plot(codLength, ellipsoidCod, 'g--', label='Rantung et al. Ellipsoid')
    ax4.plot(codLength, esAsCod, 'r:', label='Harotounian et al.')
    ax4.set_xlabel('Lengths [$cm$]')
    ax4.set_ylabel('Surface Area [$cm^2$]')
    #ax4.set_title('Cod')
    ax4.legend()
    #fig4.savefig(os.path.join(figPath, 'codWithArea.pdf'), bbox_inches='tight')
    #fig4.savefig(os.path.join(figPath, 'codNoArea.pdf'), bbox_inches='tight')

    plt.show()
    
if __name__ == '__main__':
    main()
