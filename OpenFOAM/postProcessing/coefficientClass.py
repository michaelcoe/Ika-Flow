import numpy as np

from pathlib import Path
from dataUtilities import filterData

class Coefficients:

    def __init__(self,
                 inputpath,
                 cycles = 3.0,
                 total_cycles = 4.0,
                 average = True,
                 filterForces = True,                 
                 filterType = 'flat',
                 filterWindow = 11):

        self.coeff_path = Path(inputpath).parent.joinpath('coefficient.dat')
        self.specific_case = self.coeff_path.parts[-6]
        self.parent_case = self.coeff_path.parts[-7]
        self.cycles = cycles
        self.total_cycles = total_cycles

        # all forces should be loaded by now
        # build a "nice" dict with the forces
        self.coefficients = dict()

        _rawCoefficients = self._readCoefficientFile(self.coeff_path)

        self.coefficients["time"] = _rawCoefficients[:,0]
        self.coefficientTypes = ['Cd', 'Cs', 'Cl', 'CmRoll', 'CmPitch', 'CmYaw', 'Cdf', 'Cdr', 'Csf', 'Csr', 'Clf', 'Clr']
        
        for i, coeffType in enumerate(self.coefficientTypes):
            self.coefficients[coeffType] = {}
            self.coefficients[coeffType]= _rawCoefficients[:,i+1]              
    
        if average:
            self.calculateAverageStd()
        if filterForces:
            self.filterCoefficients()
            self.calculateFilteredAverageStd()

    # function to process force.dat files
    def _readCoefficientFile(self, file_name):
        raw = np.loadtxt(file_name, comments='#', skiprows=13)
        return raw

    # Returns an indices mask based based on the number of cycles that want to be plotted
    def _getIndices(self):
        cuttoff_time = self.coefficients['time'][-1] * ((self.total_cycles-self.cycles)/self.total_cycles)
        return np.where(self.coefficients['time'] >= cuttoff_time, True, False)
    
    def _getIndicesByTime(self, dictType, startTime, endTime):
            return np.logical_and(self.coefficients['time'] >= startTime, self.coefficients['time'] <= endTime)
    
    # calculates the average and standard deviation on unfiltered data
    def calculateAverageStd(self):

        self.averageCoefficients = {}
        self.stdCoefficients = {}

        mask = self._getIndices()

        for i, coeffType in enumerate(self.coefficientTypes):
            self.averageCoefficients[coeffType] = {}
            self.stdCoefficients[coeffType] = {}
            self.averageCoefficients[coeffType] = np.average(self.coefficients[coeffType][mask])
            self.stdCoefficients[coeffType] = np.std(self.coefficients[coeffType][mask])
        
        return { 'coefficients' : { "average" : self.averageCoefficients, "std" : self.stdCoefficients}}
    
    # filters the data
    def filterCoefficients(self, filterFunction = "flat", filterWindow = 11):
        if filterWindow % 2 == 0:
            raise Exception("filterWindow needs to be an uneven number!")

        mask = self._getIndices()
        endTimeIndex = int(len(self.coefficients["time"][mask]) - ((filterWindow - 1)/2))

        self.filteredCoefficients = {}
        self.filteredCoefficients["time"] =  self.coefficients["time"][int((filterWindow - 1)/2):endTimeIndex]

        for i, coeffType in enumerate(self.coefficientTypes):
            self.filteredCoefficients[coeffType] = {}
            self.filteredCoefficients[coeffType]= filterData(self.coefficients[coeffType][mask], filterWindow, filterFunction)

        return self.filteredCoefficients

    # Calculates the average and standard deviation on filtered data
    def calculateFilteredAverageStd(self):

        if hasattr(self, "filteredCoefficients") == False:
            raise Exception("missing attribute filteredForces. Please run filterForces prior to calculateFilteredAveragesStd!")
        
        self.averageFilteredCoefficients = {}
        self.stdFilteredCoefficients = {}

        for i, coeffType in enumerate(self.coefficientTypes):
            self.averageFilteredCoefficients[coeffType] = {}
            self.stdFilteredCoefficients[coeffType] = {}
            # calculate average forces
            self.averageFilteredCoefficients[coeffType] = np.average(self.filteredCoefficients[coeffType])
            self.stdFilteredCoefficients[coeffType] = np.std(self.filteredCoefficients[coeffType])

        return { 'filteredCoefficients' : { "average" : self.averageFilteredCoefficients, "std" : self.stdFilteredCoefficients}}

    def getCoefficientsMinTime(self):
        print("min time is {}".format(self.coefficients["time"][0]))
        return self.coefficients["time"][0]

    ## define a method for getting forces by time
    def getCoefficientsByTime(self,  startTime = 0, endTime = 0, coeffType = "Cd"):
        mask = self._getIndicesByTime(startTime, endTime)
        return self.coefficients[coeffType][mask]