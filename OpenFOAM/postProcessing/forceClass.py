import numpy as np

from pathlib import Path
from dataUtilities import filterData

class Forces:

    def __init__(self,
                 inputpath,
                 cycles = 2.0,
                 total_cycles = 3.0,
                 average = True,
                 filterForces = True):

        self.force_path = Path(inputpath).parent.joinpath('force.dat')
        self.moment_path = Path(inputpath).parent.joinpath('moment.dat')
        self.specific_case = self.force_path.parts[-6]
        self.parent_case = self.force_path.parts[-7]
        self.cycles = cycles
        self.total_cycles = total_cycles

        # all forces should be loaded by now
        # build a "nice" dict with the forces                    
        pos = iter(range(1,10))
        self.forces = dict()
        self.moments = dict()

        _rawForces = self._readForceFile(self.force_path)
        _rawMoments = self._readForceFile(self.moment_path)
        self.forces["time"] = self.moments["time"] = _rawForces[:,0]
        for forceType in ("total", "pressure", "viscous"):
            self.forces[forceType] = {}
            self.moments[forceType] = {}
            for component in "x", "y", "z":
                currentPos = next(pos)
                self.forces[forceType][component] = _rawForces[:,currentPos]
                self.moments[forceType][component] = _rawMoments[:,currentPos]               
    
        if average:
            self.calculateAverageStd()
        if filterForces:
            self.filterForcesMoments()
            self.calculateFilteredAverageStd()

    # function to process force.dat files
    def _readForceFile(self, file_name):
        raw = []

        with open(file_name, 'r') as f:
            for line in f:
                tmp = [x.strip('(').strip(')') for x in line.split()]
                if len(tmp) == 0:
                    continue
                elif tmp[0] == '#':
                    continue
                else:
                    try:
                        raw.append([ float(i) for i in tmp ])
                    except:
                        print("could not convert string to float in line:")
                        print("\t" + line)
                        print("in file:")
                        print("\t" + file_name)

        raw = np.array(raw)

        return raw

    # Returns an indices mask based based on the number of cycles that want to be plotted
    def _getIndices(self):
        cuttoff_time = self.forces['time'][-1] * ((self.total_cycles-self.cycles)/self.total_cycles)
        return np.where(self.forces['time'] >= cuttoff_time, True, False)
    
    def _getIndicesByTime(self, dictType, startTime, endTime):
        if dictType == 'forces':
            return np.logical_and(self.forces['time'] >= startTime, self.forces['time'] <= endTime)
        else:
            return np.logical_and(self.moments['time'] >= startTime, self.moments['time'] <= endTime)
    
    # calculates the average and standard deviation on unfiltered data
    def calculateAverageStd(self):

        self.averageForces = {}
        self.stdForces = {}

        self.averageMoments = {}
        self.stdMoments = {}

        mask = self._getIndices()

        for forceType in ("total", "pressure", "viscous"):
            self.averageForces[forceType] = {}
            self.averageMoments[forceType] = {}
            self.stdForces[forceType] = {}
            self.stdMoments[forceType] = {}
            for component in ("x", "y", "z"):
                self.averageForces[forceType][component] = np.average(self.forces[forceType][component][mask])
                self.averageMoments[forceType][component] = np.average(self.moments[forceType][component][mask])
                self.stdForces[forceType][component] = np.std(self.forces[forceType][component][mask])
                self.stdMoments[forceType][component] = np.std(self.moments[forceType][component][mask])

        return {"forces" : { "average" : self.averageForces, "std" : self.stdForces },
                "moments" : { "average" : self.averageMoments, "std" : self.stdMoments} }

    # filters the data
    def filterForcesMoments(self, filterFunction = "flat", filterWindow = 11):
        if filterWindow % 2 == 0:
            raise Exception("filterWindow needs to be an uneven number!")

        mask = self._getIndices()
        endTimeIndex = int(len(self.forces["time"][mask]) - ((filterWindow - 1)/2))

        self.filteredForces = {}
        self.filteredMoments = {}
        self.filteredForces["time"] =  self.forces["time"][int((filterWindow - 1)/2):endTimeIndex]
        self.filteredMoments["time"] =  self.moments["time"][int((filterWindow - 1)/2):endTimeIndex]

        for forceType in ("total", "pressure", "viscous"):
            self.filteredForces[forceType] = {}
            self.filteredMoments[forceType] = {}
            for component in ("x", "y", "z"):
                self.filteredForces[forceType][component] = filterData(self.forces[forceType][component][mask],                                                                          filterWindow, filterFunction)
                
                self.filteredMoments[forceType][component] = filterData(self.moments[forceType][component][mask],                                                                        filterWindow, filterFunction)

        return (self.filteredForces, self.filteredMoments)

    # Calculates the average and standard deviation on filtered data
    def calculateFilteredAverageStd(self):

        if hasattr(self, "filteredForces") == False:
            raise Exception("missing attribute filteredForces. Please run filterForces prior to calculateFilteredAveragesStd!")
        
        self.averageFilteredForces = {}
        self.stdFilteredForces = {}

        self.averageFilteredMoments = {}
        self.stdFilteredMoments = {}

        for forceType in ("total", "pressure", "viscous"):
            self.averageFilteredForces[forceType] = {}
            self.averageFilteredMoments[forceType] = {}
            self.stdFilteredForces[forceType] = {}
            self.stdFilteredMoments[forceType] = {}
            for component in ("x", "y", "z"):
                # calculate average forces
                self.averageFilteredForces[forceType][component] = np.average(self.filteredForces[forceType][component])
                self.stdFilteredForces[forceType][component] = np.std(self.filteredForces[forceType][component])
                # calculate average moments
                self.averageFilteredMoments[forceType][component] = np.average(self.filteredMoments[forceType][component])
                self.stdFilteredMoments[forceType][component] = np.std(self.filteredMoments[forceType][component])


        return { "forces" : { "average" : self.averageFilteredForces, "std" : self.stdFilteredForces},
                 "moments": { "average" : self.averageFilteredMoments, "std" : self.stdFilteredMoments}}

    def convertToCoefficient(self):
        pass

    def getForcesMinTime(self):
        print("min time is {}".format(self.forces["time"][0]))
        return self.forces["time"][0]

    def getMomentsMinTime(self):
        print("min time is {}".format(self.moments["time"][0]))
        return self.moments["time"][0]

    ## define a method for getting forces by time
    def getForcesByTime(self,  startTime = 0, endTime = 0, forceType = "total", forceComponent = "x"):
        mask = self._getIndicesByTime('forces', startTime, endTime)
        return self.forces[forceType][forceComponent][mask]

    ## define a method for getting moments by time
    def getMomentsByTime(self,  startTime = 0, endTime = 0, forceType = "total", forceComponent = "x"):
        mask = self._getIndicesByTime('moments', startTime, endTime)
        return self.moments[forceType][forceComponent][mask]