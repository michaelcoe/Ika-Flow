import numpy as np

from pathlib import Path
from dataUtilities import filterData

class ForceBins:

    def __init__(self,
                 inputpath,
                 cycles = 3.0,
                 total_cycles = 3.0,
                 average = True,
                 filterForces = True):

        self.force_path = Path(inputpath).parent.joinpath('forceBin.dat')
        self.moment_path = Path(inputpath).parent.joinpath('momentBin.dat')
        self.specific_case = self.force_path.parts[-6]
        self.parent_case = self.force_path.parts[-7]
        self.cycles = cycles
        self.total_cycles = total_cycles

        # all forces should be loaded by now
        # build a "nice" dict with the forces                    
   
        self.forces = dict()
        self.moments = dict()

        self.forceCoord_x, self.forceCoord_y, self.forceCoord_z, _rawForces = self._readForceBinFile(self.force_path)
        self.momentCoord_x, self.momentCoord_y, self.momentCoord_z, _rawMoments = self._readForceBinFile(self.moment_path)

        self.number_coordinates = len(self.forceCoord_x)
        pos = iter(range(1,10*self.number_coordinates))
        
        self.forces["time"] = self.moments["time"] = _rawForces[:,0]
        for num in range(self.number_coordinates):
            self.forces[num] = {}
            self.moments[num] = {}
            for forceType in ("total", "pressure", "viscous"):
                self.forces[num][forceType] = {}
                self.moments[num][forceType] = {}
                for component in "x", "y", "z":
                    currentPos = next(pos)
                    self.forces[num][forceType][component] = _rawForces[:,currentPos]
                    self.moments[num][forceType][component] = _rawMoments[:,currentPos]               
        if average:
            self.calculateAverageStd()
        if filterForces:
            self.filterForcesMoments()
            self.calculateFilteredAverageStd()

    # function to process force.dat files
    def _readForceBinFile(self, file_name):
        raw = []
        x_coords = []
        y_coords = []
        z_coords = []

        with open(file_name, 'r') as f:
            for line in f:
                tmp = [x.strip('(').strip(')') for x in line.split()]
                if len(tmp) == 0:
                    continue
                elif tmp[0] == '#' and len(tmp) == 1:
                    continue
                elif tmp[1] == "x":
                    data = tmp[4:]
                    x_coords.append([ float(i) for i in data ])
                elif tmp[0] == '#' and tmp[1] == 'y':
                    data = tmp[4:]
                    y_coords.append([ float(i) for i in data ])
                elif tmp[0] == '#' and tmp[1] == 'z':
                    data = tmp[4:]
                    z_coords.append([ float(i) for i in data ])
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
        x_coords = np.array(x_coords)
        y_coords = np.array(y_coords)
        z_coords = np.array(z_coords)

        return x_coords[0], y_coords[0], z_coords[0], raw

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

        for num in range(self.number_coordinates):
            self.averageForces[num] = {}
            self.averageMoments[num] = {}
            self.stdForces[num] = {}
            self.stdMoments[num] = {}
            for forceType in ("total", "pressure", "viscous"):
                self.averageForces[num][forceType] = {}
                self.averageMoments[num][forceType] = {}
                self.stdForces[num][forceType] = {}
                self.stdMoments[num][forceType] = {}
                for component in ("x", "y", "z"):
                    self.averageForces[num][forceType][component] = np.average(self.forces[num][forceType][component][mask])
                    self.averageMoments[num][forceType][component] = np.average(self.moments[num][forceType][component][mask])
                    self.stdForces[num][forceType][component] = np.std(self.forces[num][forceType][component][mask])
                    self.stdMoments[num][forceType][component] = np.std(self.moments[num][forceType][component][mask])

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
        
        for num in range(self.number_coordinates):
            self.filteredForces[num] = {}
            self.filteredMoments[num] = {}
            for forceType in ("total", "pressure", "viscous"):
                self.filteredForces[num][forceType] = {}
                self.filteredMoments[num][forceType] = {}
                for component in ("x", "y", "z"):
                    self.filteredForces[num][forceType][component] = filterData(self.forces[num][forceType][component][mask], filterWindow, filterFunction)
                    self.filteredMoments[num][forceType][component] = filterData(self.moments[num][forceType][component][mask], filterWindow, filterFunction)

        return self.filteredForces, self.filteredMoments

    # Calculates the average and standard deviation on filtered data
    def calculateFilteredAverageStd(self):

        if hasattr(self, "filteredForces") == False:
            raise Exception("missing attribute filteredForces. Please run filterForces prior to calculateFilteredAveragesStd!")
        
        self.averageFilteredForces = {}
        self.stdFilteredForces = {}

        self.averageFilteredMoments = {}
        self.stdFilteredMoments = {}

        for num in range(self.number_coordinates):
            self.averageFilteredForces[num] = {}
            self.stdFilteredForces[num] = {}
            self.averageFilteredMoments[num] = {}
            self.stdFilteredMoments[num] = {}
            for forceType in ("total", "pressure", "viscous"):
                self.averageFilteredForces[num][forceType] = {}
                self.averageFilteredMoments[num][forceType] = {}
                self.stdFilteredForces[num][forceType] = {}
                self.stdFilteredMoments[num][forceType] = {}
                for component in ("x", "y", "z"):
                    # calculate average forces
                    self.averageFilteredForces[num][forceType][component] = np.average(self.filteredForces[num][forceType][component])
                    self.stdFilteredForces[num][forceType][component] = np.std(self.filteredForces[num][forceType][component])
                    # calculate average moments
                    self.averageFilteredMoments[num][forceType][component] = np.average(self.filteredMoments[num][forceType][component])
                    self.stdFilteredMoments[num][forceType][component] = np.std(self.filteredMoments[num][forceType][component])


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