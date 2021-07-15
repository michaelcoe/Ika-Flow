import numpy as np
from os.path import isfile
import subprocess as sb
import re
from os import listdir
# import matplotlib.pyplot as plt
import math

toCoefficient = lambda rho, u, A: 0.5*rho*(u**2)*A

def isNumber(s):
    try:
        complex(s) # for int, long, float and complex
    except ValueError:
        return False
    return True

isTimeDir = re.compile(r'(^\d+\.?\d*$)')

def getTimeDirs(path):
    timeDirs = []
    for dir in listdir(path):
        if isTimeDir.search(dir):
            timeDirs.append(dir)
    if len(timeDirs) == 0:
        return None
    else:
        return timeDirs

            

def getIndices(data, start = 0, end = 0):
    if start > end:
        startIndex = 0
        while data[startIndex] < start:
            startIndex += 1
        return [startIndex, len(data)]
        
    elif start == 0 and end == 0:
        return [0, len(data)]
    else:   
        startIndex = 0
        endIndex = 0
        while data[startIndex] < start: 
            startIndex += 1
        while data[endIndex] < end and endIndex < len(data)-1:
            endIndex += 1        
        return [int(startIndex), int(endIndex + 1)]
        
### define a simple function that reads a line and tries to convert to an array of floats
def readLine(line):
    
    try:
        elements = line.split()
        if elements[0] == '#':
            return None
        if len(elements) > 1:
            elements = [float(x) for x in elements]
        else:
            return float(elements[0])
    except:
        return None

    return elements
    
def readFile(filepath):
    try:
        with open(filepath, 'r') as f:
            raw = []
            for line in f:
                data = readLine(line)
                if data != None and len(data) > 1:
                    raw.append(data)
            return np.array(raw)
            
    except (IOError, OSError) as e:
        print("An Error occured while trying to read file {}".format(filepath))
        print("error({0}): {1}".format(e.errno, e.strerror))
    

def readTimeFile(filepath, startTime = 0.1, endTime = 10, timeAxis = 0):
    
    try:
        filehandle = open(filepath, 'r')
    except (IOError, OSError) as e:
        print ("could not open file: ", filepath)
        print ("Exception: ", e)
        
    raw = []
    
    for line in filehandle:
        data = readLine(line)
        if data != None:
            if len(data) > 1:
                if data[timeAxis] > startTime and data[timeAxis] <= endTime:
                    raw.append(data)
            else:
                raw.append(data)
    
    filehandle.close()
    
    return np.array(raw)

    

# define a simple function to read in file linewise
# and return force coefficients
# if no specific fluid & geometry parameters are provided, the original forces are returned
###
# OpenFOAM returns forces in the following order:
# Time       forces(pressure viscous porous) moment(pressure viscous porous)
def readForceFile (filepath, startTime = 0.1, endTime = 10, rho = 1.0, u = 1.0, A = 2):
    raw = []
    
    filehandle = open(filepath, 'r')
    for line in filehandle:
        tmp = [x.strip('(').strip(')') for x in line.split()]
        if len(tmp) == 0:
            continue
        elif tmp[0] == '#':
            continue
        elif len(tmp) != 19:
            continue
        else:
            if float(tmp[0]) < startTime:
                continue
            elif float(tmp[0]) > endTime:
                break
            else:
                # [float(i) for i in lst]
                # raw.append(list(map(float, tmp)))
                try:
                    raw.append([ float(i) for i in tmp ])
                except:
                    print("could not convert string to float in line:")
                    print("\t" + line)
                    print("in file:")
                    print("\t" + filepath)
    
    filehandle.close()
    raw = np.array(raw)
    return np.array([ raw[:,0], ( raw[:,1]  / (0.5 * rho * u**2 * A )), (raw[:,2]  / (0.5 * rho * u**2 * A )), (raw[:,3]  / (0.5 * rho * u**2 * A )), (raw[:,4]  / (0.5 * rho * u**2 * A )), (raw[:,5]  / (0.5 * rho * u**2 * A )),(raw[:,6]  / (0.5 * rho * u**2 * A )) ])


def readForceFile2 (filepath, startTime = 0, endTime = 0):
    raw = []

    with open(filepath, 'r') as filehandle:
        for line in filehandle:
            tmp = [x.strip('(').strip(')') for x in line.split()]
            if len(tmp) == 0:
                continue
            elif tmp[0] == '#':
                continue
            elif len(tmp) != 10:
                continue
            else:
                if float(tmp[0]) < startTime:
                    continue
                elif float(tmp[0]) > endTime:
                    break
                else:
                # [float(i) for i in lst]
                # raw.append(list(map(float, tmp)))
                    try:
                        raw.append([ float(i) for i in tmp ])
                    except ValueError:
                        print("could not convert string to float in line:")
                        print("\t" + line)
                        print("in file:")
                        print("\t" + filepath)

    filehandle.close()
    raw = np.array(raw).T
    return raw



## simple hann window filter
# TODO does not work properly on numpy arrays, needs to be improved
def hann_filter (g, n):
    for i in range(0,len(g)):
        g[i] = g[i] * 0.5 * (1.0 - math.cos((2.0 * math.pi * i) / (n - 1)))

# simple function that finds the power of 2 that is next to the provided value
def nextpow2 (i):
    n = 1
    while n < i: n *= 2
    return n
    
### runs a simple fft analysis on arbitrary data
def fftAnalysis (time, data, printPeaks=False):

    N = len(time)
   
    if N != len(data):
        print("number of elements and time steps must match!")
        return
    elif N <= 1:
        return

    T = (max(time) - min(time))/(N-1)
    f = np.linspace(0.0, 1.0/(2*T), N/2)
    average = np.average(data)

    filteredData = []

    for i in range(0, N):
        filteredData.append((data[i] - average) * 0.5 * (1.0 - np.cos((2.0*np.pi*i) / (N - 1))))

    fftDistribution = np.abs(np.fft.fft(filteredData))[0:int(N/2)]
    peaks = [f[index] for index in fftDistribution.argsort()[-20:]]
    peakDistro = [fftDistribution[index] for index in fftDistribution.argsort()[-20:]]

    peaks.reverse()
    peakDistro.reverse()
   
    if printPeaks == True:
        print("# Frequency | Amplitude")
        for freq, amp in zip(peaks, peakDistro):
            print(freq, "|", amp)

#    return np.array([np.array(peaks), np.array(peakDistro), fftDistribution])
    return np.array([np.array(peaks), np.array(peakDistro), np.array(f), fftDistribution])
    
# TODO es fehlt scalar probes und vector probes!
# for probes
def readProbesFile (filepath, start_time = 1, end_time = 10):
    raw = []
    probes = {}

    x_c = []
    y_c = []

    
    f = open(filepath, 'r')
    print("opened file", f)
    for line in f:
        tmp = line.split()
        if len(tmp) <= 1:
            continue
        if tmp[0] == '#' and tmp[1] == 'Probe' and len(tmp) == 6:       # header of the probes file
            # probes[tmp[2]] = [float(tmp[3].strip('(')), float(tmp[4]), float(tmp[5].strip(')'))]
            x_c.append(float(tmp[3].strip('(')))
            y_c.append(float(tmp[4]))
        elif tmp[0] == '#':
            continue
        else:
            raw.append([float(i) for i in tmp])
    
    print("length raw: ", len(raw))
    
    f.close()    
    
    return x_c, y_c, raw
        # else:
  
# FIX add tempfile support
# FIX parse the output file in python (not grep and awk)  
def extractPressureGradient(filepath):
    if isfile(filepath):
        extractTimeCommand = "cat " + str(filepath) + " | grep \"^Time\" | awk \'{print $3}\' > timeSteps"
        extractPressureCommand = "cat " + str(filepath) + " | grep \"pressure gradient\" | awk \'{print $11}\' > rawPressureGradient"
        
        try:
            print (extractTimeCommand)
            sb.call(extractTimeCommand, shell = True)
        except:
            print ("could not execute extractTimeCommand")
            print (extractTimeCommand)
        try:
            print (extractPressureCommand)
            sb.call(extractPressureCommand, shell = True)
        except:
            print ("could not execeute extractPressureCommand")
            print (extractPressureCommand)
    else:
        raise Exception("no such file or directory: ", filepath)

    if isfile("./timeSteps") and isfile("./rawPressureGradient"):
        
        times = readFile("./timeSteps")
        pressure = readFile("./rawPressureGradient")
        
        if (len(pressure) % len(times) != 0):
            raise Exception("time steps and pressure gradient steps do not match")
        
        nPressurePerTime = len(pressure) / len(times)
        averagedPressure = []        
        
        for i in np.arange(nPressurePerTime, len(pressure) + 1, nPressurePerTime):
            averagedPressure.append(np.average(pressure[i-6:i]))

        try:
            fOut = open("pressureGradient", 'w')
        except:
            raise Exception("could not open out file")

        fOut.writelines("# timestep aveargedPressureGradient")
        for t, p in zip(times, averagedPressure):
            fOut.writelines(str(t) + " " + str(p) + "\n")
        
        fOut.close()
            
### simple fliter function
def filterData(x, kernelLength = 11, kernelFunction = 'flat'):

    if len(x) < kernelLength:
        raise ValueError("kernel length > data")
    
    if not kernelFunction in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("kernel function available are: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    # x = np.r_[x[kernelLength-1:0:-1], x, x[-2:-kernelLength-1:-1]]
    
    # flat corresponds to a moving average filter
    if kernelFunction == "flat":
        kernel = np.ones(kernelLength, 'd')
    else:
        kernel = eval('np.' + kernelFunction + '(kernelLength)')
    
    return np.convolve(kernel / kernel.sum(), x, mode='valid')




        
