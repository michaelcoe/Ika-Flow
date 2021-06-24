import os
import sys
import re
import numpy as np
import pandas as pd

from pathlib import Path

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