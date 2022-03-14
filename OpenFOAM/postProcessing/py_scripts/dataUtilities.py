from logging import Filter
import numpy as np

### simple high pass fliter function
def filterData(x, kernelLength = 11, kernelFunction = 'flat'):

    if len(x) < kernelLength:
        raise ValueError("kernel length > data")
    
    if not kernelFunction in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("kernel function available are: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
        
    # flat corresponds to a moving average filter
    if kernelFunction == "flat":
        kernel = np.ones(kernelLength)
    else:
        kernel = eval('np.' + kernelFunction + '(' + str(kernelLength) + ')')
        
        kernel /= np.sum(kernel)
    
    return np.convolve(x, kernel, mode='valid')