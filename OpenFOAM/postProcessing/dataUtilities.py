import numpy as np

### simple high pass fliter function
def filterData(x, kernelLength = 11, kernelFunction = 'flat'):

    if len(x) < kernelLength:
        raise ValueError("kernel length > data")
    
    if not kernelFunction in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("kernel function available are: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    n = np.arange(kernelLength)
    h = np.sinc(2 * 0.1 * (n - (kernelLength - 1) / 2))
    
    # flat corresponds to a moving average filter
    if kernelFunction == "flat":
        kernel = np.ones(kernelLength)
    else:
        kernel = eval('np.' + kernelFunction + '(' + str(kernelLength) + ')')

    h = h * kernel
    h = h / np.sum(h)
    h = -h
    h[(kernelLength - 1) // 2] += 1
    
    return np.convolve(x, h, mode='valid')