
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

def computeError(cPoints, dPoints):

    # Local Variables: i, cPoints, j, dPoints, error
    # Function calls: computeError, sqrt
    #% Compute L-2 norm error between cornerDetector and data set points
    #% The input arguments are:
    #%    cPoints = (x,y) points from cornerDetector
    #%    dPoints = (x,y) points from dataset
    #% The output argument is:
    #%    error = L-2 norm error sqrt((x_c - x_d)^2 + (y_c - y_d)^2)
    error = 0.
    for i in np.arange(1., 5.0):
        for j in np.arange(1., 3.0):
            error = error+(cPoints[int(i)-1,int(j)-1]-dPoints[int(i)-1,int(j)-1])**2.
            
        
    error = np.sqrt(error)
    return [error]