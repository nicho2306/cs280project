
import numpy as np
import scipy
import matcompat

# if available import pylab (from matlibplot)
try:
    import matplotlib.pylab as plt
except ImportError:
    pass

#% Author: Nicholas Sunjaya
#% CS280
#% HW1
#% Problem 3.4
def homography_solve(u, v):

    # Local Variables: A, i, H, v, N, u, V, h
    # Function calls: reshape, homography_solve, size
    N = matcompat.size(v, 2.)
    V = np.reshape(v, np.array(np.hstack((2.*N, 1.))))
    A = np.array([])
    for i in np.arange(1., (N)+1):
        A = np.array(np.vstack((np.hstack((A)), np.hstack((u[0,int(i)-1], u[1,int(i)-1], 1., 0., 0., 0., np.dot(-u[0,int(i)-1], v[0,int(i)-1]), np.dot(-u[1,int(i)-1], v[0,int(i)-1]))), np.hstack((0., 0., 0., u[0,int(i)-1], u[1,int(i)-1], 1., np.dot(-u[0,int(i)-1], v[1,int(i)-1]), np.dot(-u[1,int(i)-1], v[1,int(i)-1]))))))
        
    h = linalg.solve(A, V)
    H = np.array(np.vstack((np.hstack((h[0], h[1], h[2])), np.hstack((h[3], h[4], h[5])), np.hstack((h[6], h[7], 1.)))))
    return [H]