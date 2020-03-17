import numpy as np

def interpolate(X,Y,x):

    idx = (np.abs(X - x)).argmin()
    if X[idx]==x:
        return Y[idx]
    elif X[idx]<x:
        x0 = X[idx]
        x1 = X[idx+1]
        y0 = Y[idx]
        y1 = Y[idx+1]
    else:
        x0 = X[idx-1]
        x1 = X[idx]
        y0 = Y[idx-1]
        y1 = Y[idx]
    return y0 + (y1-y0)/(x1-x0)*x
