def tred(t,x,vred):
    '''
    VJS 9/2014
    Apply a reduction velocity to the time axis
    Usage: t_red = tred(t,x,vred)
    Input:
        t: time
        x: source-receiver ranges
        vred: reduction velocity to apply
    Output:
        t_red:  time with reduced velocity
    '''
    
    import numpy as np

    t_red = t - np.abs(x)/vred
    return t_red