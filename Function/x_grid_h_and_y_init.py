###############################################################
##             Define an initializaton Function              ##
###############################################################

def x_grid_h_and_y_init(y0,interval,N):
    """
    Solve the IVP y' = A y + b, y(0) = y0, in the interval,
    This function is used to create the x grid with step h and initialize matrix y to store the result
    Parameters
    ----------
    y0 : vector
        Initial data
    interval : vector
        Interval on which solution is required
    N : int
        Number of steps

    Returns
    -------
    x : array of float
        Coordinate locations of the approximate solution
    h : int
        Stepsize
    y : array of float
        Values of approximate solution at locations x
    """
    
    #Start
    
    #To solve an IVP numerically, we require an x-grid of N+1 points in the interval
    x=np.linspace(interval[0],interval[1],N+1)

    #Cal the stepsize h
    h=x[1]-x[0]

    #We can create a solution matrix to record the whole process
    y=np.zeros([N+1,len(y0)])
    #The col of y corresponds to a single y ODE
    #The row of y corresponds to a x gridpoint

    return x,h,y
    #End