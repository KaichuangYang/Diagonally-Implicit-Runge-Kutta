def dirk3(A, bvector, y0, interval, N):
    """
    Solve the IVP y' = A y + b, y(0) = y0, in the interval,
    using N steps of DIRK3.

    Parameters
    ----------
    A : matrix
        Partially defines ODE according to (4)
    bvector : function name returning vector
        Completes definition of ODE according to (4)
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
    y : array of float
        Values of approximate solution at locations x
    """
    
    # Add DIRK3 algorithm implementation here according to Task 2

    #Start

    #Initialization
    x,h,y=x_grid_h_and_y_init(y0,interval,N)
    y[0,:]=y0 #Update y with given initial val y0
    mu=(1/2)*(1-1/np.sqrt(3))
    nu=(1/2)*(np.sqrt(3)-1)
    gamma=3/(2*(3+np.sqrt(3)))
    lamb_da=(3*(1+np.sqrt(3)))/(2*(3+np.sqrt(3)))
    I=np.identity(np.shape(A)[0])

    #Implicit algoritm
    for i in range(0,N):
        x_n=x[i]
        y_n=y[i,:]
        #solve AX=b to get y_1
        A_1=I-h*mu*A
        b_1=y_n+h*mu*bvector(x_n+h*mu)
        y_1=np.linalg.solve(A_1,b_1)
        #solve AX=b to get y_2
        A_2=I-h*mu*A
        b_2=y_1+h*nu*(np.dot(A,y_1)+bvector(x_n+h*mu))+h*mu*bvector(x_n+h*nu+2*h*mu)
        y_2=np.linalg.solve(A_2,b_2)
        #update y
        y[i+1,:]=(1-lamb_da)*y_n+lamb_da*y_2+h*gamma*(np.dot(A,y_2)+bvector(x_n+h*nu+2*h*mu))
    
    #Since we require shape of y be n x (N+1)
    #Notice that we can def y in fuc x_grid_h_and_y_init as n x (N+1), so that we do not need to transpose here
    #In this case, please also change for loop update!
    #But I follow the 1 st lab session's slides to def y as in the slides
    y=y.T

    # The first output argument x is the locations x_j at which the solution is evaluated;
    # this should be a real vector of length N + 1 covering the required interval.
    # The second output argument y should be the numerical solution approximated at the
    # locations x_j, which will be an array of size n × (N + 1).
    return x, y
#End