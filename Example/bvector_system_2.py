def bvector_system_2(x):
    """
    b vector definition in Task 4.

    Parameters
    ----------
    x : float
        Coordinate

    Returns
    -------
    b : array of float
        b as given by equation (13)
    """

    # Define vector b according to equation (13)
    b=np.array([
        np.cos(10*x)-10*np.sin(10*x),
        199*np.cos(10*x)-10*np.sin(10*x),
        208*np.cos(10*x)+10000*np.sin(10*x)
    ])
    # Return vector b
    return b
