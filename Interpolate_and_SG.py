"""
Last modified date 11/16/2018

@author: Peco
"""
#SG smoothening function
def sgolay2d ( z, window_size, order, derivative=None):
    """
    This code applies SG in 2D images.
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return sp.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return sp.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return sp.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return sp.signal.fftconvolve(Z, -r, mode='valid'), sp.signal.fftconvolve(Z, -c, mode='valid')

#sgolay2d(Zn, window_size=boxsize, order=smoothening_order)
def eliminate_module_beamstop(detector_img):
    import pandas
    full_image = np.empty((516,1556))
    full_image[:] = np.nan
    delpixel = 6 #this is the width of a module. six pixels
    for i in range(6):
        #x scan -- bottom row
        if i == 0:
            full_image[0:261-delpixel,30+260*i:261-delpixel+260*i] = detector_img[0:261-delpixel,30+260*i:261-delpixel+260*i]
        else:
            full_image[0:261-delpixel,0+260*i:261-delpixel+260*i] = detector_img[0:261-delpixel,0+260*i:261-delpixel+260*i]

        #x scan -- top row
        if i == 0:
            full_image[260:517,30+260*i:261-delpixel+260*i] = detector_img[260:517,30+260*i:261-delpixel+260*i]
        else:
            full_image[260:517,0+260*i:261-delpixel+260*i] = detector_img[260:517,0+260*i:261-delpixel+260*i]
    image = pandas.DataFrame(full_image) #pandas is used to do interpoloation
    #now let's do interpoloation along vertical and horizontal. Horizontal needs two sweeping right to left and left to right.
    #interpolation is done such that intersity goes from high to low as you go right to left
    #except at the very left where intensity is kept constant.
    b = image.interpolate(method='linear', limit_direction='forward', axis=0) #along vertical
    c = b.interpolate(method='linear', limit_direction='forward', axis=1) #sweep to right along horizontal
    full_padded = c.interpolate(method='linear', limit_direction='backward', axis=1) #sweep to left along horizontal
    return(full_image, full_padded.as_matrix(), np.ndarray.tolist(full_image == 0))

def eliminate_module_beamstop_smootheachmodule(detector_img):
    """This isn't used. Different from the definition above in that each module is smoothened first
    and then put together into one detector image. Doesn't work well. Doesn't have nice continuity"""
    import pandas
    full_image = np.empty((516,1556))
    full_image[:] = np.nan
    delpixel = 6
    for i in range(6):
        #x scan -- bottom row
        if i == 0:
            full_image[0:261-delpixel,30+260*i:261-delpixel+260*i] = sgolay2d(detector_img[0:261-delpixel,30+260*i:261-delpixel+260*i],window_size=boxsize, order=smoothening_order, derivative=set_der)
        else:
            full_image[0:261-delpixel,0+260*i:261-delpixel+260*i] = sgolay2d(detector_img[0:261-delpixel,0+260*i:261-delpixel+260*i],window_size=boxsize, order=smoothening_order, derivative=set_der)

        #x scan -- top row
        if i == 0:
            full_image[260:517,30+260*i:261-delpixel+260*i] = sgolay2d(detector_img[260:517,30+260*i:261-delpixel+260*i],window_size=boxsize, order=smoothening_order, derivative=set_der)
        else:
            full_image[260:517,0+260*i:261-delpixel+260*i] = sgolay2d(detector_img[260:517,0+260*i:261-delpixel+260*i],window_size=boxsize, order=smoothening_order, derivative=set_der)
    return(full_image,np.ndarray.tolist(full_image == 0))

def correct_andSG(data_frame,window_size, order):
    """This is get rids of beam stops. Do interpoloation and then apply SG."""
    return(sgolay2d(eliminate_module_beamstop(data_frame)[1],window_size, order))