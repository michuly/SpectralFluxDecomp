import numpy as np
from scipy import ndimage
from scipy.signal import butter, lfilter, sosfiltfilt

COMPENSATING_FACTOR = np.sqrt(8 / 3)  # for hanning window


def reflect_signal(data, axes, flags, normalize=False, verbose=False):

    shape = np.shape(data)
    dim = np.ndim(data)

    if len(flags)!=len(axes):
        print('!!!not enough flags for given axes!!!')
    if verbose:
        print('reflecting for fft... ')
    for i in range(len(axes)):
        if verbose and normalize:
            print(0, np.sum(data * data))

        flag = flags[i]
        axis = axes[i]
        N = shape[axis]
        if axis==2 and flag!=0:
            print('Is this axis anf flag correct?', axis, flag)

        if flag == 1 or flag == -1: # -1 for sine (antisymmetric), 1 for cosine (symmetric)

            slicer1 = [slice(None)] * dim  # building slicer
            slicer2 = [slice(None)] * dim
            slicer1[axis] = slice(0, N-1)
            slicer2[axis] = slice(1, N)
            # slicer2[axis] = slice(0, N)
            data = np.concatenate((data[tuple(slicer1)], flag * np.flip(data[tuple(slicer2)], axis=axis)), axis=axis)
            # data = data/np.sqrt(2)  # NORMALIZING FACTOR, to keep energy. need to check for cosine.

        elif flag == 0:

            data = np.concatenate((data, data), axis=axis) # because of n=512 instead of 513
            # data = np.concatenate((data, data, data[:,0:1]), axis=axis) # because of n=512 instead of 513
            # data = data/np.sqrt(2)  # NORMALIZING FACTOR, to keep energy. need to check for cosine.

        elif flag == 2:
            slicer1 = [slice(None)] * dim  # building slicer
            slicer2 = [slice(None)] * dim
            slicer1[axis] = slice(0, N-1)
            slicer2[axis] = slice(1, N)
            data = np.concatenate((data[tuple(slicer1)], data[tuple(slicer2)]), axis=axis)

        # normalize=True
        # print('ERROR - I AM NOW NORMALIZING ALTHOUGH I HAVE NEVER DONE THIS. WE WILL SEE IF THIS WORKS.')
        if normalize:
            if verbose:
                print(1, np.sum(data*data))
            data = data / np.sqrt(2)  # NORMALIZING FACTOR, to keep energy. need to check for cosine.
            if verbose:
                print(2, np.sum(data * data))

    return data


def derivative(y, L, dx=1, n=1, axis=1, flag=1, to_slice=True):
    """
    gets an matrix y of physical length L_i, physical interval dx_i (and N_i points) for each dimension.
    the function returns the nth derivative (DD), for the axis given.
    flag is 1 if the function is even (derivative using cosine transform) and
    -1 if it's odd (derivative using sine transform). If flag=2 the derivative
    is done using np.gradient. If flag=0 derivative is done for a
    periodic function using simple fft. Make sure in this case that the L you
    pass is that of the entire domain.
    """
    N=np.shape(y)[axis]

    # I can write flags 0,1,-1,2 together.
    if flag == 0:
        N_new = np.shape(y)[axis]
        y = reflect_signal(y, axes=(axis,), flags=(flag,)) # is it necessary to reflect?
        N2 = np.shape(y)[axis]
        y = spectral_derivative(y, axis, n, L, N2, N_new)

    elif flag == 1 or flag == -1:
        N_new = np.shape(y)[axis]
        y = reflect_signal(y, axes=(axis,), flags=(flag,))
        N2 = np.shape(y)[axis]
        y = spectral_derivative(y, axis, n, L, N2, N_new)

    elif flag == 2:
        N_new1 = np.shape(y)[axis]
        y = reflect_signal(y, axes=(axis,), flags=(flag,)) # is it necessary to reflect?
        N2 = np.shape(y)[axis]
        N_new = int(N2 / 2)+1 if to_slice else N2
        if N_new==N_new1:
            print('Can change code for flag ', flag)
        y = spectral_derivative(y, axis, n, L, N2, N_new)

    elif flag == 3:
        dim = np.ndim(y)
        for i in range(n):
            slicer1 = [slice(None)] * dim  # building slicer
            slicer2 = [slice(None)] * dim
            slicer1[axis] = slice(2, N)
            slicer2[axis] = slice(0, N - 2)
            y = (y[tuple(slicer1)] - y[tuple(slicer2)]) / (2 * dx)
            N = np.shape(y)[axis]
        # DD=np.concatenate(([np.nan]*n,DD))

    elif flag == 4:
        for i in range(n):
            y = np.gradient(y, axis=axis, edge_order=2) / dx

    else:
        print("NO FLAG GIVEN")
        return

    return y


def gauss_filtering(data, filter_widths, axes=(0,)):
    """
    using ndimage.gaussian_filter_ to filter the data,
    sigma array must be the same length as axes
    """
    if filter_widths==(0,):  # no filtering
        return data
    sigmas = np.array(filter_widths)/2
    dim = np.ndim(data)
    _sigmas = list(sigmas)
    _sigmas = [_sigmas.pop(0) if _ax in axes else 0 for _ax in range(dim)]
    return ndimage.gaussian_filter(data, sigma=_sigmas)


def butter_ba_filter(data, filter_width, dt, axis=0):
    if filter_width==0:  # no filtering
        return data
    order = 8
    highcut = 1/filter_width
    fs = 1/dt
    nyq = 0.5 * fs
    high = highcut / nyq
    print(high)
    b, a = butter(order, high)
    y = lfilter(b, a, data, axis=axis)
    return y


def butter_sos2_filter(data, filter_width, dt, axis=0, filter_order=3):
    """
    data: signal to filter
    filter_width: the width of the filter in units of dt
    dt: step in time of the signal
    """

    fs = 1 / dt
    if type(filter_width)==tuple:
        _btype='bandpass'
        f_cutoff = 1 / np.array(filter_width)
        print('Filter used: time=%.2g,%.2g, freq=%.2g,%.2g' % (*filter_width*dt, *f_cutoff))
    elif type(filter_width) == int:
        _btype='lowpass'
        f_cutoff = 1 / filter_width
        print('Filter used: time=%.2g, freq=%.2g' % (filter_width*dt, f_cutoff))
    else:
        print(type(filter_width))
        print(type(filter_width)==int)

    sos = butter(filter_order, f_cutoff, btype=_btype, output='sos', fs=fs)
    data_filt = sosfiltfilt(sos, data, axis=axis)
    return data_filt


def spectral_derivative(y, axis, n, L, N2, N_new):
    dk = np.pi / L  # maybe L=N2*dx
    k = np.fft.fftfreq(N2, 1 / N2) * dk
    k = k.reshape([-1 if _ax == axis else 1 for _ax in range(np.ndim(y))])
    y = np.fft.fft(y, axis=axis)
    y = ((1j * k) ** n) * y
    y = np.fft.ifft(y, axis=axis)
    if N2!=N_new:
        dim = np.ndim(y)
        slicer3 = [slice(None)] * dim
        slicer3[axis] = slice(0, N_new)
        y = y[tuple(slicer3)]
    y = np.real(y)
    return y


def poisson_eq(y, L, axes=(1,2), flags=(0,-1)):
    """
    gets an matrix y of physical length L_i, physical interval dx_i (and N_i points) for each dimension.
    the function returns the solution for the laplasian equation for 2 dim:
    vort = (d^2/dx^2 + d^2/dy^2)psi
    the function returns psi.
    flag is 1 if the function is even (derivative using cosine transform) ands
    -1 if it's odd (derivative using sine transform). If flag=2 the derivative
    is done using np.gradient. If flag=0 derivative is done for a
    periodic function using simple fft. Make sure in this case that the L you
    pass is that of the entire domain.
    """

    dim = np.ndim(y)
    freqs= [0]*dim
    slicer3 = [slice(None)] * dim

    for i in range(len(axes)):
        flag=flags[i]
        axis=axes[i]
        N=np.shape(y)[axis]

        if flag == 0:
            # dk = 2 * np.pi / L
            # k = dk * np.fft.fftfreq(N, 1 / N)
            # k = k.reshape([-1 if _ax == axis else 1 for _ax in range(dim)])

            y = reflect_signal(y, axes=(axis,), flags=(flag,))

            N2 = np.shape(y)[axis]
            dk = np.pi / L
            k = np.fft.fftfreq(N2, 1 / N2) * dk  # changed it to fftfreq. # only true when N2 is even
            k = k.reshape([-1 if _ax == axis else 1 for _ax in range(dim)])
            slicer3[axis] = slice(0, int(N2 / 2))

        elif flag == 1 or flag == -1:
            """
            is -1 only for sine? and it find cosine?
            """
            y = reflect_signal(y, axes=(axis,), flags=(flag,))

            N2 = np.shape(y)[axis]
            dk = np.pi / L
            k = np.fft.fftfreq(N2, 1 / N2) * dk  # changed it to fftfreq. # only true when N2 is even
            k = k.reshape([-1 if _ax == axis else 1 for _ax in range(dim)])

            slicer3[axis] = slice(0, int(N2 / 2) + 1)

        freqs[axis]=k

    yyt = np.fft.fftn(y, axes=axes)
    denom = -((freqs[axes[0]])**2+(freqs[axes[1]])**2)
    yyt = yyt / denom
    slicer4 = [slice(0,1) if _ax in axes else slice(None) for _ax in range(dim)]
    yyt[tuple(slicer4)] = 0 # yyt[:,0,0]=0 subtrackting the constant velocity
    yy = np.fft.ifftn(yyt, axes=axes)
    field = np.real(yy[tuple(slicer3)])

    return field
