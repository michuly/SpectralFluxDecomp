import numpy as np
import scipy as sp
from tools.data_analysis_tools import reflect_signal, gauss_filtering, butter_sos2_filter, COMPENSATING_FACTOR


def processing_and_fft(data, axes, minus_mean_axes=(), slice_axes=(), window_axes=(), reflect_axes=(), reflect_flags=(),
                       reflection_normalize=False, filter_flag=0, filter_widths=(), filter_axes=(), hours_cut=(),
                       data_type=32):
    """
    fft for data
    """
    # is minus axes the same as fft axes?
    data = data_processing(data, minus_mean_axes=minus_mean_axes, window_axes=window_axes, reflect_axes=reflect_axes,
                           reflect_flags=reflect_flags, reflection_normalize=reflection_normalize, filter_flag=filter_flag,
                           filter_widths=filter_widths, filter_axes=filter_axes, hours_cut=hours_cut)

    shape = np.shape(data)
    data = fft_and_normalizing(data, axes, data_type)

    if slice_axes:
        data = slice_and_multiply(data, axes=slice_axes)

    return data, shape


def fft_and_normalizing(data, axes, data_type):
    """
    fft for data
    """
    print('calculating fft... ')
    shape = np.shape(data)
    # factor to normalize data after fft: 1/N. np.sum(SPD)*dk=np.sum(np.abs(x)**2))/N
    # np.sum(fft(x)*conj(fft(x)))/N=np.sum(np.abs(x)**2))
    norm_factor = np.prod([shape[ax] for ax in axes])
    # norm_factor = np.prod([np.sqrt(shape[ax]) for ax in axes])
    # data_tf = np.fft.fftn(data, axes=axes) / norm_factor # NORMALIZING FACTOR 1/sqrt(N) for Variances
    if data_type==32:
        return (sp.fft.fftn(data, axes=axes, overwrite_x=True) / norm_factor).astype(np.complex64)
    elif data_type==64:
        return sp.fft.fftn(data, axes=axes, overwrite_x=True) / norm_factor


def freq_for_fft(N1, D1, N2=None, D2=None):
    """
    using np.fft to calculate the freq array of a given axis.
    it cuts half of the array, due to the symmetrical nature of the fft.
    """
    if N2 is None and D2 is None:
        freq = np.fft.fftfreq(N1, d=D1)[0:int(np.floor(N1 / 2) + 1)]  # cut the freq array in the middle
        if N1 % 2 == 0:
            freq[-1] = -freq[-1]  # if N is even, the last freq is negative, and need to be positive.
        return freq
    else:
        df_1, df_2 = 1 / (N1 * D1), 1 / (N2 * D2)
        # not sure if this should be mean, or mean of sqr. Roy code dell with same length.
        dk_h = np.sqrt(np.mean([df_1**2, df_2**2]))

        kx = np.abs(np.fft.fftfreq(N1, D1))
        ky = np.abs(np.fft.fftfreq(N2, D2))
        kx_mat, ky_mat = np.meshgrid(kx, ky, sparse=True, indexing='ij')
        kh_mat = np.hypot(kx_mat, ky_mat)  # r = sqrt(kx^2+ky^2)
        kh_array_normalized = (kh_mat/dk_h).round().astype(int).ravel()
        kh_array = np.unique(kh_array_normalized) * dk_h  # give kh right scale

        # NORMALIZING FACTOR, find factor to get density.
        # a_cor = 2 * np.pi * kh_array * dk_h
        # a_cor[0]=1
        # data_h /= a_cor # normalizing df_h
        return kh_array


def slice_and_multiply(data, axes=(1,)):
    """
    due to the symmetric nature of fft, we can slice the data to 2.
    after slicing we need to multiply the energy by 2, to maintain total energy.
    """
    print('Slicing... ')
    # slicing array to half
    shape = np.shape(data)
    dim = np.ndim(data)
    slicer = [slice(None)] * dim  # building slicer
    for ax in axes:  # only for the axes that will be cut.
        slicer[ax] = slice(0, int(np.floor(shape[ax] / 2) + 1))
    data = data[tuple(slicer)]

    # multiply the amplitudes (multiply velocity by sqrt(2) = multiply energy by 2)
    slicer = [slice(None)] * dim  # building slicer
    for ax in axes:  # only for the axes that will be cut.
        if np.mod(shape[ax], 2) == 0:  # even
            slicer[ax] = slice(1, int(np.floor(shape[ax] / 2)))
        else:  # odd
            slicer[ax] = slice(1, int(np.floor(shape[ax] / 2) + 1))
    data[tuple(slicer)] = data[tuple(slicer)] * np.sqrt(2) # NORMALIZING FACTOR, to keep energy.

    return data


def data_processing(data, minus_mean_axes=(), window_axes=(), reflect_axes=(), reflect_flags=(), filter_flag=0,
                    filter_widths=(), filter_axes=(), hours_cut=(), reflection_normalize=False, verbose=True):
    """
    data processing before fft: subtracting the mean value and applying windowing.


    according to the dimensions of the fft.
    """
    if reflect_axes:
        if verbose:
            print('reflecting signal...', reflect_axes)
        data = reflect_signal(data, axes=reflect_axes, flags=reflect_flags, normalize=reflection_normalize)

    if minus_mean_axes:
        if verbose:
            print('subtracting mean...', minus_mean_axes)
        data = data - np.mean(data, axis=minus_mean_axes, keepdims=True)

    if filter_flag!=0:
        if verbose:
            print('applying filter u...', filter_flag)
        if filter_flag==1:
            data = gauss_filtering(data, filter_widths=filter_widths, axes=filter_axes)
        elif filter_flag==2:
            data = butter_sos2_filter(data, filter_width=filter_widths[0], dt=1)

    if hours_cut:
        if verbose:
            print('cutting hours: ', hours_cut)
        dim = np.ndim(data)
        slicer = [slice(None)] * dim  # building slicer
        for ax in range(dim):  # only for the axes that will be cut.
            if hours_cut[ax]!=0:
                slicer[ax] = slice(hours_cut[ax], -hours_cut[ax])
        data = data[tuple(slicer)]

    if window_axes:
        if verbose:
            print('applying window...', window_axes)
        shape = np.shape(data)
        dim = np.ndim(data)
        window = 1
        for ax in range(dim):
            if ax in window_axes:
                window = window * np.hanning(shape[ax]).reshape([-1 if _ax == ax else 1 for _ax in range(dim)])
                # window = window * sp.signal.hanning(shape[ax]).reshape([-1 if _ax == ax else 1 for _ax in range(dim)])
        data = window * data * COMPENSATING_FACTOR ** len(window_axes) # NORMALIZING FACTOR compensating for windowing

        del window

    return data


def radial_profile(data, shape, dx_shape, axes=(1, 2), verbose=True):
    """
    this function calculates the psd for horizontal wavenumber (kh), from psd for cartesian wavenumbers (kx and ky).
    compared to Roy's code.
    kh_array construction is different, and fits different n's (and so different df's).
    for the same df's - the same as Roy's code.
    """
    if verbose:
        print("calculating radial profile")
    # calculating and normalizing kh
    n1, n2 = [shape[ax] for ax in axes]
    dx1, dx2 = [dx_shape[ax] for ax in axes]
    dim = np.ndim(data)
    df_1, df_2 = 1 / (n1 * dx1), 1 / (n2 * dx2)
    # not sure if this should be mean, or mean of sqr. Roy code dell with same length.
    dk_h = np.sqrt(np.mean([df_1**2, df_2**2]))

    kx = np.abs(np.fft.fftfreq(n1, dx1))
    ky = np.abs(np.fft.fftfreq(n2, dx2))
    kx_mat, ky_mat = np.meshgrid(kx, ky, sparse=True, indexing='ij')
    kh_mat = np.hypot(kx_mat, ky_mat)  # r = sqrt(kx^2+ky^2)

    # using kh matrix to get kh vector with unique values:
    kh_array_normalized = (kh_mat/dk_h).round().astype(int).ravel()
    kh_array = np.unique(kh_array_normalized) * dk_h  # give kh right scale

    # calculating psd
    if dim == 2:
        data = np.array([data])
    nt = np.shape(data)[0]
    data_h_i = np.zeros((nt, len(kh_array)))
    data_h_r = np.zeros((nt, len(kh_array)))
    for t_ind in range(nt):
        # data_h[i, :] = np.bincount(kh, weights=np.abs(data[i, :, :].ravel()))  # why absolute??
        # count = np.bincount(kh_array_normalized)
        data_h_i[t_ind, :] = np.bincount(kh_array_normalized, weights=np.imag(data[t_ind, :, :]).ravel())
        data_h_r[t_ind, :] = np.bincount(kh_array_normalized, weights=np.real(data[t_ind, :, :]).ravel())
    data_h = data_h_i + data_h_r
    del data_h_i, data_h_r
    if dim == 2:
        data_h = data_h[0]

    return kh_array, data_h