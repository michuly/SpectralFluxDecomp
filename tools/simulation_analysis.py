import sys

sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
import numpy as np
from tools.PSD_calculations import psd_1d, psd_2d, psd_3d, psd_enstrophy_1d
from tools.data_analysis_tools import derivative, poisson_eq
from tools.flow_analysis_tools import processing_and_fft, data_processing, freq_for_fft
from tools.flow_analysis_tools2 import *

F_cor = 1.2e-4
N = 10e-2
rho_0 = 1e3


def Horizontal_SF_netcdf(data_sets, flags, vel_str, axes=(1,2), ref_flag=None):
    """
    Pi(k)=Int{ Real[ F{u}_k' F{u_m du_dx_m}_k' ] }dk'
    """
    print('calc SF...')
    if axes==(1,2):
        reflect_axes=(1, 2)
        reflect_flags=(ref_flag, 0)
        window_axes = ()
        adv_u = data_sets.return_data('u', flags[1]) * data_sets.return_data('d%sx' % vel_str, flags[2]) + \
                data_sets.return_data('v', flags[1]) * data_sets.return_data('d%sy' % vel_str, flags[2])
        adv_u = data_processing(adv_u, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        adv_u, shape = processing_and_fft(adv_u, axes=axes, data_type=32)

    elif axes==(0,):
        reflect_axes=()
        reflect_flags=()
        window_axes = (0, )
        u = data_processing(data_sets.return_data('u', flags[1]), window_axes=window_axes, reflect_axes=reflect_axes,
                            reflect_flags=reflect_flags)
        dvx = data_processing(data_sets.return_data('d%sx' % vel_str, flags[2]), window_axes=window_axes,
                              reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        v = data_processing(data_sets.return_data('v', flags[1]), window_axes=window_axes, reflect_axes=reflect_axes,
                            reflect_flags=reflect_flags)
        duy = data_processing(data_sets.return_data('d%sy' % vel_str, flags[2]), window_axes=window_axes,
                              reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        adv_u, shape = processing_and_fft(u*dvx+v*duy, axes=axes, data_type=32)

    u = data_sets.return_data(vel_str, flags[0])
    u = data_processing(u, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    u, shape = processing_and_fft(u, axes=axes, data_type=32)  # flags=(1, 0)

    SF = np.real(u * np.conj(adv_u))

    del adv_u

    return shape, SF, u


def Vertical_SF_netcdf(data_sets, flags, vel_str, axes, ref_flag, u=None):
    """
    Pi(k)=Int{ Real[ F{u}_k' F{u_m du_dx_m}_k' ] }dk'
    """
    print('calc SF...')
    if axes==(1,2):
        reflect_axes = (1, 2)
        reflect_flags = (ref_flag, 0)
        window_axes = ()
        adv_u = data_sets.return_data('w', flags[1]) * data_sets.return_data('d%sz' % vel_str, flags[2])
        adv_u = data_processing(adv_u, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        adv_u, shape = processing_and_fft(adv_u, axes=axes, data_type=32)
    elif axes == (0,):
        reflect_axes = ()
        reflect_flags = ()
        window_axes = (0,)
        w = data_processing(data_sets.return_data('w', flags[1]), window_axes=window_axes, reflect_axes=reflect_axes,
                            reflect_flags=reflect_flags)
        duz = data_processing(data_sets.return_data('d%sz' % vel_str, flags[2]), window_axes=window_axes,
                              reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        adv_u, shape = processing_and_fft(w*duz, axes=axes, data_type=32)

    if u is None:
        u = data_sets.return_data(vel_str, flags[0])
        u = data_processing(u, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
        u, shape = processing_and_fft(u, axes=axes, data_type=32)

    SF = np.real(u * np.conj(adv_u))
    del u, adv_u

    return shape, SF


def Enstrophy_flux_netcdf(data_sets, flags, axes=(1,2)):
    """
    Pi(k)=Int{ Real[ F{u}_k' F{u_m du_dx_m}_k' ] }dk'
    """
    print('calc SF...')
    reflect_axes=(1, 2)
    window_axes = ()

    vort = data_sets.return_data('dvx', flags[2]) - data_sets.return_data('duy', flags[2])
    reflect_flags = (-1, 0)

    dvort=derivative(vort, L=data_sets.L, dx=data_sets.DX, n=1, axis=data_sets.AX_X, flag=0)
    adv_vort = data_sets.return_data('u', flags[1]) * dvort
    dvort=derivative(vort, L=data_sets.L, dx=data_sets.DX, n=1, axis=data_sets.AX_Y, flag=-1)
    adv_vort += data_sets.return_data('v', flags[1]) * dvort
    del dvort
    adv_vort = data_processing(adv_vort, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    adv_vort, shape = processing_and_fft(adv_vort, axes=axes, data_type=32)

    if flags[0]!=flags[2]:
        vort = data_sets.return_data('dvx', flags[0]) - data_sets.return_data('duy', flags[0])

    vort = data_processing(vort, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    vort, shape = processing_and_fft(vort, axes=axes, data_type=32)  # flags=(1, 0)

    SF_hor = np.real(vort * np.conj(adv_vort))

    adv_vort = derivative(data_sets.return_data('dvz', flags[2]), L=data_sets.L, dx=data_sets.DX, n=1,
                          axis=data_sets.AX_X, flag=0)
    adv_vort -= derivative(data_sets.return_data('duz', flags[2]), L=data_sets.L, dx=data_sets.DX, n=1,
                           axis=data_sets.AX_Y, flag=1)
    adv_vort *= data_sets.return_data('w', flags[1])
    adv_vort = data_processing(adv_vort, window_axes=window_axes, reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    adv_vort, shape = processing_and_fft(adv_vort, axes=axes, data_type=32)

    SF_ver = np.real(vort * np.conj(adv_vort))

    return shape, SF_hor, SF_ver


def psd_1d_netcdf(data_set, dx_shape, axis=0, sub_mean=False, window=False, hours_cut=0):
    """
    this functions takes a velocity matrix and calculates its energy spectrum Vs. frequency (time).
    """
    _, psd1 = psd_1d(data_set.return_data('u'), dx_shape=dx_shape, axis=axis, sub_mean=sub_mean, hours_cut=hours_cut,
                     window=window)
    freq, psd2 = psd_1d(data_set.return_data('v'), dx_shape=dx_shape, axis=axis, sub_mean=sub_mean, hours_cut=hours_cut,
                        window=window)
    return freq, psd1 + psd2


def psd_2d_netcdf(data, dx_shape, axes=(1, 2), sub_mean=True, reflecting=True, window=False):
    """
    this functions takes a velocity matrix and calculates its energy spectrum Vs. horizontal wavenumber (space).
    using 2d fft for each time step, then averaging over them. (gives the same result as 'psd_horizontal').
    """
    reflect_axes = ()
    reflect_flags = ()

    if reflecting:
        reflect_axes = (1, 2)
        reflect_flags = (1, 0)
    _, psd1 = psd_2d(data.return_data('u'), dx_shape=dx_shape, axes=axes, sub_mean=sub_mean, window=window,
                     reflect_axes=reflect_axes, reflect_flags=reflect_flags)

    if reflecting:
        reflect_axes = (1, 2)
        reflect_flags = (-1, 0)
    freq, psd2 = psd_2d(data.return_data('v'), dx_shape=dx_shape, axes=axes, sub_mean=sub_mean, window=window,
                        reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    return freq, psd1 + psd2


def psd_2d_xt_netcdf(data, dx_shape, axes=(0, 2), sub_mean=False, window=False):
    """
    this functions takes a velocity matrix and calculates its energy spectrum Vs. horizontal wavenumber (space).
    using 2d fft for each time step, then averaging over them. (gives the same result as 'psd_horizontal').
    """
    reflect_axes = ()
    reflect_flags = ()

    _, psd1 = psd_2d(data.return_data('u'), dx_shape=dx_shape, axes=axes, sub_mean=sub_mean, window=window,
                     reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    freq, psd2 = psd_2d(data.return_data('v'), dx_shape=dx_shape, axes=axes, sub_mean=sub_mean, window=window,
                        reflect_axes=reflect_axes, reflect_flags=reflect_flags)
    return freq, psd1 + psd2


def psd_3d_netcdf(data, dx_shape, sub_mean=True, reflecting=True, window=False, hours_cut=0):
    """
    this functions takes a velocity matrix and calculates its energy spectrum Vs. horizontal wavenumber (space).
    using 2d fft for each time step, then averaging over them. (gives the same result as 'psd_horizontal').
    """

    _, _, psd1 = psd_3d(data.return_data('u'), dx_shape=dx_shape, sub_mean=sub_mean, window=window,
                        reflecting=reflecting, hours_cut=hours_cut)
    freq_t, freq_h, psd2 = psd_3d(data.return_data('v'), dx_shape=dx_shape, sub_mean=sub_mean, window=window,
                                  reflecting=reflecting, hours_cut=hours_cut)
    return freq_t, freq_h, psd1 + psd2


def helmholz_decomposition_netcdf(data_set, use_derivatives=True, vertical_derivative=False):
    """
    check for 2dim!
    """
    vort = vorticity_netcdf(data_set, use_derivatives=use_derivatives, vertical_derivative=vertical_derivative)
    return helmholz_decomposition_vorticity(vort, data_set.L, data_set.DX, ax_x=data_set.AX_X, ax_y=data_set.AX_Y)


def dbz_decomposition_netcdf(data_set):
    """
    NEED explanation
    """
    dbz = data_set.return_data('dbz')
    dbz[np.where(dbz < F_cor ** 2)] = F_cor ** 2
    N2 = data_set.return_data('N2')  # this is essentially dbz, need to build a different data for N2
    dbz_N2 = dbz / N2 - 1
    dbz_N2 *= F_cor
    psi = poisson_eq(dbz_N2, L=data_set.L, axes=(data_set.AX_X, data_set.AX_Y), flags=(0, 1))
    u_dbz = -derivative(y=psi, L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y, flag=1)
    v_dbz = derivative(y=psi, L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X, flag=0)

    return u_dbz, v_dbz


def vorticity_netcdf(data_set, use_derivatives=True, vertical_derivative=False, check_existing=True):
    """
    vort = dv/dx - du/dy
    dx has to be the same for x and y.
    """
    if vertical_derivative:
        str_format = 'd%sz'
    else:
        str_format = '%s'

    if check_existing and data_set.check_variable_exist(var_str='Vort'):
        return data_set.return_data('Vort', float_type=64)
    elif not vertical_derivative and use_derivatives:
        return data_set.return_data('dvx', float_type=64) - data_set.return_data('duy', float_type=64)
    else:
        return vorticity(data_set.return_data(str_format % 'u', float_type=64),
                         data_set.return_data(str_format % 'v', float_type=64), data_set.L, data_set.DX, data_set.AX_X,
                         data_set.AX_Y)


def divergence_netcdf(data_set, use_derivatives=True):
    """
    div = du/dx + dv/dy
    dx has to be the same for x and y.
    """
    if use_derivatives:
        return data_set.return_data('dux', float_type=64) + data_set.return_data('dvy', float_type=64)
    else:
        return divergence(data_set.return_data('u', float_type=64), data_set.return_data('v', float_type=64),
                          data_set.L, data_set.DX, data_set.AX_X, data_set.AX_Y)


def rossby_number_netcdf(data_set, use_derivatives=True):
    """
    rossby = vort/coriolis_frequency
    """
    return vorticity_netcdf(data_set, use_derivatives) / F_cor


def Richardson_number_netcdf(data_set):
    """
    Richardson number = N^2 / duz
    """
    helm = data_set.helm_domain
    data_set.helm_domain = None
    N2 = np.mean(data_set.return_data('dbz'), axis=0, keepdims=True)
    data_set.helm_domain = helm
    return Richardson_number(duz=data_set.return_data('duz'), dvz=data_set.return_data('dvz'), N2=N2)


def Sigma_netcdf(data_set):
    """
    Based on Baldawa 2021
    """
    return Sigma_calc(data_set.return_data('dux'), data_set.return_data('duy'), data_set.return_data('dvx'),
                      data_set.return_data('dvy'))


def Ertel_PV_ver_netcdf(data_set):
    """
    rossby = vort/coriolis_frequency
    """
    vort = vorticity_netcdf(data_set)
    filter_width = data_set.filter_width
    freq_domain = data_set.freq_domain
    helm_domain = data_set.helm_domain
    geostrophic = data_set.geostrophic
    data_set.filter_width, data_set.freq_domain, data_set.helm_domain = 0, None, None
    data_set.geostrophic = False
    N2 = np.mean(data_set.return_data('dbz'), axis=0, keepdims=True)
    data_set.filter_width, data_set.freq_domain, data_set.helm_domain = filter_width, freq_domain, helm_domain
    data_set.geostrophic = geostrophic
    N2[N2 < 0] = np.NaN
    return Ertel_PV_ver(vort, dbz=N2)  # make sure with Roy that this is N2 and not dbz


def Ertel_PV_hor_netcdf(data_set, TWB=False):
    """
    rossby = vort/coriolis_frequency
    """
    if TWB:
        return TWB_Ertel_PV_hor(data_set.return_data('duz'), data_set.return_data('dvz'))
    else:
        helm = data_set.helm_domain
        geostrophic = data_set.geostrophic
        data_set.helm_domain = None
        data_set.geostrophic = False
        dbx = data_set.return_data('dbx')
        dby = data_set.return_data('dby')
        data_set.helm_domain = helm
        data_set.geostrophic = geostrophic
        return Ertel_PV_hor(data_set.return_data('duz'), data_set.return_data('dvz'), dbx, dby)


def Curvature_netcdf(data_set):
    filter_width = data_set.filter_width
    freq_domain = data_set.freq_domain
    helm_domain = data_set.helm_domain
    geostrophic = data_set.geostrophic
    u = data_set.return_data('u')
    v = data_set.return_data('v')
    duy = data_set.return_data('duy')
    dvx = data_set.return_data('dvx')
    # calculating V as the geostrophic velocity, using the vorticity to get the sign.
    if filter_width == 24 and freq_domain == 'LF' and helm_domain == 'rot' and not geostrophic:
        V = np.sqrt(u ** 2 + v ** 2)
        sign = np.sign(dvx - duy)
    else:
        data_set.filter_width, data_set.freq_domain, data_set.helm_domain = 24, 'LF', 'rot'
        data_set.geostrophic = False
        V = np.sqrt(data_set.return_data('u') ** 2 + data_set.return_data('v') ** 2)
        sign = np.sign(data_set.return_data('dvx') - data_set.return_data('duy'))
        data_set.filter_width, data_set.freq_domain, data_set.helm_domain = filter_width, freq_domain, helm_domain
        data_set.geostrophic = geostrophic
    V = sign * V

    return Curvature(u, v, data_set.return_data('dux'), duy, dvx, data_set.return_data('dvy'), V)


def Curvature_PV_hor_netcdf(data_set):
    """
    based on Baladawa 2021
    Curvature PV
    V geostrophic velocity
    """
    # Sigma = Sigma_calc(dux, duy, dvx, dvy)
    return Curvature_PV_hor(data_set)


def Curvature_PV_ver_netcdf(data_set):
    """
    based on Baladawa 2021
    Curvature PV
    V geostrophic velocity
    """
    # Sigma = Sigma_calc(dux, duy, dvx, dvy)
    return Curvature_PV_ver(data_set)


def get_coordinate_shift(exp):
    if exp == 'Stochastic':
        U = 2.289  # unit dx (390) to dt (3600)
    elif exp == 'Steady':
        U = 2.43
    return U

