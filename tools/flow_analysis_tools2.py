from tools.data_analysis_tools import derivative, poisson_eq
from tools.flow_analysis_tools import *
from tools.flow_analysis_tools import processing_and_fft
from tools.PSD_calculations import psd_1d

F_cor = 1.2e-4
N = 10e-2


def vorticity(data_u, data_v, L, dx, ax_x=1, ax_y=0):
    """
    vort = dv/dx - du/dy
    dx has to be the same for x and y.
    check for 2dim!
    """
    return derivative(data_v, L=L, dx=dx, n=1, axis=ax_x, flag=0) - \
        derivative(data_u, L=L, dx=dx, n=1, axis=ax_y, flag=1)


def rossby_number(data_u, data_v, dx, ax_x=1, ax_y=0):
    """
    rossby = vort/coriolis_frequency
    """
    return vorticity(data_u, data_v, dx, ax_x, ax_y) / F_cor


def divergence(data_u, data_v, L, dx, ax_x=1, ax_y=0):
    """
    div = du/dx + dv/dy
    dx has to be the same for x and y.
    check for 2dim!
    """
    return derivative(data_u, L=L, dx=dx, n=1, axis=ax_x, flag=0) + \
        derivative(data_v, L=L, dx=dx, n=1, axis=ax_y, flag=-1)


def helmholz_decomposition_vorticity(vort, L, dx, ax_x=1, ax_y=0):
    """
    check for 2dim!
    """
    psi = poisson_eq(vort, L=L, axes=(ax_x, ax_y), flags=(0, -1))
    u_rot = -derivative(y=psi, L=L, dx=dx, n=1, axis=ax_y, flag=-1)
    v_rot = derivative(y=psi, L=L, dx=dx, n=1, axis=ax_x, flag=0)

    return u_rot, v_rot


def helmholz_decomposition(u, v, L, dx, ax_x=1, ax_y=0):
    """
    check for 2dim!
    """
    return helmholz_decomposition_vorticity(vorticity(u, v, L, dx, ax_x, ax_y), L, dx, ax_x, ax_y)


def psd_1d_uv(u, v, dx_shape, axis=0, sub_mean=True, hours_cut=0):
    """
    this functions takes a velocity matrix and calculates its energy spectrum Vs. frequency (time).
    """
    _, psd1 = psd_1d(u, dx_shape=dx_shape, axis=axis, sub_mean=sub_mean, hours_cut=hours_cut)
    freq, psd2 = psd_1d(v, dx_shape=dx_shape, axis=axis, sub_mean=sub_mean, hours_cut=hours_cut)
    return freq, psd1 + psd2


def KE_TOT(u, v):
    """
    0.5*u_i*u_i
    """
    return 0.5 * (u ** 2 + v ** 2)


def KE_TOT_LF(u, v, sigma, ax_t):
    """
    0.5*tilde{u_i*u_i}
    """
    u2_tild = gauss_filtering(u ** 2, sigmas=(sigma,), axes=(ax_t,))
    v2_tild = gauss_filtering(v ** 2, sigmas=(sigma,), axes=(ax_t,))
    return 0.5 * (u2_tild + v2_tild)


def KE_LF(u, v, sigma, ax_t):
    """
    0.5*tilde{u_i}*tilde{u_i}
    """
    u_tild = gauss_filtering(u, sigmas=(sigma,), axes=(ax_t,))
    v_tild = gauss_filtering(v, sigmas=(sigma,), axes=(ax_t,))
    return 0.5 * (u_tild ** 2 + v_tild ** 2)


def KE_HF(u, v, sigma, ax_t, KE_tot_lf=None, KE_lf=None):
    """
    0.5*tilde{u_i*u_i}-KE_LF
    """
    if KE_tot_lf is None:
        KE_tot_lf = KE_TOT_LF(u, v, sigma, ax_t)
    if KE_lf is None:
        KE_lf = KE_LF(u, v, sigma, ax_t)
    return KE_tot_lf - KE_lf


def APE_TOT(b):
    """
    0.5*b*b
    """
    return 0.5 / N ** 2 * b ** 2


def APE_TOT_LF(b, sigma, ax_t):
    """
    0.5*tilde{b*b}
    """
    b2_tild = gauss_filtering(b ** 2, sigmas=(sigma,), axes=(ax_t,))
    return 0.5 / N ** 2 * b2_tild


def APE_LF(b, sigma, ax_t):
    """
    0.5*tilde{b}*tilde{b}
    """
    b_tild = gauss_filtering(b, sigmas=(sigma,), axes=(ax_t,))
    return 0.5 / N ** 2 * b_tild ** 2


def APE_HF(b, sigma, ax_t, APE_tot_lf=None, APE_lf=None):
    """
    APE_TOT_LF-APE_LF
    """
    if APE_tot_lf is None:
        APE_tot_lf = APE_TOT_LF(b, sigma, ax_t)
    if APE_lf is None:
        APE_lf = APE_LF(b, sigma, ax_t)
    return APE_tot_lf - APE_lf


def Richardson_number(duz, dvz, N2):
    """
    Richardson number = N^2 / duz^2
    """
    N2[np.where(N2 < 0)] = np.NaN  # OR =0
    Ri = N2 / (duz ** 2 + dvz ** 2)
    return Ri


def Sigma_calc(dux, duy, dvx, dvy):
    """
    based on Baladawa 2021
    """
    Sigma_n = dux - dvy
    Sigma_s = dvx + duy
    Sigma = np.sqrt(Sigma_n ** 2 + Sigma_s ** 2)
    return Sigma


def Curvature(u, v, dux, duy, dvx, dvy, V):
    """
    based on Buckingham 2021
    """
    R = (u ** 2 + v ** 2) ** (3 / 2) / (u ** 2 * dvx - v ** 2 * duy + u * v * (dvy - dux))
    return 2 * V / R / F_cor


def Baldawa_Instability_Criteria(u, v, dux, dvx, duy, dvy, duz, dvz, Vort, N2, V=None):
    """
    based on Baladawa 2021
    Curvature PV
    V geostrophic velocity
    """
    Ro = Vort / F_cor
    Sigma = Sigma_calc(dux, duy, dvx, dvy)
    Cu = Curvature(u, v, dux, duy, dvx, dvy, V)
    Ri = Richardson_number(duz, dvz, N2)
    Phi = (1 + Cu) * (1 + Ro) - (1 + Cu) ** 2 * Ri - 1
    return Ro, Sigma / F_cor, Phi


def Normalized_Ertel_PV_ver(vort):
    """
    fq = (f + dv/dx - du/dy) * db/dz
    """
    return 1 + vort / F_cor


def Ertel_PV_ver(vort, dbz):
    """
    fq = (f + dv/dx - du/dy) * db/dz
    """
    return F_cor * (F_cor + vort) * dbz


def Curvature_PV_ver(data_set):
    """
    based on Baladawa 2021
    Barotropic part of Curvature PV
    """
    Cu = data_set.return_data('Cu')
    Ro = data_set.return_data('Vort') / F_cor

    return (1 + Cu) * (1 + Ro)


def Normalized_Ertel_PV_hor(duz, dvz, dbz):
    """
    fq = (du/dz - dw/dx) * db/dy + (dw/dy - dv/dz) * db/dx
    """
    return -(duz ** 2 + dvz ** 2) / dbz


def Ertel_PV_hor(duz, dvz, dbx, dby):
    """
    fq = (du/dz - dw/dx) * db/dy + (dw/dy - dv/dz) * db/dx
    """
    return F_cor * (duz * dby - dvz * dbx)


def TWB_Ertel_PV_hor(duz, dvz):
    """
    fq = f [ (du/dz - dw/dx) * db/dy + (dw/dy - dv/dz) * db/dx ]
    """
    return -(F_cor ** 2) * (duz ** 2 + dvz ** 2)


def Curvature_PV_hor(data_set):
    """
    based on Baladawa 2021
    Baroclinic oart of Curvature PV
    """
    Cu = data_set.return_data('Cu')
    Ri = data_set.return_data('Ri')
    return (1 + Cu) ** 2 / Ri


def Instability_Criteria(Q_vert, Q_bc, N2, Ro=None):
    """
    building a matrix of instabilty
    No Instability = NaN  [ (Q_vert + Q_bc) > 0 ]

    Cyclonic > 0  [ Ro > 0 ]
    AntiCyclonic < 0  [ Ro <= 0 ]

    |Gravitational Instability| = 1 [ (Q_vert + Q_bc) < 0 && N2 < 0 && Q_vert < 0 ]
    |Gravitational/Symmetric Instability| = 2 [ N2 > 0 && Q_vert > 0 && Q_bc < 0 && |Q_bc| > Q_vert ] OR
                                [ (Q_vert + Q_bc) < 0 &&  N2 > 0 && Q_vert > 0]

    |Inertial Instability| = 3  [ (Q_vert + Q_bc) < 0 && N2 > 0 && Q_vert < 0 ]
    |Symmetric Instability| = 4  [ N2 > 0 && Q_vert > 0 && Q_bc < 0 && |Q_bc| > Q_vert ] OR
                                [ (Q_vert + Q_bc) < 0 &&  N2 > 0 && Q_vert > 0]

    check that [ (Q_vert + Q_bc) < 0 && N2 < 0 ] == |Gravitational Instability| or |Gravitational/Symmetric Instability|
    check that [ (Q_vert + Q_bc) < 0 && N2 > 0 ] == |Inertial Instability| or |Symmetric Instability|


    *** We assume f > 0 !! ***
    """
    Instability = np.full(np.shape(Q_vert), np.NaN)

    Crit_Qtot = np.full(np.shape(Q_vert), False)
    Crit_Qvert = np.full(np.shape(Q_vert), False)
    Crit_N2_neg = np.full(np.shape(N2), False)
    Crit_N2_pos = np.full(np.shape(N2), False)

    Crit_Qtot[(Q_vert + Q_bc) < 0] = True
    Crit_Qvert[Q_vert < 0] = True
    Crit_N2_neg[N2 <= 0] = True
    Crit_N2_pos[N2 > 0] = True

    Instability[np.logical_and.reduce(np.array((Crit_Qtot, Crit_N2_neg, Crit_Qvert)))] = 1
    Instability[np.logical_and.reduce(np.array((Crit_Qtot, Crit_N2_neg, np.logical_not(Crit_Qvert))))] = 2
    Instability[np.logical_and.reduce(np.array((Crit_Qtot, Crit_N2_pos, Crit_Qvert)))] = 3
    Instability[np.logical_and.reduce(np.array((Crit_Qtot, Crit_N2_pos, np.logical_not(Crit_Qvert))))] = 4

    if Ro is not None:
        Instability[Ro < 0] = - Instability[Ro < 0]

    return Instability


def rms(y):
    return np.sqrt(np.mean(y ** 2))
