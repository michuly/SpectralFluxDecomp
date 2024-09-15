#!/usr/bin/env python
#################################################
# Load Modules
#################################################
from memory_profiler import profile
from open_data import *
import time
import sys
import os
# import R_tools_fort as tF
from scipy.fft import rfft2, irfft2
from scipy.ndimage import uniform_filter


def Forder(var):
    return np.asfortranarray(var.T, dtype=np.float32)
    # return var


def ncopen(filename, depth, opt='r'):
    return Dataset(filename + '.{0:04}'.format(depth) + '.nc', opt)


######################################
def makenan(x, val=np.nan):
    x[x == 0] = val
    x[x < -100000] = val
    x[x > 1000000] = val


######################################
@profile
def loadvd(data_sets, flag=0):
    vd = {}
    for var in ['u', 'v', 'w', 'dux', 'dvx', 'duy', 'dvy', 'duz', 'dvz']:
        vd[var] = np.squeeze(data_sets.return_data(var, flag))
        # print(vd[var].shape)
    return vd


######################################
Nx, Ny = 512, 513
dx = 390
kx = np.fft.fftfreq(Nx, d=1. / Nx)
ky = np.fft.fftfreq(Ny, d=1. / Ny)
nk = Nx // 2 + 1
KY, KX = np.meshgrid(ky[:nk], kx)

# @profile
def filt_spectral(ur):
    # n = lenth of wave in units (possible 2-512)
    print(ur.shape)
    print('CHECK FILTERING> IT DOES NOT MAKE SENSE TO ME)')
    uf = rfft2(ur, axes=(0,1))
    nt, nx, nk = uf.shape
    nc = nx / nfilt
    K2 = KX ** 2 + KY ** 2
    mask = K2 > nc ** 2
    uf[mask] = 0
    return irfft2(uf, axes=(0,1))


# @profile
def filt_spatial(ur):
    # n = lenth of wave in units (possible 2-512)
    print('Uniform filtering:', ur.shape)
    u_filt = uniform_filter(ur, (int(512 / nfilt), int(512 / nfilt), 0), None, 'mirror')
    return u_filt


def filtr(ur):
    if filter_type=='Spatial':
        return filt_spatial(ur)
    elif filter_type=='Spectral':
        return filt_spectral(ur)


######################################
# @profile
def prodf(u, v):
    return filtr(u * v) - filtr(u) * filtr(v)


######################################
def loadvar(var_str, _flag=''):
    if _flag=='w':
        flag=2
    elif _flag=='e':
        flag=1
    elif _flag=='':
        flag=0

    return np.squeeze(Forder(data_sets.return_data(var_str, flag)))


def tau(fields, flags=None):
    # print conv[dims[0]], conv[dims[1]]
    if flags is not None:
        flag1,flag2,flag3=flags

    def filt_t(data):
        print('time filtering:', data.shape)
        if flag3=='E':
            return butter_sos2_filter(data, filter_width=24, dt=1, axis=2)
        elif flag3 == 'W':
            return data - butter_sos2_filter(data, filter_width=16, dt=1, axis=2)

    if flags is None:
        return prodf(loadvar(fields[0]), loadvar(fields[1]))

    elif locality==1:
        return prodf(loadvar(fields[0], flag1), loadvar(fields[1], flag2))

    elif locality==2:
        u = loadvar(fields[0], flag1)
        v = loadvar(fields[1], flag2)
        return filt_t(filtr(u * v) - filtr(u) * filtr(v))

    elif locality==3:
        u = loadvar(fields[0], flag1)
        v = loadvar(fields[1], flag2)
        return filt_t(filtr(u * v)) - filt_t(filtr(u))*filt_t(filtr(v))


def rem(fields, flags):
    # print conv[dims[0]], conv[dims[1]]
    if flags is not None:
        flag1,flag2,flag3=flags
    def filt_t(data):
        print('time filtering:', data.shape)
        if flag3=='E':
            return butter_sos2_filter(data, filter_width=24, dt=1, axis=2)
        elif flag3 == 'W':
            return data - butter_sos2_filter(data, filter_width=16, dt=1, axis=2)

    if flags is None:
        return prodf(loadvar(fields[0]), loadvar(fields[1]))

    elif locality==1:
        return filtr(loadvar(fields[0], flag1)) * filtr(loadvar(fields[1], flag2))

    elif locality==2:
        u=loadvar(fields[0], flag1)
        v=loadvar(fields[1], flag2)
        return filt_t(filtr(u) * filtr(v))

    elif locality==3:
        u = loadvar(fields[0], flag1)
        v = loadvar(fields[1], flag2)
        return filt_t(filtr(u))*filt_t(filtr(v))



######################################
# @profile
def mapfilt(varlist):
    return [filtr(var) for var in varlist]
    # return [filt_spectral(var, n) for var in varlist]


def filters_to_var_name(direction, flags):
    var_str = 'SP'
    if direction == 'Vertical':
        var_str += 'V'
    elif direction == 'Horizontal':
        var_str += 'H'
    if flags == (0, 0, 0):
        return var_str + 'tot'
    elif flags == (3, 3, 3):
        return var_str
    for i in [2, 1]:
        if flags[i] == 1:
            var_str += 'e'
        elif flags[i] == 2:
            var_str += 'w'
    if flags[0] == 1:
        var_str += 'E'
    elif flags[0] == 2:
        var_str += 'W'
    return var_str


# w->div, e->rot
# @ray.remote
@profile
def tauijVij(nc_file_name, itime, ftime):
    # w->div, e->rot
    # calculate all flags togerther,, and save them differently.
    # calculate as is, just make saving different,s
    # how do they save it for different times??s
    # r = n * param.DX
    # print(r)

    def save_and_return(var_str, SP_tmp):
        print('Saving ', var_str)
        finalize_2D_coarse_graining_ray(nc_file_name, var_str, np.squeeze(Forder(SP_tmp)), itime, ftime)
        return SP_tmp

    ux, uy, vx, vy = mapfilt([loadvar('dux'), loadvar('duy'), loadvar('dvx'), loadvar('dvy')])

    print('SHAPE after filtering', ux.shape)
    # print 'Done smoothing'
    tauxx, tauyy, tauxy, tauyx = [tau('uu'), tau('vv'), tau('uv'), tau('vu')]
    save_and_return('SPH', -(tauxx * ux + tauxy * uy + tauyx * vx + tauyy * vy))

    del tauxx, ux, tauxy, uy, tauyx, vx, tauyy, vy

    if filter_type == 'Spatial':
        SPHtot = np.zeros((512, 513, ftime - itime))
    elif filter_type == 'Spectral':
        SPHtot = np.zeros((512, 512, ftime - itime))
    if SPHtot.shape[-1]==1:
        SPHtot=SPHtot[:,:,0]

    uxe, uye, vxe, vye = mapfilt([loadvar('dux','e'), loadvar('duy','e'), loadvar('dvx','e'), loadvar('dvy','e')])
    #
    tauxxee, tauyyee, tauxyee, tauyxee = [tau('uu', 'eeE'), tau('vv', 'eeE'), tau('uv', 'eeE'), tau('vu', 'eeE')]
    SPHtot += save_and_return('SPHeeE',  -(tauxxee * uxe + tauxyee * uye + tauyxee * vxe + tauyyee * vye))
    del tauxxee, tauyyee, tauxyee, tauyxee

    tauxxwe, tauyywe, tauxywe, tauyxwe = [tau('uu', 'weE'), tau('vv', 'weE'), tau('uv', 'weE'), tau('vu', 'weE')]
    SPHtot += save_and_return('SPHweE',  -(tauxxwe * uxe + tauxywe * uye + tauyxwe * vxe + tauyywe * vye))
    del tauxxwe, tauyywe, tauxywe, tauyxwe

    tauxxew, tauyyew, tauxyew, tauyxew = [tau('uu', 'ewE'), tau('vv', 'ewE'), tau('uv', 'ewE'), tau('vu', 'ewE')]
    SPHtot += save_and_return('SPHewE',  -(tauxxew * uxe + tauxyew * uye + tauyxew * vxe + tauyyew * vye))
    del tauxxew, tauyyew, tauxyew, tauyxew

    tauxxww, tauyyww, tauxyww, tauyxww = [tau('uu', 'wwE'), tau('vv', 'wwE'), tau('uv', 'wwE'), tau('vu', 'wwE')]
    SPHtot += save_and_return('SPHwwE',  -(tauxxww * uxe + tauxyww * uye + tauyxww * vxe + tauyyww * vye))
    del tauxxww, tauyyww, tauxyww, tauyxww
    del uxe, uye, vxe, vye
    #
    uxw, uyw, vxw, vyw = mapfilt([loadvar('dux','w'), loadvar('duy','w'), loadvar('dvx','w'), loadvar('dvy','w')])

    tauxxww, tauyyww, tauxyww, tauyxww = [tau('uu', 'wwW'), tau('vv', 'wwW'), tau('uv', 'wwW'), tau('vu', 'wwW')]
    SPHtot += save_and_return('SPHwwW',  -(tauxxww * uxw + tauxyww * uyw + tauyxww * vxw + tauyyww * vyw))
    del tauxxww, tauyyww, tauxyww, tauyxww

    tauxxee, tauyyee, tauxyee, tauyxee = [tau('uu', 'eeW'), tau('vv', 'eeW'), tau('uv', 'eeW'), tau('vu', 'eeW')]
    SPHtot += save_and_return('SPHeeW',  -(tauxxee * uxw + tauxyee * uyw + tauyxee * vxw + tauyyee * vyw))
    del tauxxee, tauyyee, tauxyee, tauyxee

    tauxxwe, tauyywe, tauxywe, tauyxwe = [tau('uu', 'weW'), tau('vv', 'weW'), tau('uv', 'weW'), tau('vu', 'weW')]
    SPHtot += save_and_return('SPHweW',  -(tauxxwe * uxw + tauxywe * uyw + tauyxwe * vxw + tauyywe * vyw))
    del tauxxwe, tauyywe, tauxywe, tauyxwe

    tauxxew, tauyyew, tauxyew, tauyxew = [tau('uu', 'ewW'), tau('vv', 'ewW'), tau('uv', 'ewW'), tau('vu', 'ewW')]
    SPHtot += save_and_return('SPHewW',  -(tauxxew * uxw + tauxyew * uyw + tauyxew * vxw + tauyyew * vyw))
    del tauxxew, tauyyew, tauxyew, tauyxew
    del uxw, uyw, vxw, vyw

    save_and_return('SPHtot', SPHtot)

    uz, vz = mapfilt([loadvar('duz'), loadvar('dvz')])
    taux, tauy = [tau('uw'), tau('vw')]
    save_and_return('SPV', -(taux * uz + tauy * vz))

    del uz, vz, taux, tauy

    uze, vze, uzw, vzw = mapfilt([loadvar('duz','e'), loadvar('dvz','e'), loadvar('duz','w'), loadvar('dvz','w')])

    if filter_type == 'Spatial':
        SPVtot = np.zeros((512, 513, ftime - itime))
    elif filter_type == 'Spectral':
        SPVtot = np.zeros((512, 512, ftime - itime))
    if SPVtot.shape[-1]==1:
        SPVtot=SPVtot[:,:,0]

    tauxee, tauyee = [tau('uw', 'eeE'), tau('vw', 'eeE')]
    SPVtot += save_and_return('SPVeeE', -(tauxee * uze + tauyee * vze))
    del tauyee, tauxee

    tauxew, tauyew = tau('uw', 'ewE'), tau('vw', 'ewE')
    SPVtot += save_and_return('SPVewE',  -(tauxew * uze + tauyew * vze))
    del tauxew, tauyew

    tauxwe, tauywe = tau('uw', 'weE'), tau('vw', 'weE')
    SPVtot += save_and_return('SPVweE',  -(tauxwe * uze + tauywe * vze))
    del tauxwe, tauywe

    tauxww, tauyww = [tau('uw', 'wwE'), tau('vw', 'wwE')]
    SPVtot += save_and_return('SPVwwE',  -(tauxww * uze + tauyww * vze))
    del tauxww, tauyww
    del uze, vze

    tauxww, tauyww = [tau('uw', 'wwW'), tau('vw', 'wwW')]
    SPVtot += save_and_return('SPVwwW',  -(tauxww * uzw + tauyww * vzw))
    del tauxww, tauyww

    tauxee, tauyee = [tau('uw', 'eeW'), tau('vw', 'eeW')]
    SPVtot += save_and_return('SPVeeW',  -(tauxee * uzw + tauyee * vzw))
    del tauyee, tauxee

    tauxew, tauyew = tau('uw', 'ewW'), tau('vw', 'ewW')
    SPVtot += save_and_return('SPVewW', -(tauxew * uzw + tauyew * vzw))
    del tauxew, tauyew

    tauxwe, tauywe = tau('uw', 'weW'), tau('vw', 'weW')
    SPVtot += save_and_return('SPVweW', -(tauxwe * uzw + tauywe * vzw))
    del tauxwe, tauywe
    del uzw, vzw

    save_and_return('SPVtot', SPVtot)


######################################
#############MAIN#####################
######################################


try:
    depth = int(sys.argv[1])
    exp = str(sys.argv[2])
    nfilt = int(sys.argv[3])

except:
    depth = 128
    nfilt = 10
    exp='Stochastic'
    # exp='Steady'
    # sys.exit(0)


# data processing
flux_remainder = 'CG_ray'  # 'vertical' or 'galilean'
directions = ['Horizontal', 'Vertical']
filter_width1 = 24
freq_domain1 = 'LF'
helm_domain1 = 'rot'
filter_width2 = 16
freq_domain2 = 'HF'  # HF
helm_domain2 = None  # divs
coor_shift = True
reduce_mean = False
Temporal = False
filter_type = 'Spatial'
locality=2

param = SimulationParameter(exp)

# Output
varlist = ['SPHeeE', 'SPHeeW', 'SPHwwW', 'SPHwwE', 'SPHweW', 'SPHewW', 'SPHewE', 'SPHweE', 'SPVeeE', 'SPVeeW',
           'SPVwwW', 'SPVwwE', 'SPVweW', 'SPVewW', 'SPVewE', 'SPVweE', 'SPH', 'SPHtot', 'SPV', 'SPVtot']


print('Process pid: ', os. getpid())
print('nfilt:', nfilt)
Nt = param.Nt
Nt0 = 0
NtEW = 722
print(NtEW)

######################################
print(os.uname()[1])
print('n filter:', nfilt * param.DX / 1000)

if exp == 'Stochastic' and depth == 77:
    pass
else:
    data_sets = TwoDataSets(filt_width1=filter_width1, freq_domain1=freq_domain1, helm_domain1=helm_domain1,
                            filt_width2=filter_width2, freq_domain2=freq_domain2, helm_domain2=helm_domain2,
                            coor_shift=coor_shift, exp=exp, depth=depth, reduce_mean=reduce_mean)

    start_time_ini = time.time()
    file_name = get_CG_file_name(depth, exp, filt1=filter_width1, filt2=filter_width2,
                                 domain1=freq_domain1, domain2=freq_domain2, helm1=helm_domain1,
                                 helm2=helm_domain2, coor_shift=coor_shift, reduce_mean=reduce_mean,
                                 temporal=Temporal, flux_remainder=flux_remainder)
    # file_name = file_name.replace('XY','CG_remainder_%d/XY_n%d' % (locality, nfilt))
    file_name = file_name.replace('XY','CG_local_%d/XY_n%d' % (locality, nfilt))

    if filter_type=='Spatial':
        shape=(722, 513, 512)
    elif filter_type=='Spectral':
        shape = (722, 512, 512)

    initialize_2D_coarse_graining_ray(file_name, shape)

    # time arrays options:
    # time_array = np.linspace(Nt0,NtEW,10).astype(int)
    time_array = [0,722]
    # time_array = np.arange(NtEW)

    start_time = time.time()
    for itime, ftime in zip(time_array[:-1],time_array[1:]):
        start_time_step = time.time()
        print(f'itime: {itime}, {ftime}')
        data_sets.data_set1.ind_min = itime
        data_sets.data_set1.ind_max = ftime
        data_sets.data_set2.ind_min = itime
        data_sets.data_set2.ind_max = ftime

        start_time_step = time.time()
        start_time_load = time.time()
        print()
        print('TIME for load %.1fsec' % (time.time() - start_time_load))
        print()

        start_time_calc = time.time()
        # remainderijVij(file_name, itime, ftime)
        tauijVij(file_name, itime, ftime)

        print()
        print('TIME for Calculation %.1fsec' % (time.time() - start_time_calc))
        print()

        start_time_save = time.time()
        print()
        print('TIME for saving: %.1fsec' % (time.time() - start_time_save))
        print()

        print()
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        print('TIME for step: %.1fsec' % (time.time() - start_time_step))
        print()

        print()
        print('TIME for initializing: %.1fsec' % (time.time() - start_time_ini))
        print()
