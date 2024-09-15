#!/usr/bin/env python
#################################################
# Load Modules
#################################################
from open_data import *
import time
import sys
import os
# import R_tools_fort as tF
import ray
from scipy.fft import rfft2, irfft2
from pathlib import Path
import shutil
import logging


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
def loadvd(data_sets, flag=0):
    vd = {}
    for var in ['u', 'v', 'w', 'dux', 'dvx', 'duy', 'dvy', 'duz', 'dvz']:
        vd[var] = np.squeeze(Forder(data_sets.return_data(var, flag)))
        # print(vd[var].shape)
    return vd


######################################
Nx, Ny = 512, 513
dx = 390
kx = np.fft.fftfreq(Nx, d=1. / Nx)
ky = np.fft.fftfreq(Ny, d=1. / Ny)
nk = Nx // 2 + 1
KY, KX = np.meshgrid(ky[:nk], kx)


def filtr(ur, n):
    # n = lenth of wave in units (possible 2-512)
    print(ur.shape)
    uf = rfft2(ur)
    if uf.ndim==2:
        nx, nk = uf.shape
        nc = nx / n
        K2 = KX ** 2 + KY ** 2
        mask = K2 > nc ** 2
        uf[mask] = 0
        return irfft2(uf)
    elif uf.ndim==3:
        nt, nx, nk = uf.shape
        nc = nx / n
        K2 = KX ** 2 + KY ** 2
        mask = K2 > nc ** 2
        uf[np.repeat(mask[np.newaxis,:,:],nt,0)] = 0
        return irfft2(uf, axis=(1,2))


######################################
def prodf(u, v, n):
    return filtr(u * v, n) - filtr(u, n) * filtr(v, n)


######################################
def tau(dims, vd, n, vd2=None):
    # print conv[dims[0]], conv[dims[1]]
    if vd2 is None:
        return prodf(vd[dims[0]], vd[dims[1]], n)
    else:
        return [prodf(vd[dims[0]], vd2[dims[1]], n), prodf(vd2[dims[0]], vd[dims[1]], n)]


######################################
def tau2(dims, vd, n, vd2=None):
    # print conv[dims[0]], conv[dims[1]]
    if vd2 is None:
        return prodf(vd[dims[0]], vd[dims[1]], n)
    else:
        return prodf(vd[dims[0]], vd2[dims[1]], n)


######################################

def mapfilt(varlist, n):
    return [filtr(var, n) for var in varlist]


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


def redirecting_logs():
    # Ray stores logs in `/tmp/ray/session_latest/logs` by default.
    ray_log_dir = Path("/tmp/ray/session_latest/")
    if ray_log_dir.exists():
        for log_file in ray_log_dir.iterdir():
            if log_file.is_file():
                destination_file = Path(destination_log_dir) / log_file.name

                if destination_file.exists():
                    new_destination_file = destination_file
                    counter = 1
                    while new_destination_file.exists():
                        new_destination_file = destination_file.with_name(
                            f"{log_file.stem}_{counter}{log_file.suffix}")
                        counter += 1
                    shutil.move(str(log_file), new_destination_file)
                else:
                    shutil.move(str(log_file), destination_file)
        print(f"All log files moved to {destination_log_dir}")
    else:
        print(f"Ray log directory does not exist: {ray_log_dir}")


###########################################
# Roy modified
# Narray = list(range(1,12))
# Narray.extend(list(range(12,int(110/scalefac),int(6/scalefac))))
# Narray = np.array(Narray)[0:-1:2]

# def nanmeanzones(varval, Nbands):

# w->div, e->rot
@ray.remote
def tauijVij(vd, vdw, vde, param, varlist, n):
    # w->div, e->rot
    # calculate all flags togerther,, and save them differently.
    # calculate as is, just make saving different,s
    # how do they save it for different times??s
    od = {}
    r = n * param.DX
    print(r)

    ux, uy, vx, vy = mapfilt([vd['dux'], vd['duy'], vd['dvx'], vd['dvy']], n)
    uxe, uye, vxe, vye = mapfilt([vde['dux'], vde['duy'], vde['dvx'], vde['dvy']], n)
    uxw, uyw, vxw, vyw = mapfilt([vdw['dux'], vdw['duy'], vdw['dvx'], vdw['dvy']], n)

    tauxx, tauyy, tauxy, tauyx = [tau('uu', vd, n), tau('vv', vd, n), tau('uv', vd, n), tau('vu', vd, n)]
    tauxxww, tauyyww, tauxyww, tauyxww = [tau('uu', vdw, n), tau('vv', vdw, n), tau('uv', vdw, n), tau('vu', vdw, n)]
    tauxxee, tauyyee, tauxyee, tauyxee = [tau('uu', vde, n), tau('vv', vde, n), tau('uv', vde, n), tau('vu', vde, n)]
    tauxxew, tauyyew, tauxyew, tauyxew = [tau2('uu', vde, n, vdw), tau2('vv', vde, n, vdw), tau2('uv', vde, n, vdw),
                                          tau2('vu', vde, n, vdw)]
    tauxxwe, tauyywe, tauxywe, tauyxwe = [tau2('uu', vdw, n, vde), tau2('vv', vdw, n, vde), tau2('uv', vdw, n, vde),
                                          tau2('vu', vdw, n, vde)]

    SPH = -(tauxx * ux + tauxy * uy + tauyx * vx + tauyy * vy)
    SPHeeE = -(tauxxee * uxe + tauxyee * uye + tauyxee * vxe + tauyyee * vye)
    SPHwwW = -(tauxxww * uxw + tauxyww * uyw + tauyxww * vxw + tauyyww * vyw)
    SPHeeW = -(tauxxee * uxw + tauxyee * uyw + tauyxee * vxw + tauyyee * vyw)
    SPHwwE = -(tauxxww * uxe + tauxyww * uye + tauyxww * vxe + tauyyww * vye)
    SPHweE = -(tauxxwe * uxe + tauxywe * uye + tauyxwe * vxe + tauyywe * vye)
    SPHewE = -(tauxxew * uxe + tauxyew * uye + tauyxew * vxe + tauyyew * vye)
    SPHweW = -(tauxxwe * uxw + tauxywe * uyw + tauyxwe * vxw + tauyywe * vyw)
    SPHewW = -(tauxxew * uxw + tauxyew * uyw + tauyxew * vxw + tauyyew * vyw)

    uze, vze, uzw, vzw, uz, vz = mapfilt([vde['duz'], vde['dvz'], vdw['duz'], vdw['dvz'], vd['duz'], vd['dvz']], n)
    taux, tauy = [tau('uw', vd, n), tau('vw', vd, n)]
    tauxww, tauyww = [tau('uw', vdw, n), tau('vw', vdw, n)]
    tauxee, tauyee = [tau('uw', vde, n), tau('vw', vde, n)]
    tauxew, tauxwe = tau('uw', vde, n, vdw)
    tauyew, tauywe = tau('vw', vde, n, vdw)

    SPVeeE = -(tauxee * uze + tauyee * vze)
    SPVwwW = -(tauxww * uzw + tauyww * vzw)
    SPVeeW = -(tauxee * uzw + tauyee * vzw)
    SPVwwE = -(tauxww * uze + tauyww * vze)
    SPVewE = -(tauxew * uze + tauyew * vze)
    SPVweE = -(tauxwe * uze + tauywe * vze)
    SPVewW = -(tauxew * uzw + tauyew * vzw)
    SPVweW = -(tauxwe * uzw + tauywe * vzw)
    SPV = -(taux * uz + tauy * vz)

    SPVtot = SPVeeE + SPVwwW + SPVeeW + SPVwwE + SPVewW + SPVewE + SPVweW + SPVweE
    SPHtot = SPHeeE + SPHwwW + SPHeeW + SPHwwE + SPHweE + SPHewE + SPHweW + SPHewW

    for var in varlist:
        od[var] = np.mean((locals()[var]))
        # od[var] = np.squeeze((locals()[var]))
    # Print the CPU IDs 1
    ray_nodes = ray.nodes()
    for node in ray_nodes:
        print("Node", node["NodeID"], "is using CPU IDs:", node["Resources"]["CPU"])
    memory_info = ray.cluster_resources()["memory"]
    print("Memory used by Ray: %.1f GB" % (memory_info / 1e9))
    return od


def tau_scales(vd, vdw, vde, param, varlist):
    r = np.zeros(N_filters, dtype=np.float32)
    n_scales = np.zeros(N_filters, dtype=np.float32)
    od = {}

    print('Starting flux loop with %d steps' % (N_filters))
    print('Filter order:', end=' ')
    # Print the CPU IDs 1
    ray_nodes = ray.nodes()
    for node in ray_nodes:
        print("Node", node["NodeID"], "is using CPU IDs:", node["Resources"]["CPU"])

    vd_r = ray.put(vd)
    vdw_r = ray.put(vdw)
    vde_r = ray.put(vde)
    gd_r = ray.put(param)
    memory_info = ray.cluster_resources()["memory"]
    print("Memory used by Ray: %.1f GB" % (memory_info / 1e9))

    dictarr = ray.get([tauijVij.remote(vd_r, vdw_r, vde_r, gd_r, varlist, Narray[i]) for i in range(num_cpus)])

    del vd_r, vde_r, vdw_r, gd_r
    for I, valdict in enumerate(dictarr):
        od[I] = {}
        for var in varlist:
            print(I, var)
            sys.stdout.flush()
            od[I][var] = valdict[var]
            r[I] = Narray[I] * param.DX
            n_scales[I] = Narray[I]
    return n_scales, r, od


######################################
#############MAIN#####################
######################################

# varlist = [ 'SPHscatt', 'SPVscatt', 'SPH', 'SPV']
####################################################################################
try:
    depth = int(sys.argv[1])
    nfilt = list(map(int, sys.argv[2:]))
except:
    depth = 21
    # nfilt = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 112, 144, 176, 240, 368, 496])
    nfilt = np.array([32])
    # sys.exit(0)

# data processing
flux_remainder = 'CG_ray'  # 'vertical' or 'galilean'
directions = ['Horizontal', 'Vertical']
experiments = ['Stochastic']
# experiments = ['Steady']
filter_width1 = 24
filter_width2 = 16
freq_domain1 = 'LF'  # LF
freq_domain2 = 'HF'  # HF
helm_domain1 = None  # rot
helm_domain2 = None  # divs
coor_shift = True
reduce_mean = False
Temporal = False

param = SimulationParameter('Stochastic')

# Output
varlist = ['SPHeeE', 'SPHeeW', 'SPHwwW', 'SPHwwE', 'SPHweW', 'SPHewW', 'SPHewE', 'SPHweE', 'SPVeeE', 'SPVeeW',
           'SPVwwW', 'SPVwwE', 'SPVweW', 'SPVewW', 'SPVewE', 'SPVweE', 'SPH', 'SPHtot', 'SPV', 'SPVtot']

print('nfilt:', nfilt)
Nt = param.Nt
Nt0 = 0
NtEW = 722
print(NtEW)

######################################
print(os.uname()[1])
Narray = nfilt
N_filters = len(Narray)
print('Narray:', np.array(Narray) * param.DX / 1000, len(Narray))
num_cpus = N_filters
destination_log_dir = "/analysis/michalshaham/loggings/"
Path(destination_log_dir).mkdir(parents=True, exist_ok=True)


for exp in experiments:

    if exp == 'Stochastic' and depth == 77:
        continue

    data_sets = TwoDataSets(filt_width1=filter_width1, freq_domain1=freq_domain1, helm_domain1=helm_domain1,
                            filt_width2=filter_width2, freq_domain2=freq_domain2, helm_domain2=helm_domain2,
                            coor_shift=coor_shift, exp=exp, depth=depth, reduce_mean=reduce_mean)

    start_time_ini = time.time()
    file_name = get_CG_file_name(depth, exp, filt1=filter_width1, filt2=filter_width2,
                                 domain1=freq_domain1, domain2=freq_domain2, helm1=helm_domain1,
                                 helm2=helm_domain2, coor_shift=coor_shift, reduce_mean=reduce_mean,
                                 temporal=Temporal, flux_remainder=flux_remainder)
    # file_name = file_name.replace('XY','Kh_TEST/XY')
    file_name = file_name.replace('XY','Kh/XY')
    initialize_coarse_graining_ray_kh(file_name)
    root_grp = open_netcdf(file_name)

    for var_str in varlist:
        print(var_str)
        if not variable_exist(file_name, var_str):
            root_grp.createVariable(var_str, 'f8', ('t', 'k'))
        else:
            print('Already Exists: {}, {}'.format(file_name, var_str))
    root_grp.close()

    print()
    print('TIME for initializing: %.1fsec' % (time.time() - start_time_ini))
    print()
    start_time_ray_init = time.time()
    init = ray.init(num_cpus=num_cpus, object_store_memory=6e9, _memory=6e9,logging_level=logging.ERROR, log_to_driver=False)
    print(f'Number of usable cores: {num_cpus}')
    print('TIME for ray init%.1fsec' % (time.time() - start_time_ray_init))
    print()

    t_jump = 1
    start_time = time.time()
    for itime in range(Nt0, NtEW, t_jump):
        print('STARTING TIME ', itime)

        start_time_step = time.time()
        data_sets.data_set1.ind_min = itime
        data_sets.data_set1.ind_max = itime + t_jump
        data_sets.data_set2.ind_min = itime
        data_sets.data_set2.ind_max = itime + t_jump

        start_time_load = time.time()
        vd = loadvd(data_sets, flag=0)
        vde = loadvd(data_sets, flag=1)
        vdw = loadvd(data_sets, flag=2)
        print('TIME for load %.1fsec' % (time.time() - start_time_load))
        print()

        start_time_calc = time.time()
        n_scales, filtscales, od = tau_scales(vd, vdw, vde, param, varlist)
        print('TIME for calc %.1fsec' % (time.time() - start_time_calc))
        print()

        start_time_save = time.time()
        root_grp = open_netcdf(file_name)
        for I in range(num_cpus):
            for var_str in varlist:
                root_grp.variables[var_str][itime, nfilt[I]] = od[I][var_str]
        root_grp.close()
            # finalize_2D_coarse_graining_ray(file_name, var_str, od[I][var_str], itime)
        print('TIME for saving: %.1fsec' % (time.time() - start_time_save))
        print()

        print(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
        print('TIME from start: %.1fsec' % (time.time() - start_time))
        print('TIME for step: %.1fsec' % (time.time() - start_time_step))
        print()

    ray.shutdown()
    # Wait for a moment to ensure all logs are written
    time.sleep(5)

    # Move the log files to the destination directory
    redirecting_logs()