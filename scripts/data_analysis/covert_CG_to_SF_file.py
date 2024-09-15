#!/usr/bin/env python
#################################################
# Load Modules
#################################################
from open_data import *

def Forder(var):
    return np.asfortranarray(var.T, dtype=np.float32)
    # return var


def ncopen(filename, depth, opt='r'):
    return Dataset(filename + '.{0:04}'.format(depth) + '.nc', opt)


def flags_to_var_name(direction, flags):
    var_str = 'SP'
    if direction == 'Vertical':
        var_str += 'V'
    elif direction == 'Horizontal':
        var_str += 'H'
    if flags == (0, 0, 0):
        return var_str
    elif flags == (3, 3, 3):
        return var_str + 'tot'
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


######################################
#############MAIN#####################
######################################
# Output
varlist = ['SPHeeE', 'SPHeeW', 'SPHwwW', 'SPHwwE', 'SPHweW', 'SPHewW', 'SPHewE', 'SPHweE', 'SPVeeE', 'SPVeeW',
           'SPVwwW', 'SPVwwE', 'SPVweW', 'SPVewW', 'SPVewE', 'SPVweE', 'SPH', 'SPHtot', 'SPV', 'SPVtot']

# varlist = ['SPHeeW', 'SPHwwE', 'SPHewW', 'SPHweE', 'SPVeeW', 'SPVwwE', 'SPVewW', 'SPVweE', 'SPH', 'SPV']

# varlist = [ 'SPHscatt', 'SPVscatt', 'SPH', 'SPV']
####################################################################################

depths = np.arange(1,130)
nfilt = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 112, 144, 176, 240, 368, 496])
# data processing
flux_remainder = 'CG_ray'  # 'vertical' or 'galilean'
directions = ['Horizontal', 'Vertical']
experiments = ['Stochastic']
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

flags_list = [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 2, 2), (2, 2, 1), (2, 1, 1), (2, 1, 2), (0, 0, 0),
              (3, 3, 3)]
# flags_list = [(1, 2, 2), (2,2,1)]
flags_list = [(2,2,2)]



######################################
Narray = nfilt
N_filters = len(Narray)
print('Narray:', np.array(Narray) * param.DX / 1000, len(Narray))


for exp in experiments:
    for direction in directions:
        for flags in flags_list:
            var_str = flags_to_var_name(direction, flags)
            file_name_SF = get_SF_file_name(flags, exp, direction, filt1=filter_width1, filt2=filter_width2,
                                            domain1=freq_domain1, domain2=freq_domain2, helm1=helm_domain1,
                                            helm2=helm_domain2, coor_shift=coor_shift, reduce_mean=reduce_mean,
                                            temporal=Temporal, flux_type=flux_remainder)
            print(file_name_SF)

            for depth in depths:
                if exp == 'Stochastic' and depth == 77:
                    continue
                initialize_coarse_graining(file_name_SF, axes_keys=('freq',), axes_values=(), shape=(512,),
                                           depth=depth, wavelength_str='n', freq_str='freq_h', wavelength_max=512)

                file_name_CG = get_CG_file_name(depth, exp, filt1=filter_width1, filt2=filter_width2,
                                         domain1=freq_domain1, domain2=freq_domain2, helm1=helm_domain1,
                                         helm2=helm_domain2, coor_shift=coor_shift, reduce_mean=reduce_mean,
                                         temporal=Temporal, flux_remainder=flux_remainder)
                file_name_CG = file_name_CG.replace('XY','Kh/XY')
                if type(depth) == str:
                    depth_str = depth
                else:
                    depth_str = f'{depth:03d}'

                try:
                    with open_netcdf(file_name_CG, 'r') as dat:
                        data_new = np.nanmean(dat.variables[var_str][:], 0)
                    root_grp = open_netcdf(file_name_SF)
                    if variable_exist(file_name_SF, depth_str, check_nan=False):
                        data = root_grp.variables[depth_str]
                        data[:] = data_new
                    else:
                        create_variable(file_name_SF, var_str, data_new, axes=('freq',))
                    root_grp.close()
                except:
                    print('There is no ', exp, direction, depth_str, var_str)





