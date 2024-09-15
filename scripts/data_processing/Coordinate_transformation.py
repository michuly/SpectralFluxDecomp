import sys
import numpy as np

from data_file_names import get_depths_run, print_start

sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
from open_data import *

# data processing
to_save = True
# experiments=['Stochastic']
experiments=['Stochastic']
filter_widths = [0]
freq_domains = [None]
helm_domains = ['rot']
reduce_mean = False
geostrophic = False
dim = 2
# var_strs = ['w']
var_strs = ['u', 'v', 'w', 'duz', 'dvz', 'dvy', 'dux', 'dvx', 'duy', 'dwz']
# var_strs = ['duz', 'dvz']
# var_strs = ['dwz', 'w']


def Coordinate_shift(data_set, var_str, U, file_name=''):

    if initialize_data(file_name, var_str, description='coor shift U=%.3f' % U, dim=dim):
        print('Transformaing Coor {}...'.format(var_str))

        data_tmp = data_set.return_data(var_str)
        data = np.zeros(np.shape(data_tmp))
        for t in range(data_tmp.shape[0]):
            # if data_tmp.shape[0]==721 and t>=498:
            #     x_0 = int(np.round(U * ((t+1) - 100)))
            # else:
            #     x_0 = int(np.round(U * (t - 100)))  # U [ind/t] [DX/DT] [390m/3600s]s
            x_0 = int(np.round(U * (t - 100)))
            x_0 = x_0 % 512
            print(t, x_0)
            if x_0 > 0:
                data[t, :, 0:(512 - x_0)] = data_tmp[t, :, x_0:512]
                data[t, :, (512 - x_0):512] = data_tmp[t, :, 0:x_0]
            elif x_0 == 0:
                data[t, :, :] = data_tmp[t, :, :]

        finalize_data(file_name, var_str=var_str, data=data, dim=dim)
    # change_variable(file_name, var_str=var_str, value=data)


depths=get_depths_run(sys.argv, dim=dim)
for depth in np.flip(depths):
    for filter_width in filter_widths:
        for exp in experiments:
            for freq_domain in freq_domains:
                for helm_domain in helm_domains:
                    print_start(depth, sys.argv)
                    if exp=='Stochastic' and depth==77 and dim==2:
                        continue
                    for var_str in var_strs:
                        data_set = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm=helm_domain,
                                           geostrophic=geostrophic, dim=dim, reduce_mean=reduce_mean, exp=exp, depth=depth)
                        if helm_domain.lower()=='rot' and (var_str=='dwz' or var_str=='w'):
                            continue
                        file_name = get_data_file_name(exp, filter_width, freq_domain, helm_domain, geostrophic, dim=dim,
                                                       roy_folder=False, coor_shift=True).format(f'{depth:03d}')
                        print(file_name)
                        sys.stdout.flush()
                        U = get_coordinate_shift(exp)
                        Coordinate_shift(data_set, var_str, U, file_name=file_name)
                        sys.stdout.flush()

                if filter_width == 0:
                    break

print('BYE BYE')
sys.stdout.flush()
