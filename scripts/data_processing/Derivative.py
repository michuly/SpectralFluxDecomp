import sys

from data_file_names import get_depths_run, print_start

sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
from open_data import *

# data processing
to_save = True
experiments = ['Stochastic', 'Steady']
experiments = ['Stochastic']
filter_widths = [0]
freq_domains = [None]
helm_domains = ['rot']
reduce_mean = False
geostrophic = False
dim = 2
coor_shift = False
var_strs = ['dvy', 'dux', 'dvx', 'duy']
# var_strs = ['dpy', 'ddpyz']
if len(sys.argv) > 3:
    experiments = [str(sys.argv[3])]


def derivative_and_save(data_set, var_str, file_name=''):
    check_var = 'phi' if var_str in ['dpy', 'dpx'] else var_str[1:-1]
    if not data_set.check_variable_exist(check_var):
        print('NO VARIABLE TO DERIVE')
        return False

    if initialize_data(file_name, var_str, description='rotational velocities', dim=dim, exp=exp):
        print(file_name)
        print('Calc derivatives {}...'.format(var_str))
        if var_str == 'dux':
            data = derivative(data_set.return_data('u'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X, flag=0)
        elif var_str == 'dvx':
            data = derivative(data_set.return_data('v'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X, flag=0)
        elif var_str == 'dwx':
            data = derivative(data_set.return_data('w'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X, flag=0)
        elif var_str == 'dvy':
            if data_set.geostrophic:
                flag = 1
            else:
                flag = -1
            data = derivative(data_set.return_data('v'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y,
                              flag=flag)
        elif var_str == 'duy':
            if data_set.geostrophic:
                flag = -1
            else:
                flag = 1
            data = derivative(data_set.return_data('u'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y,
                              flag=flag)
        elif var_str == 'dwy':
            data = derivative(data_set.return_data('w'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y, flag=1)
        elif var_str == 'dbx':
            data = derivative(data_set.return_data('b'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X, flag=0)
        elif var_str == 'dby':
            data = derivative(data_set.return_data('b'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y, flag=1)
        elif var_str == 'dpx':
            if dim == 2:
                data = derivative(data_set.return_data('p'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_X,
                                  flag=0)
            elif dim == 3:
                data = derivative(data_set.return_data('phi'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_3d_X,
                                  flag=0)
        elif var_str == 'dpy':
            print('IS THE FLAG 1 OR -1?')
            if dim == 2:
                data = derivative(data_set.return_data('p'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_Y,
                                  flag=1)
            elif dim == 3:
                data = derivative(data_set.return_data('phi'), L=data_set.L, dx=data_set.DX, n=1, axis=data_set.AX_3d_Y,
                                  flag=1)
        elif var_str == 'ddpxz':
            data = derivative(data_set.return_data('dpx'), L=data_set.H, dx=data_set.DZ, n=1, axis=data_set.AX_3d_Z,
                              flag=1)
        elif var_str == 'ddpyz':
            data = derivative(data_set.return_data('dpx'), L=data_set.H, dx=data_set.DZ, n=1, axis=data_set.AX_3d_Z,
                              flag=1)
        if np.ndim(data) == 4 and dim == 3:
            data = data[:, 0, :, :]
        print(np.shape(data))
        finalize_data(file_name, var_str=var_str, data=data, dim=dim)
        # change_variable(file_name, var_str=var_str, value=data)


depths = get_depths_run(sys.argv, dim=dim)
for depth in np.flip(depths):
    for filter_width in filter_widths:
        for exp in experiments:
            for freq_domain in freq_domains:
                for helm_domain in helm_domains:
                    print_start(depth, sys.argv)
                    if exp == 'Stochastic' and depth == 77 and dim == 2:
                        continue
                    for var_str in var_strs:
                        data_set = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm=helm_domain,
                                           geostrophic=geostrophic, dim=dim, reduce_mean=reduce_mean, exp=exp,
                                           depth=depth, coor_shift=coor_shift)
                        file_name = get_data_file_name(exp, filter_width, freq_domain, helm_domain, geostrophic,
                                                       dim=dim,
                                                       roy_folder=False, coor_shift=coor_shift).format(f'{depth:03d}')

                        print(file_name)
                        # print(data_set.find_self_file_name(var_str[1]))
                        sys.stdout.flush()

                        derivative_and_save(data_set, var_str, file_name=file_name)
                        sys.stdout.flush()

                if filter_width == 0:
                    break

print('BYE BYE')
sys.stdout.flush()
