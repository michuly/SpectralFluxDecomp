# import sys
# sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
from data_file_names import get_depths_run, print_start
from open_data import *
from tools.simulation_analysis import helmholz_decomposition_netcdf

dbz=False

# data processing
experiments = ['Stochastic']
filter_widths = [0]
freq_domains = [None]
helms=['rot']
vertical_derivatives = [False, True]
sys.stdout.flush()

def helmholz_and_save(data_set, file_name='', vertical_derivative=False, helm='rot'):
    if vertical_derivative:
        str_format = 'd%sz'
    else:
        str_format = '%s'

    if (initialize_data(file_name, str_format % 'v', description='rotational velocities') or
            initialize_data(file_name, str_format % 'u', description='rotational velocities')):

        # print('CALCulating vort or div')
        # if helm.lower() == 'div':
        #     print('Divergence')
        #     data_set2 = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm='div', exp=exp,
        #                        depth=depth)
        #     data = vorticity_netcdf(data_set2, check_existing=False)
        # elif helm.lower() == 'rot':
        #     print('Rotation')
        #     data_set2 = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm='rot', exp=exp,
        #                        depth=depth)
        #     data = divergence_netcdf(data_set2)
        # print(np.mean(np.abs(data)))
        # if np.mean(np.abs(data))>1e-14:
        #     print('CHECK IF BUG!!!')

        print('Calc helm...')
        u_rot, v_rot = helmholz_decomposition_netcdf(data_set, vertical_derivative=vertical_derivative)
        if helm.lower() == 'rot':
            finalize_data(file_name, var_str=str_format % 'u', data=u_rot)
            finalize_data(file_name, var_str=str_format % 'v', data=v_rot)
            # data = divergence(u_rot, v_rot, data_set.L, data_set.DX, data_set.AX_X, data_set.AX_Y)
            # print(np.mean(np.abs(data)))
        elif helm.lower() == 'div':
            u = data_set.return_data('u')
            v = data_set.return_data('v')
            finalize_data(file_name, var_str=str_format % 'u', data=u-u_rot)
            finalize_data(file_name, var_str=str_format % 'v', data=v-v_rot)
            # data = vorticity(u-u_rot, v-v_rot, data_set.L, data_set.DX, data_set.AX_X, data_set.AX_Y)
            # print(np.mean(np.abs(data)))


def dbz_and_save(data_set, file_name=''):

    if initialize_data(file_name, 'u', description='rotational velocities'):
        print('Calc helm...')
        u_dbz, v_dbz = dbz_decomposition_netcdf(data_set)
        finalize_data(file_name, var_str='u', data=u_dbz)
        finalize_data(file_name, var_str='v', data=v_dbz)


depths=get_depths_run(sys.argv, experiments[0])
for vertical_derivative in vertical_derivatives:
    for depth in np.flip(depths):
        for exp in experiments:
            for helm, filter_width, freq_domain in zip(helms, filter_widths, freq_domains):
                if exp=='Stochastic' and depth==77:
                    continue
                file_name = get_data_file_name(exp, filter_width, freq_domain, helm=helm, dbz_source=dbz).format(f'{depth:03d}')

                print_start(depth, sys.argv)
                data_set = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm=None, exp=exp,
                                   depth=depth)
                print(file_name)
                sys.stdout.flush()

                helmholz_and_save(data_set, file_name=file_name, vertical_derivative=vertical_derivative, helm=helm)
                sys.stdout.flush()

                if filter_width == 0:
                    break

print('BYE BYE')
sys.stdout.flush()
