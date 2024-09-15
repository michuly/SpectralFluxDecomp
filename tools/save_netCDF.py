from netCDF4 import Dataset
import numpy as np
import os
import time
import sys


def create_netcdf(name='../../Analysis/example.nc', description='', axes_keys=('time', 'y', 'x'), axes_values=None,
                  shape=(3, 3, 4), fields_keys=('u', 'v'),
                  fields_values=(np.random.random((3, 3, 4)), np.random.random((3, 3, 4))),
                  verbose=True):
    root_grp = open_netcdf(name, verbose=verbose)
    if description:
        root_grp.description = description

    # axes
    for i in range(len(axes_keys)):
        if axes_keys[i] not in root_grp.dimensions.keys():
            root_grp.createDimension(axes_keys[i], shape[i])
            if not axes_values is None:
                if verbose:
                    print('Saving axes {}, {}'.format(name, axes_keys[i], shape[i]))
                variable = root_grp.createVariable(axes_keys[i], 'f4', (axes_keys[i],))
                variable[:] = axes_values[i]
            # elif dimensions is not None:
            #     if verbose:
            #         print('Saving axes {}, {}'.format(name, dimensions[i], shape[i]))
            #     variable = root_grp.createVariable(axes_keys[i], 'f4', (axes_keys[i],))
            #     variable[:] = np.linspace(0, dimensions[i], shape[i])

    # fields
    for j in range(len(fields_keys)):
        if fields_keys[j] not in root_grp.variables.keys():
            if verbose:
                print('Saving field {}, {}'.format(name, fields_keys[j]))
            data = root_grp.createVariable(fields_keys[j], 'f8', axes_keys)
            data[:] = fields_values[j]

    root_grp.close()


def open_netcdf(name, verbose=True):
    max_tries = 3
    try_num = 0
    while try_num < max_tries:
        try:
            print('Accessing file {}, try_num {} ...'.format(name, try_num))
            if os.path.exists(name):
                root_grp = Dataset(name, 'a')
            else:
                if verbose:
                    print('Creating NCfile {}'.format(name))
                path_name = os.path.dirname(name)
                if not os.path.isdir(path_name):
                    print('Creating Directory: ', path_name)
                    os.makedirs(path_name)
                root_grp = Dataset(name, 'w', format='NETCDF4')
            break
        except PermissionError as err:
            print('CANNOT OPEN NETCDF')
            print(err)
            time.sleep(10)
            try_num += 1
        except OSError as err:
            print('CANNOT OPEN NETCDF')
            print(err)
            time.sleep(10)
            try_num += 1
    if try_num == 3:
        raise OSError('CANNOT OPEN NETCDF {}'.format(name))
    return root_grp


def read_netcdf(name, verbose=True):
    max_tries = 3
    try_num = 0
    while try_num < max_tries:
        try:
            print('Accessing file {}, try_num {} ...'.format(name, try_num))
            root_grp = Dataset(name, 'r')
            break
        except PermissionError as err:
            print('CANNOT OPEN NETCDF')
            print(err)
            time.sleep(10)
            try_num += 1
        except OSError as err:
            print('CANNOT OPEN NETCDF')
            print(err)
            time.sleep(10)
            try_num += 1
    if try_num == 3:
        raise OSError('CANNOT OPEN NETCDF {}'.format(name))
    return root_grp


def variable_exist(file_name, var_str, check_nan=True):
    possible_machine = ['/atlantic/', '/atlantic2/']

    for machine in possible_machine:
        file_name = file_name.replace(machine,'{}')
    for machine in possible_machine:
        file_name_tmp = file_name.format(machine)
        if os.path.exists(file_name_tmp):
            root_grp = read_netcdf(file_name_tmp)
            if var_str in root_grp.variables.keys():
                # if not root_grp.variables[var_str][:]:
                #     root_grp.close()
                if check_nan:
                    if np.any(np.isnan(root_grp.variables[var_str][:])):
                        root_grp.close()
                        return False
                root_grp.close()
                return file_name_tmp
            else:
                root_grp.close()
    return False


def file_exist(file_name):
    possible_machine = ['/atlantic/', '/atlantic2/']

    for machine in possible_machine:
        file_name = file_name.replace(machine,'{}')
    for machine in possible_machine:
        file_name_tmp = file_name.format(machine)
        if os.path.exists(file_name_tmp):
            return file_name_tmp
    return False


def create_variable(file_name, var_str, value, axes, verbose=True):
    if os.path.exists(file_name):
        if variable_exist(file_name, var_str):
            if verbose:
                print('Already Exists: {}, {}'.format(file_name, var_str))
        else:
            if verbose:
                print('Saving {}, {}'.format(file_name, var_str))
            root_grp = open_netcdf(file_name)
            data = root_grp.createVariable(var_str, 'f8', axes)
            data[:] = value
            root_grp.close()
            return True

    return False


def change_variable(file_name, var_str, value, slicer=None, verbose=True):
    if verbose:
        print('Changing Variable: ', file_name, var_str)
    root_grp = open_netcdf(file_name)
    if slicer is None:
        root_grp.variables[var_str][:] = value
    else:
        root_grp.variables[var_str][tuple(slicer)] = value
    root_grp.close()


def initialize_netcdf(file_name, description, axes_keys, axes_values, shape):
    if not os.path.exists(file_name):
        create_netcdf(file_name, description, axes_keys, axes_values, shape, fields_keys=(), fields_values=())
    for ax_key, ax_value in zip(axes_keys, axes_values):
        if not variable_exist(file_name, ax_key):
            create_variable(file_name=file_name, var_str=ax_key, value=ax_value, axes=ax_key)
        # if freq is not None:
        #     if type(freq) != tuple:
        #         freq = (freq,)
        #         freq_str = (freq_str,)
        #     variable_axes = freq_str
        #     #$$$ fix it - it doesn't use the freq I give it!!!
        #     for i in range(len(freq)):
        #         if not variable_exist(file_name, freq_str[i]):
        #             create_netcdf(name=file_name, fields_keys=(freq_str[i],), fields_values=(freq[i],),
        #                           variable_axes=(variable_axes[i],), axes_keys=())


def initialize_spectral_flux(file_name, description, axes_keys, axes_values, shape, depth, verbose=True):
    var_str = f'{depth:03d}'
    tested_file_name = variable_exist(file_name, var_str)
    if tested_file_name:
        print('already exist: {} {}'.format(tested_file_name, var_str))
        return False
    initialize_netcdf(file_name, description, axes_keys, axes_values, shape)
    if verbose:
        print('Initialzing Variable: ', file_name, var_str)
    return True


def initialize_data(file_name, var_str, description, dim=2, exp='Stochastic'):
    if not os.path.exists(file_name):
        if dim==3:
            if exp != 'Stochastic_only':
                shape = (129, 513, 512)
                dimensions = (2000, 200000, 200000)
            axes_keys = ('z', 'y', 'x')
        else:
            if exp == 'Steady':
                shape = (722, 513, 512)
                dimensions = (722, 200000, 200000)
            elif exp == 'Stochastic':
                shape = (722, 513, 512)
                dimensions = (722, 200000, 200000)
            elif exp == 'Stochastic_only':
                shape = (2223, 257, 256)
                dimensions = (2222, 200000, 200000)
            axes_keys = ('t', 'y', 'x')
        axes_values = [np.linspace(0, dimensions[i], shape[i]) for i in range(3)]
        initialize_netcdf(file_name, description, axes_keys, axes_values, shape)
    if var_str:
        if variable_exist(file_name, var_str, check_nan=True):
            print('already exist: {} {}'.format(file_name, var_str))
            return False

    return True


def initialize_2D_coarse_graining_ray(file_name, shape=(722, 512, 512)):
    if not os.path.exists(file_name):
        dimensions = (722, 200000, 200000)
        axes_keys = ('t', 'y', 'x')
        axes_values = [np.linspace(0, dimensions[i], shape[i]) for i in range(3)]
        initialize_netcdf(file_name, '', axes_keys, axes_values, shape)

    return True


def initialize_coarse_graining_ray_kh(file_name, shape=(722, 512)):
    if not os.path.exists(file_name):
        dimensions = (722, 200000)
        axes_keys = ('t', 'k')
        axes_values = [np.linspace(0, dimensions[i], shape[i]) for i in range(2)]
        initialize_netcdf(file_name, '', axes_keys, axes_values, shape)

    return True

def finalize_data(file_name, var_str, data, dim=2):
    if dim==3:
        create_variable(file_name=file_name, var_str=var_str, value=data, axes=('z', 'y', 'x'))
    else:
        create_variable(file_name=file_name, var_str=var_str, value=data, axes=('t', 'y', 'x'))


def finalize_spectral_flux(file_name, data, depth_, axes_keys=('freq',)):
    if type(depth_) == str:
        var_str = depth_
    else:
        var_str = f'{depth_:03d}'
    create_variable(file_name=file_name, var_str=var_str, value=data, axes=axes_keys)


def initialize_coarse_graining(file_name, axes_keys, axes_values, shape, depth, wavelength=None, wavelength_str='n',
                               freq_str='freq_h', wavelength_max=512):
    if not os.path.exists(file_name):
        initialize_netcdf(file_name, '', axes_keys, axes_values, shape)
    if not variable_exist(file_name, var_str=wavelength_str, check_nan=False):
        create_variable(file_name, var_str=wavelength_str, value=np.arange(1,wavelength_max+1), axes=axes_keys)
    if not variable_exist(file_name, var_str=freq_str, check_nan=False):
        create_variable(file_name, var_str=freq_str, value=1/np.arange(1,wavelength_max+1), axes=axes_keys)

    if type(depth) == str:
        var_str = depth
    else:
        var_str = f'{depth:03d}'
    confirmed_file_name = variable_exist(file_name, var_str, check_nan=False)
    if confirmed_file_name:
        root_grp = read_netcdf(confirmed_file_name)
        if wavelength is not None and not np.isnan(root_grp.variables[var_str][wavelength-1]):
            print('Already Exist: ', file_name, var_str, wavelength)
            root_grp.close()
            return False
        else:
            root_grp.close()
            return True
    else:
        nan_vec = np.empty(wavelength_max)
        nan_vec.fill(np.nan)
        create_variable(file_name, var_str, nan_vec, axes_keys)
        print('Initialzing Variable: ', file_name, var_str, wavelength)
        return True


def finalize_coarse_graining(file_name, CG, depth, wavelength):
    if type(depth) == str:
        var_str = depth
    else:
        var_str = f'{depth:03d}'

    if os.path.exists(file_name):
        root_grp = open_netcdf(file_name)
        if variable_exist(file_name, var_str, check_nan=False):
            data = root_grp.variables[var_str]
            if np.isnan(data[wavelength-1]):
                data[wavelength-1] = CG
                root_grp.close()
                return True

        root_grp.close()
    print('Could NOT save CG')
    return False


def initialize_coarse_graining_ray(file_name, axes_keys, axes_values, shape, depth, wavelength_str='n',
                               freq_str='freq_h', wavelength_max=512, time_max=722):
    if not os.path.exists(file_name):
        initialize_netcdf(file_name, '', axes_keys, axes_values, shape)
    if not variable_exist(file_name, var_str=wavelength_str, check_nan=False):
        create_variable(file_name, var_str=wavelength_str, value=np.arange(1,wavelength_max+1), axes=(axes_keys[1],))
    if not variable_exist(file_name, var_str=freq_str, check_nan=False):
        create_variable(file_name, var_str=freq_str, value=1/np.arange(1,wavelength_max+1), axes=(axes_keys[1],))
    if not variable_exist(file_name, var_str='time', check_nan=False):
        create_variable(file_name, var_str='time', value=np.arange(time_max), axes=(axes_keys[0],))

    if type(depth) == int:
        var_str = f'{depth:03d}'
    else: var_str=depth
    if not variable_exist(file_name, var_str, check_nan=False):
        nan_vec = np.empty((time_max, wavelength_max))
        nan_vec.fill(np.nan)
        create_variable(file_name, var_str, nan_vec, axes_keys)
        print('Initialzing Variable: ', file_name, var_str)
        return True


def finalize_coarse_graining_ray(file_name, CG_ray, depth, wavelength, itime):
    print('SAVING: ', depth, itime, wavelength, CG_ray, file_name)
    if type(depth) == str:
        var_str = depth
    elif type(depth) == int:
        var_str = f'{depth:03d}'
    root_grp = open_netcdf(file_name)
    root_grp.variables[var_str][itime, wavelength-1] = CG_ray
    root_grp.close()


def finalize_2D_coarse_graining_ray(file_name, var_str, CG_ray, itime, ftime=None):
    print('SAVING: ', var_str, file_name)
    root_grp = open_netcdf(file_name)
    print('SHAPE', var_str, CG_ray.shape)
    if not variable_exist(file_name, var_str):
        root_grp.createVariable(var_str, 'f8', ('t', 'y', 'x'))
    else:
        print('Already Exists: {}, {}'.format(file_name, var_str))
    if ftime:
        root_grp.variables[var_str][itime:ftime, :, :] = CG_ray
    else:
        root_grp.variables[var_str][itime, :, :] = CG_ray
    root_grp.close()


def time_mean_coarse_graining_ray(file_name, axes_keys, axes_values, shape, depth, wavelength_str='n',
                               freq_str='freq_h', wavelength_max=512):

    file_name_new = file_name.replace('/TIME', '')
    initialize_coarse_graining(file_name, depth, axes_keys, axes_values, shape, wavelength_str,
                               freq_str, wavelength_max)
    if type(depth) == str:
        var_str = depth
    else:
        var_str = f'{depth:03d}'

    with open_netcdf(file_name) as dat:
        data_new = np.nanmean(dat.variables[var_str][:],0)

    root_grp = open_netcdf(file_name_new)
    if variable_exist(file_name, var_str, check_nan=False):
        data = root_grp.variables[var_str]
        data[:] = data_new
    else:
        create_variable(file_name, var_str, data_new, axes_keys)

    root_grp.close()




