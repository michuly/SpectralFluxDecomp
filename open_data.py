from scipy.interpolate import interp1d

from data_file_names import *
from tools.save_netCDF import *
from tools.simulation_analysis import *

class SimulationParameter(object):
    def __init__(self, exp):
        if exp in ['Stochastic', 'Steady']:
            self.SHAPE = np.array([722, 513, 512])  # time, y, x
            self.DX = 390
            self.DZ = 15
        elif exp == 'Stochastic_only':
            self.SHAPE = np.array([2223, 257, 256])  # time, y, x
            self.DX = 390 * 2
        self.DT = 3600
        self.DX_SHAPE = np.array([self.DT, self.DX, self.DX])
        self.DX_SHAPE_PSD = np.array([1, self.DX / 1000, self.DX / 1000])
        self.L = 200000
        self.H = 2000
        self.AX_T = 0
        self.AX_X = 2
        self.AX_Y = 1
        self.AX_Z = 0
        self.AX_3d_X = 2
        self.AX_3d_Y = 3
        self.AX_3d_Z = 0
        self.Nt = 722
        self.Ny = 513
        self.Nx = 512
        self.NZ = 129


class DataSet(SimulationParameter):
    def __init__(self, filter_width=None, freq_domain=None, helm=None, geostrophic=False, exp=None, depth=None,
                 dbz_source=False, filter_order=3, reduce_mean=False, verbose=True, dim=2, coor_shift=False, ind_min=0,
                 ind_max=722):
        SimulationParameter.__init__(self, exp)
        # print('DEFINING DATASET: ', filt_width, freq_domain, helm_domain, exp, depth, reduce_mean, ind_min, ind_max)
        self.ind_min = ind_min
        if exp == 'Stochastic_only' and ind_max == 722:
            self.ind_max = 2223
        else:
            self.ind_max = ind_max
            self.SHAPE[0] = ind_max-ind_min
        self.filter_width = filter_width
        self.freq_domain = freq_domain
        self.helm_domain = helm
        self.dbz_source = dbz_source  # this is the point source of the vorticity in the PV=0 (dbz/N^2)
        self.reduce_mean = reduce_mean
        self.original = True if (not filter_width and not helm and not geostrophic and coor_shift) else False
        self.dim = dim
        self.exp = exp
        self.depth = depth
        self.filter_order = filter_order
        self.geostrophic = geostrophic
        self.U = get_coordinate_shift(self.exp) if coor_shift else 0

        print('coor shift is: ', self.U)
        if exp not in ['Stochastic', 'Steady', 'Stochastic_only']:
            print('EXPERIMENT NAME IS NOT CORRECT')
        if verbose:
            print('Defining dataset :', exp, depth)

    def get_dataset(self, filter_width, freq_domain, rotational, direction, dbz_source=False, verbose=True,
                    dim=2, filter_order=3, roy_folder=True):
        if dim==2 or not roy_folder:
            num_str = '%03d' % self.depth
        elif dim==3:
            if self.exp=='Stochastic':
                num_str = '%06d' % (self.depth * 180)
            elif self.exp=='Steady':
                num_str = '%06d' % (self.depth * 90)
        file_name = get_data_file_name(exp=self.exp, filter_width=filter_width, freq_domain=freq_domain,
                                       helm=rotational, dbz_source=dbz_source, direction=direction, dim=dim,
                                       filter_order=filter_order, roy_folder=roy_folder).format(num_str)

        if not os.path.exists(file_name):
            file_name = file_name.replace('atlantic2', 'atlantic')
            if not os.path.exists(file_name):
                print('FILE DOES NOT EXIST %s' % file_name)
                return False

        if verbose:
            print('Opening Dataset: ', file_name)
        return Dataset(file_name, 'r')

    def get_data_stochastic_only(self, var_str):
        if var_str == 'u':
            var_str = 'v'
        elif var_str == 'v':
            var_str = 'w'
        if self.depth == 65:
            depth_str = 'middepth'
        elif self.depth == 119:
            depth_str = 'ml'
        elif self.depth == 128:
            depth_str = 'surface'
        else:
            print('DEPTH IS NOT CORRECT')
            return
        dat = np.zeros((2223, 257, 256))
        dat.fill(np.NaN)
        print('Starting arranging velocity in only stochastic solution')
        for file_name in glob.glob(data_format_stochastic_only.format(depth_str+'_{}').format("00[0-2][0-9][0-9][0-9]")):
            dat_tmp = Dataset(file_name)
            num = int(re.findall("\d+", file_name)[-1])
            dat[num, :, :] = dat_tmp.variables[var_str][:, 0, :]
        if np.any(np.isnan(dat)):
            print('THE ARRAY', var_str, 'HAS Nan')
        else:
            print('mean', var_str, ' :', np.mean(np.abs(dat)))

        return dat

    def find_self_file_name(self, var_str, check_exist=True):
        """
        returns file_name of netcdf if the variable exist
        if not, returns False
        """
        return self.find_file_name(var_str, filter_width=self.filter_width, freq_domain=self.freq_domain,
                                   helm=self.helm_domain, geostrophic=self.geostrophic, U=self.U,
                                   dbz_source=self.dbz_source, filter_order=self.filter_order, check_exist=check_exist)

    def find_file_name(self, _var_str, filter_width=None, freq_domain=None, helm=False, geostrophic=False, U=0,
                       dbz_source=False, filter_order=3, check_exist=True):
        """
        returns file_name of netcdf if the variable exist
        if not, returns False
        """
        num_str = '%03d' % self.depth
        shift_coor = False if U == 0 else True

        direction = None
        if _var_str in ['u', 'v', 'w', 'duz', 'dvz', 'dwz', 'b']:
            direction = 'Vertical'
        elif _var_str in ['dvy', 'dux', 'dvx', 'duy', 'dwx', 'dwy']:
            direction = 'Horizontal'
        elif _var_str in ['dbz']:
            direction = 'dbz'

        if _var_str in ['Vort', 'Div', 'Cu', 'Ri', 'Sigma']:
            file_name = get_analysis_file_name(exp=self.exp, filt=filter_width, domain=freq_domain, helm=helm,
                                               geostrophic=geostrophic, coor_shift=shift_coor, dim=self.dim,
                                               reduce_mean=False).format(num_str)
        else:
            if self.exp=='Stochastic_only' and _var_str not in ['u','v','w']:
                roy_folder = False
            elif (_var_str in ['dbx', 'dby', 'p', 'dpx', 'dpy', 'ddpxz', 'ddpyz'] or self.U!=0 or
                  (self.exp=='Stochastic' and _var_str in ['dux', 'duy', 'dvx', 'dvy'])): #(maybe geostrophic also?)
                roy_folder = False
            else:
                roy_folder = True

            if self.dim == 3 and roy_folder:
                if self.exp == 'Stochastic':
                    num_str = '%06d' % (self.depth * 180)
                elif self.exp == 'Steady':
                    num_str = '%06d' % (self.depth * 90)
            if self.exp == 'Stochastic_only' and roy_folder:
                if self.depth == 65:
                    num_str = 'middepth'
                elif self.depth == 119:
                    num_str = 'ml'
                elif self.depth == 128:
                    num_str = 'surface'
                num_str += '_000000'

            file_name = get_data_file_name(exp=self.exp, filter_width=filter_width, freq_domain=freq_domain,
                                           helm=helm, geostrophic=geostrophic, coor_shift=shift_coor, dbz_source=dbz_source,
                                           direction=direction, dim=self.dim, filter_order=filter_order,
                                           roy_folder=roy_folder).format(num_str)
            # print(file_name)

        if not check_exist:
            return file_name
        file_name = file_name.replace('atlantic2', os.uname()[1])
        if os.path.exists(file_name):
            dat = Dataset(file_name, 'r')
            if _var_str in dat.variables.keys():
                dat.close()
                print('OPENING variable %s at file %s' % (_var_str, file_name))
                return file_name
            dat.close()
        file_name = file_name.replace(os.uname()[1], 'atlantic2')
        if os.path.exists(file_name):
            dat = Dataset(file_name, 'r')
            if _var_str in dat.variables.keys():
                dat.close()
                print('OPENING variable %s at file %s' % (_var_str, file_name))
                return file_name
            dat.close()
        file_name = file_name.replace('atlantic2', 'atlantic')
        if os.path.exists(file_name):
            dat = Dataset(file_name, 'r')
            if _var_str in dat.variables.keys():
                dat.close()
                print('OPENING variable %s at file %s' % (_var_str, file_name))
                return file_name
            dat.close()
        print('DATA %s DOES NOT EXIST IN FILE %s' % (_var_str, file_name))
        return False

    def var_exist(self, _var_str, filter_width=None, freq_domain=None, helm=None,
                  dbz_source=False, direction=None, filter_order=3):
        """
        returns True of netcdf if the variable exist
        if not, returns False
        """
        if type(self.find_file_name(_var_str, filter_width=filter_width, freq_domain=freq_domain,
                                    helm=helm, geostrophic=self.geostrophic, U=self.U, dbz_source=dbz_source,
                                    filter_order=filter_order)) == str:
            return True
        else:
            return False

    def check_variable_exist(self, var_str):

        direction = None
        if var_str in ['u', 'v', 'w', 'duz', 'dvz', 'dwz', 'Vort', 'div', 'b']:
            direction = 'Vertical'
        elif var_str in ['dvy', 'dux', 'dvx', 'duy', 'dwx', 'dwy']:
            direction = 'Horizontal'
        elif var_str == 'dbz':
            direction = 'dbz'

        if self.original or var_str == 'dbz' or (not self.filter_width and not self.helm_domain):
            exist = self.var_exist(var_str, direction=direction)

        elif self.dbz_source:
            if self.filter_width and self.freq_domain == 'HF':
                exist = self.var_exist(var_str, self.filter_width, 'HF', None, True)

        elif self.helm_domain:
            if self.filter_width:
                if self.freq_domain == 'LF' or self.freq_domain == 'BP':
                    if self.helm_domain.lower() == 'rot':
                        exist = self.var_exist(var_str, self.filter_width, self.freq_domain, 'Rot')
                    elif self.helm_domain.lower() == 'div':
                        exist = self.var_exist(var_str, self.filter_width, self.freq_domain, 'Div')
                        if not exist:
                            exist = self.var_exist(var_str, self.filter_width, self.freq_domain) and \
                                    self.var_exist(var_str, self.filter_width, self.freq_domain, 'Rot')

                elif self.freq_domain == 'HF':
                    if self.helm_domain.lower() == 'rot':
                        exist = self.var_exist(var_str, self.filter_width, 'HF', 'Rot')
                    elif self.helm_domain.lower() == 'div':
                        if self.U == 0:
                            exist = self.var_exist(var_str, direction=direction) and \
                                    self.var_exist(var_str, self.filter_width, 'LF') and \
                                    self.var_exist(var_str, self.filter_width, 'HF', 'Rot')

                        else:
                            exist = self.var_exist(var_str, self.filter_width, 'HF', 'Div')

            elif not self.filter_width:
                if self.helm_domain.lower() == 'rot':
                    exist = self.var_exist(var_str, helm='Rot')
                elif self.helm_domain.lower() == 'div':
                    exist = self.var_exist(var_str, direction=direction) and self.var_exist(var_str, helm='Rot')

        elif not self.helm_domain:
            if self.freq_domain == 'LF' or self.freq_domain == 'BP':
                exist = self.var_exist(var_str, self.filter_width, self.freq_domain)

            elif self.freq_domain == 'HF':
                exist = self.var_exist(var_str, direction=direction) and self.var_exist(var_str, self.filter_width,
                                                                                        'LF')

        return exist

    def open_data(self, var_str, float_type, filter_width=0, freq_domain=None, helm=None, dbz_source=False):

        if (var_str == 'w' or var_str == 'dwz') and (helm and helm.lower()=='rot'):  # w_rot=0, w_div=w_tot
            return 0
        if (var_str in ['u','v','w']) and self.exp == 'Stochastic_only' and not filter_width and not bool(helm):
            return self.get_data_stochastic_only(var_str)

        filter_order = self.filter_order
        file_name = self.find_file_name(var_str, filter_width=filter_width, freq_domain=freq_domain,
                                        helm=helm, geostrophic=self.geostrophic, dbz_source=dbz_source,
                                        filter_order=filter_order, U=self.U)
        # print('Opening file: ', file_name)
        dat = Dataset(file_name)

        # interpolation: original data (only!) is missing time = 498
        if float_type == 32:
            float_type = np.float32
        elif float_type == 64:
            float_type = np.float64

        data = dat.variables[var_str][self.ind_min:self.ind_max, :, :].astype(float_type)
        dat.close()
        return data

        # missing_ind = 498
        # t_size=dat.variables[var_str].shape[0]
        # if self.ind_max == 722:
        #     _j = t_size
        #
        # data = dat.variables[var_str][self.ind_min:_j, :, :].astype(float_type)
        # if self.exp=='Stochastic' and missing_ind > self.ind_min and missing_ind < self.ind_max and t_size == 722:
        #     data = np.concatenate((data[:498, :, :], data[499:, :, :]), axis=0)
        # dat.close()
        # return data

        # _i = self.ind_min
        # _missing_ind = 498
        # missing_ind = _missing_ind - self.ind_min
        # if 'time' in dat.dimensions.keys() and dat.dimensions['time'].size == 721 and self.ind_max == 722:
        #     _j = 721
        #     check_interpolation = True
        # else:
        #     _j = self.ind_max
        #     check_interpolation = False
        #
        # if float_type == 32:
        #     float_type = np.float32
        # elif float_type == 64:
        #     float_type = np.float64
        #
        # if check_interpolation and _i < _missing_ind < _j:
        #     data = interpolation(dat.variables[var_str][_i:_j, :, :], missing_ind).astype(float_type)
        # else:
        #     data = dat.variables[var_str][_i:_j, :, :].astype(float_type)
        # dat.close()
        #
        # return data

    def return_data(self, var_str, float_type=64):

        # print('Return Variable: ', var_str)
        dat = None

        if var_str=='u' and self.U!=0:
            v=self.U*390/3600
            print('Removing shift velocity of %f from %s, filter is: %s' %
                  (v, var_str, str((self.filter_width, self.freq_domain, self.helm_domain))))
        else:
            v = 0

        if var_str == 'N2':  # need to make N2 a different own data file!!!
            _dataset = Dataset(N2_format['time_avg'].format(self.exp))
            dat = _dataset.variables['N2'][self.depth - 1, :, :]

        elif self.original or (not self.filter_width and not self.helm_domain) or (var_str.lower() in ['div', 'rot']):
            dat = self.open_data(var_str, float_type)

        elif self.dbz_source:
            if self.filter_width and self.freq_domain == 'HF':
                dat = self.open_data(var_str, float_type, self.filter_width, 'HF', None, True)

        elif self.helm_domain:
            if self.filter_width:
                if self.freq_domain == 'LF' or self.freq_domain == 'BP':
                    if self.helm_domain.lower() == 'rot':
                        dat = self.open_data(var_str, float_type, self.filter_width, self.freq_domain, 'Rot')
                    elif self.helm_domain.lower() == 'div':
                        if self.var_exist(var_str, self.filter_width, self.freq_domain, helm='Div'):
                            dat = self.open_data(var_str, float_type, self.filter_width, self.freq_domain, 'Div')
                        else:
                            dat = self.open_data(var_str, float_type, self.filter_width, self.freq_domain) - \
                                  self.open_data(var_str, float_type, self.filter_width, self.freq_domain, 'Rot')
                        if v!=0:
                            dat = dat - v

                elif self.freq_domain == 'HF':
                    if self.helm_domain.lower() == 'rot':
                        dat = self.open_data(var_str, float_type, self.filter_width, 'HF', 'Rot')
                    elif self.helm_domain.lower() == 'div':
                        if self.var_exist(var_str, self.filter_width, 'HF', helm='Div'):
                            dat = self.open_data(var_str, float_type, self.filter_width, 'HF', helm='Div')
                        else:
                            dat = self.open_data(var_str, float_type) - \
                                  self.open_data(var_str, float_type, self.filter_width, 'LF') - \
                                  self.open_data(var_str, float_type, self.filter_width, 'HF', 'Rot')

            elif not self.filter_width:
                if self.helm_domain.lower() == 'rot':
                    dat = self.open_data(var_str, float_type, helm='Rot')
                elif self.helm_domain.lower() == 'div':
                    if self.var_exist(var_str, helm='Div'):
                        dat = self.open_data(var_str, float_type, helm='Div')
                    else:
                        dat = self.open_data(var_str, float_type) - self.open_data(var_str, float_type, helm='Rot')
                    if v != 0:
                        dat = dat - v

        elif not self.helm_domain:
            if self.freq_domain == 'LF' or self.freq_domain == 'BP':
                dat = self.open_data(var_str, float_type, self.filter_width, self.freq_domain)
                if v != 0:
                    dat = dat - v

            elif self.freq_domain == 'HF':
                if self.var_exist(var_str, self.filter_width, 'HF'):
                    dat = self.open_data(var_str, float_type, self.filter_width, 'HF')
                else:
                    dat = self.open_data(var_str, float_type) - self.open_data(var_str, float_type, self.filter_width, 'LF')

        sys.stdout.flush()
        if np.any(np.isnan(dat)):
            print('NO DATA: ', var_str)
            return np.NaN
        if self.reduce_mean and var_str != 'dux' and var_str != 'dvx':
            dat_mean = np.mean(dat, axis=(0, 2), keepdims=True)
            print('Reducing mean for: ', var_str, np.mean(dat_mean))
            return dat - dat_mean
        else:
            return dat


class TwoDataSets(object):
    def __init__(self, filt_width1, freq_domain1, helm_domain1, filt_width2, freq_domain2, helm_domain2,
                 coor_shift, exp, depth, reduce_mean, ind_min=0, ind_max=722):
        self.data_set1 = DataSet(filter_width=filt_width1, freq_domain=freq_domain1, helm=helm_domain1,
                                 exp=exp, depth=depth, coor_shift=coor_shift, reduce_mean=reduce_mean, ind_min=ind_min,
                                 ind_max=ind_max)
        self.data_set2 = DataSet(filter_width=filt_width2, freq_domain=freq_domain2, helm=helm_domain2,
                                 exp=exp, depth=depth, coor_shift=coor_shift, reduce_mean=reduce_mean, ind_min=ind_min,
                                 ind_max=ind_max)
        SimulationParameter.__init__(self, exp)

    def return_data(self, var_str, flag, float_type=64):
        if flag == 1:
            dat = self.data_set1.return_data(var_str, float_type=float_type)
        elif flag == 2:
            dat = self.data_set2.return_data(var_str, float_type=float_type)
        elif flag == 0:
            tmp = self.data_set2.original
            self.data_set2.original = True
            dat = self.data_set2.return_data(var_str, float_type=float_type)
            self.data_set2.original = tmp
        elif flag == 3:
            dat = self.data_set1.return_data(var_str, float_type=float_type) + \
                  self.data_set2.return_data(var_str, float_type=float_type)
        return dat


def interpolation(data, _ind):
    # this interpolation is for unexisted line, for adding a linesssssssss, not for an line that is wrong
    width = 10
    f2 = interp1d(np.delete(np.arange(width + 1), int(width / 2)), data[_ind - int(width / 2):_ind + int(width / 2)],
                  axis=0, kind='cubic')
    return np.concatenate((data[:_ind, :, :], f2(int(width / 2))[np.newaxis, :, :], data[_ind:, :, :]), 0)
