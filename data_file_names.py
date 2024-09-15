import glob
import os
import re
import sys
from datetime import datetime
import numpy as np
from tools.simulation_analysis import get_coordinate_shift

data_format_stochastic = {
    'hor': '/atlantic/rbarkan/flow_solve/LS_hres_stochastic_final_hfoutput/slice_dir/2D_xyt/XY_hor_{}.nc',
    'ver': '/atlantic/rbarkan/flow_solve/LS_hres_stochastic_final_hfoutput/slice_dir/2D_xyt/XY_ver_{}.nc',
    'dbz': '/atlantic/rbarkan/flow_solve/LS_hres_stochastic_final_hfoutput/slice_dir/2D_xyt/N2_{}.nc'}
data_format_steady = {'hor': '/atlantic/rbarkan/flow_solve/LS_hres_steady_hfoutput/slice_dir/2D_xyt/XY_hor_{}.nc',
                      'ver': '/atlantic/rbarkan/flow_solve/LS_hres_steady_hfoutput/slice_dir/2D_xyt/XY_ver__{}.nc',
                      'dbz': '/atlantic/rbarkan/flow_solve/LS_hres_steady_hfoutput/slice_dir/2D_xyt/N2_{}.nc'}
data_format_stochastic_only = '/atlantic/rbarkan/flow_solve/Stochastic_only_New/slice_dir/2D/XY_{}.nc'
data_format_3D_stochastic = '/atlantic/rbarkan/flow_solve/LS_hres_stochastic_final_hfoutput/3D/XYZ_{}.nc'
data_format_3D_steady = '/atlantic/rbarkan/flow_solve//LS_hres_steady_hfoutput/slice_dir/3D/XYZ_{}.nc'
Ertel_PV_format = {'hor': '/atlantic/michalshaham/Data/Ertel_PV{}{}/{}/{}/XY_hor_{}.nc',
                   'ver': '/atlantic/michalshaham/Data/Ertel_PV{}{}/{}/{}/XY_ver_{}.nc'}

if os.uname()[1]=='Michals-MacBook-Pro.local' or os.uname()[1]=='Michals-MBP':
    Data_folder = '/Users/michal/Data/'
    PSD_format = '/Users/michal/Data/PSD{}{}/{}{}{}/{}/{}.nc'
    SF_format_decomposition = {'hor': '/Users/michal/Data/Spectral_Flux{}{}/{}/{}/XY_hor_{}.nc',
                               'ver': '/Users/michal/Data/Spectral_Flux{}{}/{}/{}/XY_ver_{}.nc'}
    CG_format_decomposition = '/Users/michal/Data/Coarse_Graining{}{}/{}/{}/XY_{}.nc'
    data_format_original = '/Users/michal/Data/Filtered/{}/original/XY_{}.nc'
    data_format_decomposition = '/Users/michal/Data/{}{}{}/{}/{}/XY_{}.nc'
    data_analysis_format = '/Users/michal/Data/Analysis{}{}{}/{}/{}/XY_{}.nc'

elif os.uname()[1]=='compute-0-15' or os.uname()[1]=='compute-0-16':
    Data_folder = '/atlantic/michalshaham/Data/'
    PSD_format = '/atlantic/michalshaham/Data/PSD{}{}/{}{}{}/{}/{}.nc'
    SF_format_decomposition = {'hor': '/atlantic/michalshaham/Data/Spectral_Flux{}{}/{}/{}/XY_hor_{}.nc',
                               'ver': '/atlantic/michalshaham/Data/Spectral_Flux{}{}/{}/{}/XY_ver_{}.nc'}
    CG_format_decomposition = '/atlantic/michalshaham/Data/Coarse_Graining{}{}/{}/{}/XY_{}.nc'
    data_format_original = '/atlantic/michalshaham/Data/Filtered/{}/original/XY_{}.nc'
    data_format_decomposition = '/atlantic/michalshaham/Data/{}{}{}/{}/{}/XY_{}.nc'
    data_analysis_format = '/atlantic/michalshaham/Data/Analysis{}{}{}/{}/{}/XY_{}.nc'

else:
# elif os.uname()[1]=='atlantic.tau.ac.il' or os.uname()[1]=='southern' or os.uname()[1]=='indian':
    Data_folder = '/atlantic/michalshaham/Data/'
    PSD_format = '/atlantic/michalshaham/Data/PSD{}{}/{}{}{}/{}/{}.nc'
    SF_format_decomposition = {'hor': '/atlantic/michalshaham/Data/Spectral_Flux{}{}/{}/{}/XY_hor_{}.nc',
                               'ver': '/atlantic/michalshaham/Data/Spectral_Flux{}{}/{}/{}/XY_ver_{}.nc'}
    CG_format_decomposition = '/atlantic/michalshaham/Data/Coarse_Graining{}{}/{}/{}/XY_{}.nc'
    data_format_original = '/atlantic/michalshaham/Data/Filtered/{}/original/XY_{}.nc'
    data_format_decomposition = '/atlantic/michalshaham/Data/{}{}{}/{}/{}/XY_{}.nc'
    data_analysis_format = '/atlantic/michalshaham/Data/Analysis{}{}{}/{}/{}/XY_{}.nc'

N2_format = {'time_avg': '/atlantic/michalshaham/Data/N2/{}/N2_time_avg.nc',
             'vol_avg': '/atlantic/michalshaham/Data/N2/{}/N2_volume_avg.nc'}

exp_name = {'Stochastic': 'COMB', 'Stochastic_only': 'HFW', 'Steady': 'LFW'}
helm_name = {'div': 'Divergent', 'rot': 'Rotational'}


def get_original_data_name(exp, direction, depth=None):
    depth = {} if depth is None else depth
    if exp == 'Steady':
        return data_format_steady[direction[0:3].lower()].format(depth)
    elif exp == 'Stochastic':
        return data_format_stochastic[direction[0:3].lower()].format(depth)


def decomp_str(filt, domain, helm):
    _str = ''
    if filt:
        if domain == 'BP':
            str_tmp = 'fw%d_fw%d_BP' % (filt[0], filt[1])
        else:
            str_tmp = 'fw%d_%s' % (filt, domain)
        _str += str_tmp
    if helm:
        if _str:
            _str += '_'
        _str += helm.lower()
    if not filt and not helm:
        _str += 'original'
    return _str


def get_SF_file_name(flags, exp, direction, filt1=None, filt2=None, domain1=None, domain2=None, helm1=None, helm2=None,
                     coor_shift=False, reduce_mean=False, temporal=False, flux_type=False):

    _str = ''
    if 1 in flags or 3 in flags:
        _str += decomp_str(filt1, domain1, helm1)
    if 2 in flags or 3 in flags:
        str2 = decomp_str(filt2, domain2, helm2)
        if _str and str2:
            _str += '_'
        _str += str2
    if flags == (0, 0, 0):
        _str = 'original'

    flg_str = flags_to_str(flags)
    if flg_str == 'Eee' or flg_str == 'Www':
        flg_str = 'Aaa'

    if coor_shift:
        coor_shift = '/Coor_shift'
    else:
        coor_shift = ''

    reduce_mean = '/Reduce_mean' if reduce_mean else ''
    file_name = SF_format_decomposition[direction[0:3].lower()].format(reduce_mean, coor_shift, exp, _str, flg_str)
    if temporal:
        file_name = file_name.replace('Spectral_Flux', 'Temporal_Flux')
    if flux_type== 'vertical':
        file_name = file_name.replace('Spectral_Flux', 'Vertical_Flux_remainder')
    elif flux_type== 'galilean':
        file_name = file_name.replace('Spectral_Flux', 'Flux_remainder')
    elif flux_type == 'CG':
        file_name = file_name.replace('Spectral_Flux', 'Coarse_graining')
    elif flux_type == 'CG_ray':
        file_name = file_name.replace('Spectral_Flux', 'Coarse_graining_ray')
    elif flux_type == 'enstrophy':
        file_name = file_name.replace('Spectral_Flux', 'Enstrophy')
    if not os.path.exists(file_name):
        if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
            file_name = file_name.replace('atlantic', 'analysis')
        elif os.uname()[1] == 'southern':
            file_name = file_name.replace('atlantic', 'southern')
        else:
            file_name = file_name.replace('atlantic', 'atlantic2')
        if not os.path.exists(file_name):
            print('FILE DOES NOT EXIST %s' % file_name)

    # file_name = file_name.replace('/Coor_shift', '/Coor_shift_reduce')

    return file_name


def get_CG_file_name(depth, exp, filt1=None, filt2=None, domain1=None, domain2=None, helm1=None, helm2=None,
                     coor_shift=False, reduce_mean=False, temporal=False, flux_remainder=False):

    _str = decomp_str(filt1, domain1, helm1) + '_' + decomp_str(filt2, domain2, helm2)

    if coor_shift:
        coor_shift = '/Coor_shift'
    else:
        coor_shift = ''

    reduce_mean = '/Reduce_mean' if reduce_mean else ''
    file_name = CG_format_decomposition.format(reduce_mean, coor_shift, exp, _str, depth)
    if temporal:
        file_name = file_name.replace('Spectral_Flux', 'Temporal_Flux_TEST_2')
    if flux_remainder=='vertical':
        file_name = file_name.replace('Spectral_Flux', 'Vertical_Flux_remainder')
    elif flux_remainder == 'CG_ray':
        file_name = file_name.replace('Coarse_graining', 'Coarse_graining_ray')
        file_name = file_name.replace('XY', 'Ray/XY')
    if not os.path.exists(file_name):
        if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
            file_name = file_name.replace('atlantic', 'analysis')
        elif os.uname()[1] == 'southern':
            file_name = file_name.replace('atlantic', 'southern')
        elif os.uname()[1] == 'meddy':
            file_name = file_name.replace('atlantic', 'meddy')
        else:
            file_name = file_name.replace('atlantic', 'atlantic2')
        if not os.path.exists(file_name):
            print('FILE DOES NOT EXIST %s' % file_name)

    # file_name = file_name.replace('/Coor_shift', '/Coor_shift_reduce')

    return file_name



def get_Ertel_PV_file_name(exp, direction, filt=None, domain=None, helm=None, geostrophic=False, dim=2,
                           reduce_mean=False):
    _str = decomp_str(filt, domain, helm)

    reduce_mean = '/Reduce_mean' if reduce_mean else ''
    geostrophic = '/Geostrophic' if geostrophic else ''
    file_name = Ertel_PV_format[direction[0:3].lower()].format(reduce_mean, geostrophic, exp, str, '{}')
    if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
        file_name = file_name.replace('atlantic', 'analysis')
    else:
        file_name = file_name.replace('atlantic', 'atlantic2')
    if dim == 3:
        file_name = file_name.replace('Data', 'Data3D')
    return file_name


def get_analysis_file_name(exp, filt=None, domain=None, helm=None, geostrophic=False, dim=2, reduce_mean=False,
                           coor_shift=False):
    U = get_coordinate_shift(exp, helm) if coor_shift else 0
    _str = decomp_str(filt, domain, helm)
    reduce_mean = '/Reduce_mean' if reduce_mean else ''
    geostrophic = '/Geostrophic' if geostrophic else ''
    coor_shift = '/Coor_shift' if coor_shift else ''
    file_name = data_analysis_format.format(reduce_mean, geostrophic, coor_shift, exp, _str, '{}')
    if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
        file_name = file_name.replace('atlantic', 'analysis')
    else:
        file_name = file_name.replace('atlantic', 'atlantic2')
    if dim == 3:
        file_name = file_name.replace('Data', 'Data3D')
        file_name = file_name.replace('XY', 'XYZ')

    return file_name


def get_data_file_name(exp, filter_width=0, freq_domain=None, helm=False, geostrophic=False, dbz_source=False,
                       direction=None, dim=2, filter_order=3, roy_folder=True, coor_shift=False):
    if (not filter_width) and (not helm) and (not geostrophic) and (not coor_shift) and roy_folder:
        if dim == 3:
            if exp == 'Stochastic':
                return data_format_3D_stochastic
            elif exp == 'Steady':
                return data_format_3D_steady
        else:
            if exp == 'Stochastic':
                return data_format_stochastic[direction[0:3].lower()]
            elif exp == 'Steady':
                return data_format_steady[direction[0:3].lower()]
            elif exp == 'Stochastic_only':
                return data_format_stochastic_only
    else:
        _str = decomp_str(filter_width, freq_domain, helm)
        if helm:
            decomp_type = 'Helmholtz'
        elif dbz_source:
            decomp_type = 'Dbz_N2'
        elif filter_width:
            decomp_type = 'Filtered'
            if filter_order != 3:
                decomp_type = decomp_type.replace('Filtered', 'Filtered_order_%d' % filter_order)
        else:
            decomp_type = 'Original'

        geostrophic = '/Geostrophic' if geostrophic else ''
        coor_shift = '/Coor_shift' if coor_shift else ''
        file_name = data_format_decomposition.format(decomp_type, geostrophic, coor_shift, exp, _str, '{}')
        if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
            file_name = file_name.replace('atlantic', 'analysis')
        else:
            file_name = file_name.replace('atlantic', 'atlantic2')
        if dim == 3:
            file_name = file_name.replace('Data', 'Data3D')
            file_name = file_name.replace('XY', 'XYZ')
        return file_name


def get_PSD_file_name(data_header, hours_cut=0):
    _str = decomp_str(data_header.filter_width, data_header.freq_domain, data_header.helm)

    reduce_mean = '/Reduce_mean' if data_header.reduce_mean else ''
    geostrophic = '/Geostrophic' if data_header.geostrophic else ''
    U = get_coordinate_shift(data_header.exp) if data_header.coor_shift else 0
    # coor_shift = '/Coor_shift'
    coor_shift = '/Coor_shift' if U != 0 else ''

    dim = '%dD' % data_header.dim if type(data_header.dim) == int else data_header.dim # dim can be 'dbz_N2' and 'Ro'
    cut_hr = ('/cut_hr%d' % hours_cut) if hours_cut else ''

    file_name = PSD_format.format(reduce_mean, geostrophic, dim, coor_shift, cut_hr, data_header.exp, _str)
    if dim == 3:
        file_name = file_name.replace('Data', 'Data3D')
    if os.uname()[1] == 'compute-0-15' or os.uname()[1] == 'compute-0-16':
        possible_machine = ['/atlantic/', '/atlantic2/']
    else:
        possible_machine = ['/atlantic/', '/atlantic2/']
    for machine in possible_machine:
        file_name = file_name.replace(machine,'{}')
    print(file_name)
    for machine in possible_machine:
        file_name_tmp = file_name.format(machine)
        print(file_name_tmp)
        if os.path.exists(file_name_tmp):
            return file_name_tmp
    else:
        print('FILE DOES NOT EXIST %s' % file_name_tmp)
    return file_name_tmp


def flags_to_str(flags, ab='EW'):
    flg_str = ''
    for i in range(len(flags)):
        if flags[i] == 1:
            flg_str += ab[0].capitalize() if i == 0 else ab[0].lower()
        elif flags[i] == 2:
            flg_str += ab[1].capitalize() if i == 0 else ab[1].lower()
        elif flags[i] == 0:
            flg_str += 'O' if i == 0 else 'o'
        elif flags[i] == 3:
            flg_str += 'S' if i == 0 else 's'
    return flg_str


def filter_to_text(filter_width, freq_domain, helm, coor_shift=False, exp=None, dim=None, geostrophic=False, fw=True):

    if filter_width and helm:
        if freq_domain == 'BP':
            _str = '{}hr-{}hr/{}'.format(filter_width[0], filter_width[1], helm)
        else:
            if fw:
                _str = '%s(%dhr) %s' % (freq_domain, filter_width, helm_name[helm.lower()][:3])
            else:
                _str = '%s %s' % (freq_domain, helm_name[helm.lower()][:3])
    elif filter_width:
        if freq_domain == 'BP':
            _str = '{}hr-{}hr'.format(filter_width[0], filter_width[1])
        else:
            if fw:
                _str = '%s(%dhr)' % (freq_domain, filter_width)
            else:
                _str = '%s' % (freq_domain)
    elif helm:
        _str = '%s' % helm_name[helm.lower()][:3]
    else:
        _str = 'no filter'

    if geostrophic:
        _str = _str + ' Geo.'

    if dim == 'dbz_N2':
        _str = r'$b_z/N^2$ ' + _str
    elif dim == 'Ro':
        _str = r'$\zeta/f_{cor}$ ' + _str

    if coor_shift:
        _str += ' [coor trans.]'
    if exp:
        _str += ' ' + exp_name[exp]

    return _str


def flow_decomp_to_title(filter, ab='EW', fw=True):
    filter_width1, filter_width2, freq_domain1, freq_domain2, helm1, helm2 = filter

    title1 = '%s = ' % ab[0].capitalize()
    title2 = '%s = ' % ab[1].capitalize()

    title1 += filter_to_text(filter_width1, freq_domain1, helm1, fw=fw)
    title2 += filter_to_text(filter_width2, freq_domain2, helm2, fw=fw)

    return title1, title2


def get_filters_run(argv):
    exp_name = str(argv[3])
    filter_width1 = int(argv[4])
    filter_width2 = int(argv[5])
    freq_domain1 = str(argv[6])
    freq_domain2 = str(argv[7])
    helm_domain1 = str(argv[8])
    helm_domain2 = str(argv[9])
    coor_shift = str(argv[10])
    if helm_domain1 == 'None':
        helm_domain1 = None
    if helm_domain2 == 'None':
        helm_domain2 = None
    if coor_shift == 'True':
        coor_shift = True
    elif coor_shift == 'False':
        coor_shift = False
    return exp_name, filter_width1, filter_width2, freq_domain1, freq_domain2, helm_domain1, helm_domain2, coor_shift


def get_depths_run(argv, exp=None, dim=2):
    if exp == 'Stochastic_only':
        depths = [65, 119, 128]
    else:
        if len(argv) > 1:
            min_depth = int(argv[1])
            max_depth = int(argv[2])
            depths = range(min_depth, max_depth + 1)
        else:
            if dim == 3:
                depths = np.arange(722)
            else:
                depths = get_num_from_file_path(get_original_data_name('Steady', 'hor').format('*[0-9][0-9][0-9]'))
    print('START RUNNING: DEPTHS {}-{}'.format(np.min(depths), np.max(depths)))
    sys.stdout.flush()
    return depths


def print_start(depth, argv):
    if len(argv) > 1:
        print('\nopen files depth: {} [{}-{}]'.format(depth, argv[1], argv[2]))
    else:
        print('\nopen files depth: {}'.format(depth))
    print(datetime.now(), ', process pid: ', os.getpid())
    sys.stdout.flush()


def get_depths_from_netcdf(dat, test=False):
    depths = []
    for depth in dat.variables.keys():
        if 'freq' not in depth and depth!='n' and depth!='time':
            depths.append(depth)
    depths = np.sort(list(map(int, depths)))
    if test:
        missing_depths = [d_ for d_ in range(1, 130) if d_ not in depths]
        existing_depths = [d_ for d_ in range(1, 130) if d_ in depths]
        if len(missing_depths) > 0 and missing_depths[0] != 77:
            print('Missing Depths: ', missing_depths)
            print('Existing Depths: ', existing_depths)
        else:
            print('All depths')
    return depths


def get_numbers_from_filename(filename):
    filename = filename.split('/')[-1]
    return re.findall(r'\d+', filename)[-1]


def get_num_from_file_path(file_type):
    nums = []
    for filename in glob.glob(file_type):
        nums.append(int(get_numbers_from_filename(filename)))
    nums = sorted(nums)
    return nums
