# sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
from data_file_names import get_depths_run, print_start
from open_data import *
from tools.flow_analysis_tools import *

# data processing
experiments=['Stochastic']

filter_widths = [48,    48,    48,    24,   24,    24,   24,    16,   16,    16,   16,  16]
freq_domains =  ['LF', 'LF',  'HF',  'LF',  'LF', 'LF', 'HF',  'LF', 'HF',  'LF', 'HF', 'HF']
helm_domains =  ['rot', None, None, 'rot', 'div', None, None, 'rot', 'rot', None, None, 'div']
# filter_widths = [ 16]
# freq_domains =  ['HF']
# helm_domains =  ['div']
geostrophic = False
coor_shift = True
# # filter_widths = [11, 18, (11, 18)]
# # freq_domains=['LF', 'LF', 'BP']s
# filter_widths = [(24, 16)] *3ss
# freq_domains=['BP'] *3
filter_order = 3

if len(sys.argv)>3:
    experiments = [str(sys.argv[3])]


def filter_and_save(data_set, var_str, filter_width, freq_domain, file_name):
    if initialize_data(file_name, var_str, description='filtered velocities'):
        original_data = data_set.return_data(var_str, float_type=64)
        filtered_data = butter_sos2_filter(original_data, filter_width=filter_width, dt=1)
        if freq_domain == 'HF':
            filtered_data = original_data - filtered_data
        finalize_data(file_name, var_str, filtered_data)
    # change_variable(file_name, var_str, filtered_data)
    sys.stdout.flush()


def filter_new_and_save(data_set, var_str, file_name, filter_order):
    file_name = file_name.replace('Filtered', 'Filtered_order_%d' % filter_order)
    if initialize_data(file_name, var_str, description='filtered velocities'):
        data = butter_sos2_filter(data_set.return_data(var_str, float_type=64), filter_width=filter_width, dt=1,
                                  filter_order=filter_order)
        finalize_data(file_name, var_str, data)
        sys.stdout.flush()


depths=get_depths_run(sys.argv)
for depth in np.flip(depths):
    for filter_width, freq_domain, helm_domain in zip(filter_widths, freq_domains, helm_domains):
        for exp in experiments:

            data_set = DataSet(exp=exp, depth=depth, helm=helm_domain, geostrophic=geostrophic, coor_shift=coor_shift)
            print(data_set.find_self_file_name('u'))
            print_start(depth, sys.argv)
            if exp=='Stochastic' and depth==77:
                continue
            file_name = get_data_file_name(exp, filter_width=filter_width, freq_domain=freq_domain, helm=helm_domain,
                                           geostrophic=geostrophic, coor_shift=coor_shift, dbz_source=False).format(f'{depth:03d}')
            # file_name = file_name.replace('Helmholtz', 'Helmholtz_test')
            print(file_name)

            if filter_order==3:
                filter_and_save(data_set, 'u', filter_width, freq_domain, file_name)
                filter_and_save(data_set, 'v', filter_width, freq_domain, file_name)
                if data_set.helm_domain != 'rot':
                    filter_and_save(data_set, 'w', filter_width, freq_domain, file_name)
                filter_and_save(data_set, 'dux', filter_width, freq_domain, file_name)
                filter_and_save(data_set, 'dvx', filter_width, freq_domain, file_name)
                # # filter_and_save(data_set, 'dwx', filter_width, freq_domain, file_name)
                filter_and_save(data_set, 'duy', filter_width, freq_domain, file_name)
                filter_and_save(data_set, 'dvy', filter_width, freq_domain, file_name)
                # # filter_and_save(data_set, 'dwy', filter_width, freq_domain, file_name)
                if data_set.helm_domain != 'rot':
                    filter_and_save(data_set, 'duz', filter_width, freq_domain, file_name)
                    filter_and_save(data_set, 'dvz', filter_width, freq_domain, file_name)
                    filter_and_save(data_set, 'dwz', filter_width, freq_domain, file_name)
                # filter_and_save(data_set, 'dbx', filter_width, freq_domain, file_name)
                # filter_and_save(data_set, 'dby', filter_width, freq_domain, file_name)
                # filter_and_save(data_set, 'dbz', filter_width, freq_domain, file_name)

            sys.stdout.flush()

print('BYE BYE')
sys.stdout.flush()
