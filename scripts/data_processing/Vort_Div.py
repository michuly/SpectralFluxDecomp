import sys

from data_file_names import get_depths_run, print_start

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/analysis/michalshaham/PythonProjects/Oceanography'])
from open_data import *
from tools.simulation_analysis import vorticity_netcdf


# data processing
to_save = True
experiments=['Stochastic']
# experiments=['Stochastic_only']
filter_widths = [0]
freq_domains = [None]
helm_domains = [None]
reduce_mean = False
geostrophic = False
dim = 2
coor_shift=False


def vort_div_and_save(data_set, file_name):

    var_str = 'Vort'
    if initialize_data(file_name, var_str=var_str, description='Analysis Parameters', exp=exp):
        print('Calculating Vort...')
        data = vorticity_netcdf(data_set, check_existing=False)
        finalize_data(file_name, var_str, data)
        sys.stdout.flush()

    # var_str = 'Div'
    # if initialize_data(file_name, var_str=var_str, description='Analysis Parameters', exp=exp):
    #     print('Calculating Div...')
    #     data = divergence_netcdf(data_set)
    #     finalize_data(file_name, var_str, data)
    #     sys.stdout.flush()


depths=get_depths_run(sys.argv, dim=dim)
for depth in np.flip(depths):
    for filter_width in filter_widths:
        for exp in experiments:
            for freq_domain in freq_domains:
                for helm_domain in helm_domains:
                    print_start(depth, sys.argv)
                    if exp=='Stochastic' and depth==77 and dim==2:
                        continue
                    data_set = DataSet(filter_width=filter_width, freq_domain=freq_domain, helm=helm_domain,
                                       geostrophic=geostrophic, dim=dim, reduce_mean=reduce_mean, exp=exp,
                                       depth=depth, coor_shift=coor_shift)
                    file_name = get_analysis_file_name(exp, filt=filter_width, domain=freq_domain,
                                helm=helm_domain, geostrophic=geostrophic, coor_shift=coor_shift,
                                                       reduce_mean=False).format('%03d' % depth)
                    if exp=='Stochastic_only':
                        roy_folder=True
                    else:
                        roy_folder=False
                    get_data_file_name(exp, filter_width, freq_domain, helm_domain, geostrophic, dim=dim,
                                                   roy_folder=roy_folder, coor_shift=coor_shift).format(f'{depth:03d}')

                    print(file_name)
                    # print(data_set.find_self_file_name(var_str[1]))
                    sys.stdout.flush()

                    vort_div_and_save(data_set, file_name=file_name)
                    sys.stdout.flush()

                if filter_width == 0:
                    break

print('BYE BYE')
sys.stdout.flush()



