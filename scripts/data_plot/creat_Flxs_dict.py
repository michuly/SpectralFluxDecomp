from data_file_names import get_SF_file_name, get_depths_from_netcdf
from open_data import *

flags_list_all = [(1, 1, 1), (1, 2, 2), (1, 1, 2), (1, 2, 1), (2, 2, 1), (2, 1, 2), (2, 2, 2), (2, 1, 1)]


def concatenate_dictionaries(dic1, dic2):
    for exp in dic2:
        if exp not in dic1.keys():
            dic1[exp] = {}
        for filtering in dic2[exp].keys():
            if filtering not in dic1[exp].keys():
                dic1[exp][filtering] = {}
            for parameters in dic2[exp][filtering].keys():
                if parameters not in dic1[exp][filtering].keys():
                    dic1[exp][filtering][parameters] = {}
                for flags in dic2[exp][filtering][parameters].keys():
                    if flags not in dic1[exp][filtering][parameters].keys():
                        dic1[exp][filtering][parameters][flags] = {}
                    for direction in dic2[exp][filtering][parameters][flags].keys():
                        if direction not in dic1[exp][filtering][parameters][flags].keys():
                            dic1[exp][filtering][parameters][flags][direction] = \
                                dic2[exp][filtering][parameters][flags][direction]
    return dic1


def create_Flx_dict(flags_list, calc_directions, header, temporal=False, x_direction=False, depths=None, testing=False,
                    flux_remainder=False):
    if temporal:
        flx_len = None
    elif x_direction:
        flx_len = 513
    else:
        flx_len = 725
        # flx_len = 512
    if flux_remainder=='CG' or flux_remainder=='CG_ray':
        flx_len = 512

    if depths is None:
        depths_len = 129
        extract_depths = True
    else:
        depths_len = len(depths)
        extract_depths = False

    Flxs = {}
    Flxs = create_empty_dict(Flxs, calc_directions, depths_len, flags_list, flx_len, header)

    for exp in Flxs.keys():
        if temporal:
            if exp == 'Stochastic':
                flx_len = 361
            elif exp == 'Steady':
                flx_len = 362
        for filters in Flxs[exp].keys():
            for parameters in Flxs[exp][filters].keys():
                calculate_sum = False
                for flags in Flxs[exp][filters][parameters].keys():
                    for direction in Flxs[exp][filters][parameters][flags].keys():
                        coor_shift, geostrophic, reduce_mean = parameters
                        filter_width1, freq_domain1, helm1, filter_width2, freq_domain2, helm2 = filters
                        if testing:
                            freq=freq_for_fft(1024,390,1024,390)
                            x=np.arange(flx_len)
                            y = (np.random.rand() + 1) * np.sin((np.pi * 2 / 10) * freq * np.random.rand() + np.random.rand())
                            for depth in np.arange(depths_len):
                                Flxs[exp][filters][parameters][flags][direction][depth-1, :] = y
                            continue

                        if flags == (3, 3, 3): # check if there is a sum files
                            file_name = get_SF_file_name(flags, exp, direction, filt1=filter_width1,
                                                         filt2=filter_width2, domain1=freq_domain1,
                                                         domain2=freq_domain2, helm1=helm1, helm2=helm2,
                                                         coor_shift=coor_shift, reduce_mean=reduce_mean,
                                                         temporal=temporal, flux_remainder=flux_remainder)
                            if x_direction:
                                file_name = file_name.replace(exp, '%s/ax_x' % exp)

                            if not os.path.exists(file_name):
                                print('Sum file does not exist: ', file_name)
                                calculate_sum = True
                                # changing the dictionary after starting the loop! will it loop over everything?
                                Flxs = create_empty_dict(Flxs, calc_directions, depths_len, flags_list_all, flx_len,
                                                         header)
                                continue

                        file_name = get_SF_file_name(flags, exp, direction, filt1=filter_width1,
                                                     filt2=filter_width2, domain1=freq_domain1,
                                                     domain2=freq_domain2, helm1=helm1, helm2=helm2,
                                                     coor_shift=coor_shift, reduce_mean=reduce_mean,
                                                     temporal=temporal, flux_remainder=flux_remainder)

                        if x_direction:
                            file_name = file_name.replace(exp, '%s/ax_x' % exp)

                        is_w_rot1 = (helm1 is not None) and (helm1.lower() == 'rot') and (flags[1] == 1)
                        is_w_rot2 = (helm2 is not None) and (helm2.lower() == 'rot') and (flags[1] == 2)
                        if direction == 'Vertical' and (is_w_rot1 or is_w_rot2):  # Identity = 0
                            Flxs[exp][filters][parameters][flags][direction].fill(0)
                            continue

                        file_name = file_exist(file_name)
                        if not file_name:
                            print('File does not exist: ', file_name)
                            continue

                        dat = Dataset(file_name)
                        print(dat.filepath())
                        if extract_depths:
                            depths = get_depths_from_netcdf(dat, test=True)
                        dat.close()

                        for depth in depths:
                            dat = Dataset(file_name)
                            if temporal:
                                freq = dat.variables['freq_t'][:]
                            elif x_direction:
                                freq = dat.variables['freq_x'][:]
                            else:
                                freq = dat.variables['freq_h'][:]
                            if np.all(freq==0):
                                freq = freq_for_fft(1024,390,1024,390)
                            SF = dat.variables[f'{depth:03d}'][:]
                            dat.close()

                            if np.ndim(SF) == 2:
                                SF = np.mean(SF, 0)
                            # if np.isnan(SF[0]):
                            #     print('NAN - depth', depth, file_name)

                            depth_ind = depth-1
                            if flux_remainder in ['CG','CG_ray']:
                                Flxs[exp][filters][parameters][flags][direction][depth_ind, :]=SF
                            else:
                                Flx = np.zeros(flx_len)
                                for i in np.flip(np.arange(flx_len - 1)):
                                    Flx[i] = Flx[i + 1] - SF[i + 1]
                                Flxs[exp][filters][parameters][flags][direction][depth_ind, :] = Flx
                            # Flxs[exp][filters][parameters][flags][direction][depth_ind, :] = SF

                        print(exp, filters, parameters, flags, direction, np.sum(np.nanmean(SF,0)))

                        Flx_tmp = Flxs[exp][filters][parameters][flags][direction]
                        nan_list = np.where(np.isnan(Flx_tmp))[0]
                        if len(np.unique(nan_list)) == 1:
                            ind = nan_list[0]
                            width = 10
                            f2 = interp1d(np.delete(np.arange(width + 1), int(width / 2)),
                                          np.delete(Flx_tmp[(ind - int(width / 2)):(ind + int(width / 2)+1)],
                                                    int(width / 2), axis=0), axis=0, kind='cubic')
                            Flx_tmp = np.concatenate((Flx_tmp[:ind, :], f2(int(width / 2))[np.newaxis, :], Flx_tmp[ind+1:, :]), 0)
                            Flxs[exp][filters][parameters][flags][direction] = Flx_tmp

                # if saving_files:
                #     Flx = np.nanmean(Flxs[exp][filters][parameters][flags][direction], axis=0)
                #     file_name_for_plot = file_name.replace('Data', 'DataForPlots')
                #     create_netcdf(name=file_name_for_plot, description='', axes_names=('freq',), shape=(len(freq),),
                #                   dimensions=(np.max(freq),), fields_keys=('SF',), fields_values=(Flx,),
                #                   variable_axes=(), verbose=True)

                if calculate_sum:
                    flags = (3,3,3)
                    for flags_tmp in flags_list_all:
                        for direction in Flxs[exp][filters][parameters][flags].keys():
                            Flx = Flxs[exp][filters][parameters][flags_tmp][direction][:]
                            if np.all(np.isnan(Flxs[exp][filters][parameters][flags][direction][:])):
                                # if the entire depth is still NaN
                                Flxs[exp][filters][parameters][flags][direction][:] = Flx
                            else:
                                Flxs[exp][filters][parameters][flags][direction][:] += Flx

    return Flxs, freq


def create_empty_dict(Flxs, calc_directions, depths_len, flags_list, flx_len, header):
    flags_list = flags_list.copy()  # do not alter original list
    # for header in TwoDataHeaderList:
    if flx_len is None:
        if header.exp == 'Stochastic':
            flx_len = 361
        elif header.exp == 'Steady':
            flx_len = 362
    if header.exp not in Flxs.keys():
        Flxs[header.exp] = {}
    if header.filter not in Flxs[header.exp].keys():
        Flxs[header.exp][header.filter] = {}
    if header.parameters not in Flxs[header.exp][header.filter].keys():
        Flxs[header.exp][header.filter][header.parameters] = {}
    for flags_tmp in flags_list:
        if len(np.shape(flags_tmp)) == 1:
            flags_tmp = (flags_tmp,)
        if flags_tmp == ((0, 0, 0),):  # if flags=(0,0,0) than the filter must be original
            filters = (0, None, None, 0, None, None)
            if filters not in Flxs[header.exp].keys():
                Flxs[header.exp][filters] = {}
            if header.parameters not in Flxs[header.exp][filters].keys():
                Flxs[header.exp][filters][header.parameters] = {}
        else:
            filters = header.filter
        for flags in flags_tmp:
            if flags not in Flxs[header.exp][filters][header.parameters].keys():
                Flxs[header.exp][filters][header.parameters][flags] = {}
            for direction in calc_directions:
                if direction not in Flxs[header.exp][filters][header.parameters][flags].keys():
                    # Flxs[header.exp][filtering][header.parameters][flags][direction] = np.NaN
                    Flxs[header.exp][filters][header.parameters][flags][direction] = np.empty((depths_len, flx_len))
                    Flxs[header.exp][filters][header.parameters][flags][direction].fill(np.NaN)


    return Flxs
