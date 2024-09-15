##
from tools.plot_tools import *
from scipy.signal import savgol_filter

# data gathering objects
psd_dict = {}
freq_dict = {}

##

def PSD_plot(figures, Title, hours_cut, depth_orientation, c, line_s, s, smoothing=False, plot_f=True, fig=None, ax=None):

    # plotting initialization
    # row, col = get_row_col(len(figures.keys()))
    # fig = plt.figure()
    if fig is None:
        if len(figures.keys())==1:
            fig, axes = plt.subplots(nrows=1, ncols=1)
            fig.set_figheight(6)
        elif len(figures.keys()) == 2:
            fig, axes = plt.subplots(nrows=1, ncols=2)
            fig.set_figheight(5)
            fig.set_figwidth(13)
        elif len(figures.keys()) == 4:
            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig.set_figheight(15)
            fig.set_figwidth(15)
    # if you want to change x=0.5 - do it inside a plotting script
    fig.suptitle(Title, fontsize=sizes.title_size, fontweight=sizes.title_weight, y=0.97, x=0.5)

    # data initialization
    for m in figures.keys():
        title = figures[m]['title']
        # ax = plt.subplot(row * 100 + col * 10 + m + 1)
        if ax is None:
            ax = fig.axes[m]
        ci = 0
        for dat_head in figures[m]['filters']:
            file_name = get_PSD_file_name(dat_head, hours_cut=hours_cut)
            # file_name = file_name.replace('PSD', 'PSD/times_30_180')
            dat = Dataset(file_name)
            print(dat.filepath())
            sys.stdout.flush()
            # determine depths to plot
            depths = get_depths_from_netcdf(dat, True)
            depths = np.delete(depths, np.where(depths==129))

            # determine dimensions of psd
            if dat_head.dim == 1 or dat_head.dim in ['dbz_N2', 'Ro']:
                psd = np.zeros((len(depths), dat.variables['freq_t'].size))
            elif dat_head.dim == 2:
                psd = np.zeros((len(depths), dat.variables['freq_h'].size))
            elif dat_head.dim in [3, 4]:
                psd = np.zeros((len(depths), dat.variables['freq_t'].size, dat.variables['freq_h'].size))
            psd.fill(np.nan)

            # fill psd matrix. this is time for factors if needed!
            for ind, depth in zip(range(len(depths)), np.flip(depths)):
                if dat_head.dim in [1, 2, 'dbz_N2', 'Ro']:
                    psd[ind, :] = dat.variables[f'{depth:03d}'][:]
                elif dat_head.dim in [3, 4]:
                    psd[ind, :, :] = dat.variables[f'{depth:03d}'][:]

            # get frequencies. this is time for factors if needed!
            freq=0
            if dat_head.dim in [1, 'dbz_N2', 'Ro']:
                freq_var = [key for key in dat.variables.keys() if 'freq' in key][0]
                freq = dat.variables[freq_var][:]
            elif dat_head.dim == 2:
                freq_var = [key for key in dat.variables.keys() if 'freq' in key][0]
                freq = dat.variables[freq_var][:] * 200
            elif dat_head.dim == 3:
                freq_t = dat.variables['freq_t'][:]
                freq_h = dat.variables['freq_h'][:] * 200
            elif dat_head.dim == 4:
                freq_t = dat.variables['freq_t'][:]
                freq_h = dat.variables['freq_h'][:] * 200

            if np.any(np.isnan(psd)):
                print('TESING: there is NAN in PSD matrix')

            # if dat_head.dim == 4:
            #     psd = np.nanmean(psd, axis=2) * 512 * 722
            #     print(np.shape(psd))
            psd_dict[dat_head.header] = psd
            if dat_head.dim == 1:
                # psd = np.nanmean(psd, axis=0)
                freq_dict[dat_head.header] = freq
                psd = np.nanmean(psd, axis=0)

            if smoothing:
                if dat_head.dim==1:
                    # psd = np.concatenate((np.array(psd[0:5]), savgol_filter(psd[5:50], 11, 1),
                    #                       savgol_filter(psd[50:200], 3, 1), savgol_filter(psd[200:], 11, 1)))
                    if dat_head.exp == 'Stochastic':
                        # psd = np.concatenate((np.array(psd[0:3]), savgol_filter(psd[3:23], 3, 1),
                        #                       savgol_filter(psd[23:84], 3, 2), savgol_filter(psd[84:], 5, 1)))
                        psd = np.concatenate((savgol_filter(psd, 3, 1)[:15], savgol_filter(psd, 3, 1)[15:75],
                                              savgol_filter(psd, 11, 1)[75:200], savgol_filter(psd, 21, 1)[200:]))
                    elif dat_head.exp == 'Steady':
                        psd = np.concatenate((savgol_filter(psd, 3, 1)[:4], savgol_filter(psd, 7, 1)[4:40],
                                              savgol_filter(psd, 11, 1)[40:200], savgol_filter(psd, 21, 1)[200:]))
                        # psd = savgol_filter(psd, 9, 1)
                    # psd = np.concatenate((np.array(psd[0:5]), savgol_filter(psd[5:], 3, 1)))

                elif dat_head.dim==2:
                    psd = np.nanmean(psd, axis=0)
                    psd = np.concatenate((np.array([psd[0]]), savgol_filter(psd[1:], 15, 1)))

            # plotting
            if dat_head.dim in [1, 2, 'dbz_N2', 'Ro']:

                # plotting parameters
                _label = figures[m]['labels'][ci]
                color = c[ci]
                linestyle = s[ci]
                linewidth = line_s[ci]
                ci += 1

                if dat_head.dim in [1, 4, 'dbz_N2', 'Ro']:
                    # xlabel = 'cph'.replace('', '\hspace{0.04em}')
                    # if dat_head.dim==1:
                    #     ylabel = 'PSD/cph'.replace('', '\hspace{0.04em}')
                    xlabel = 'freq/f'.replace('', '\hspace{0.04em}')
                    if dat_head.dim==1:
                        ylabel = 'PSD '.replace('', '\hspace{0.04em}')+'\hspace{0.6em}[m^2s^{-2}/s^{-1}]'
                        # ylabel = 'PSD'.replace('', '\hspace{0.04em}')
                    elif dat_head.dim in ['dbz_N2', 'Ro']:
                        ylabel = '{(m/s)}^2/k^2'
                elif dat_head.dim == 2:
                    xlabel = 'k_{h}\hspace{0.04em}L'
                    ylabel = 'PSD '.replace('', '\hspace{0.04em}')+'\hspace{0.6em}[m^2s^{-2}/km^{-1}]'
                    # ylabel = 'PSD'.replace('', '\hspace{0.04em}')

                if dat_head.dim == 1:
                    freq = freq / 3600 / (F_cor / 2 / np.pi )
                    psd = psd * 3600
                plot_1d(freq, psd, ax=ax, title=title, xlabel=xlabel, ylabel=ylabel, xscale='log', yscale='log',
                        color=color, linestyle=linestyle, linewidth=linewidth, label=_label)

            elif dat_head.dim == 3:
                if depth_orientation:
                    xlabel = '$k_{h}L$'
                    ylabel = '$freq/f$'
                    colorbar_label = '${(cm/s)}^2$'
                    depths = -np.arange(129) * 15 / 1000
                    im = plot_3d_log(freq, depths, np.mean(psd, 0), ax=ax, fig=fig, title=title, xlabel=xlabel,
                                     ylabel=ylabel, colorbar_label=colorbar_label, vmin=1e-8, vmax=1.9, cmap='jet',
                                     semi_log=True)
                else:
                    xlabel = '$freq/f$'
                    ylabel = '$depth [km]$'
                    colorbar_label = '${(cm/s)}^2$'
                    im = plot_3d_log(freq_h, freq_t, np.mean(psd, 0), ax=ax, fig=fig, title=title, xlabel=xlabel,
                                     ylabel=ylabel, colorbar_label=colorbar_label, vmax=1e4, vmin=1e-2, cmap='jet')
                    # fig.colorbar(im, ax=axes.ravel().tolist())

            elif dat_head.dim == 4:
                xlabel = '$k_{x}L$'
                ylabel = '$cph$'
                colorbar_label = '${(cm/s)}^2$'
                im = plot_3d_log(freq_h, freq_t, np.mean(psd, 0), ax=ax, fig=fig, title=title, xlabel=xlabel,
                                 ylabel=ylabel, colorbar_label=colorbar_label, vmax=1e-9, vmin=1*10**(-3.5), cmap='jet')
                # fig.colorbar(im, ax=axes.ravel().tolist())

        if dat_head.dim in [1, 'dbz_N2', 'Ro'] and plot_f:
            # ax.axvline(x=1, color='grey', linestyle='--', linewidth=1.5,
            #            label='$f\hspace{0.05em},  (16hr)^{-1}\hspace{0.05em},  (24hr)^{-1}$', zorder=1)
            ax.axvline(x=1, color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
            ax.axvline(x=1 / 16 / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
            ax.axvline(x=1 / 24 / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
            ax.legend()

    return psd_dict, freq_dict, fig

