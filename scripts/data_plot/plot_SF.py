##
from open_data import *
from data_file_names import *
from tools.plot_tools import *
from scripts.plot_data.create_Flxs_dict import create_Flx_dict, concatenate_dictionaries
from scipy.signal import savgol_filter

# run parameters

normalize_work = 0
no_factor=True
def SF_plot(figures, title, temporal, depths=None, smoothing=False, freq=None, Flxs=None, fig=None, ax=None,
            testing=False, flux_remainder=False): # veertical or galilean

    Flxs, freq = flx_dictionary(Flxs, depths, figures, freq, temporal, testing=testing, flux_remainder=flux_remainder)
    if fig is None:
        fig, axes = plot_initializing(figures)
    fig.suptitle(title, fontsize=sizes.title_size, y=0.995, x=0.515)
    freq_len = 362 if temporal else 725
    if flux_remainder=='CG' or flux_remainder=='CG_ray':
        freq_len = 512

    # data initialization
    for m in figures.keys():
        if ax is None:
            ax_tmp = fig.axes[m]
        else:
            ax_tmp = ax
        ci = 0
        handles=[]

        for direction in figures[m]['directions']:
            for dat_head, flags_list in zip(figures[m]['data_headers'], figures[m]['flags']):

                param = SimulationParameter(dat_head.exp)
                flx_tot = np.zeros(freq_len)
                print(dat_head.exp, dat_head.filter, dat_head.parameters)

                for plot_flags in flags_list:

                    if len(np.shape(plot_flags)) == 1:
                        plot_flags = (plot_flags,)

                    if plot_flags == ((0, 0, 0),):
                        filter_tmp = (0, None, None, 0, None, None)
                    else:
                        filter_tmp = dat_head.filter

                    flx_tmp = np.zeros(freq_len)
                    for flags in plot_flags:
                        if direction == 'Total':
                            flx_tmp += np.nanmean(Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags]['Vertical'], axis=0)
                            flx_tmp += np.nanmean(Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags]['Horizontal'], axis=0)
                        else:
                            flx_tmp += np.nanmean(Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags][direction], axis=0)

                    if flags != (0, 0, 0):
                        flx_tot += flx_tmp

                    if normalize_work:
                        if dat_head.exp == 'Stochastic':
                            flx_tmp /= 23.4 / (2 * 10 ** 9)
                        elif dat_head.exp == 'Steady':
                            flx_tmp /= 20.4 / (2 * 10 ** 9)

                    if flux_remainder=='CG_ray':
                        index=np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 112, 144, 176, 240, 368, 496])-1
                        flx_tmp=np.concatenate((np.zeros(1), np.flip(flx_tmp[index])))

                    if smoothing and not np.any(np.isnan(flx_tmp)):
                        if temporal:
                            flx_tmp = np.concatenate((np.array([flx_tmp[0]]), savgol_filter(flx_tmp[1:], 5, 3)))
                            flx_tmp = np.concatenate((np.array([flx_tmp[0]]), savgol_filter(flx_tmp[1:], 15, 3)))
                            # flx_tmp = np.concatenate((np.array([flx_tmp[0]]), savgol_filter(flx_tmp[1:], 21, 3)[0:10],
                            #                           savgol_filter(flx_tmp[1:], 21, 3)[10:]))
                        elif flux_remainder=='CG_ray':
                            flx_tmp = np.concatenate((flx_tmp[0:6], savgol_filter(flx_tmp, 7, 3)[6:17], flx_tmp[17:19]))
                            flx_tmp = np.concatenate((savgol_filter(flx_tmp, 5, 3)[0:6], savgol_filter(flx_tmp, 5, 2)[6:17],
                                                      savgol_filter(flx_tmp, 5, 3)[17:19]))
                        else:
                            flx_tmp = np.concatenate((np.array([flx_tmp[0]]), savgol_filter(flx_tmp[1:], 15, 3)))

                    if temporal:
                        freq_tmp = freq / (F_cor / 2 / np.pi)
                    else:
                        freq_tmp = freq * param.L
                    if flux_remainder=='CG_ray':
                        freq_tmp=np.concatenate((np.zeros(1), np.flip(freq_tmp[index])/390))

                    _label, color, linestyle, linewidth, zorder, outline = plotting_parameters(ci, figures, m)
                    ci += 1

                    if temporal:
                        xlabel = 'freq/f'.replace('', '\hspace{0.04em}')
                        if normalize_work:
                            ylabel = r'\Pi/\tau_{work}\cdot cph'
                        else:
                            if no_factor:
                                ylabel = r'\Pi'
                            else:
                                ylabel = r'\Pi \cdot freq/f'
                    else:
                        xlabel = 'k_{h}\hspace{0.04em}L'
                        if normalize_work:
                            ylabel = r'\Pi/\tau_{work}\cdot k_{h}L'
                        else:
                            if no_factor:
                                ylabel = r'\Pi'
                            else:
                                ylabel = r'\Pi \cdot k_{h}L'

                    if flux_remainder=='CG_ray':
                        alpha=0.8
                    else:
                        alpha=1

                    if no_factor:
                        pass
                    else:
                        flx_tmp = flx_tmp * freq_tmp

                    ylim=None
                    # if flux_remainder=='CG_ray':
                    #     ax.scatter(freq_tmp, flx_tmp, color=color, marker='o', zorder=zorder)
                    #     handle = ax.lines[-1]
                    # else:
                    handle = plot_1d(freq_tmp, flx_tmp, ax=ax_tmp, title=figures[m]['title'], xlabel=xlabel,
                                     ylabel=ylabel,
                                     xscale='log', yscale='linear', ylim=ylim, color=color, linestyle=linestyle,
                                     linewidth=linewidth, label=_label, zorder=zorder, outline=outline, alpha=alpha)

                    handles.append(handle)

        # # need to figure out what to do with legends:
        # if np.any(outline==True):
        #     ax.legend(handles=handles, labels=figures[m]['labels'], handlelength=1.6)
        # else:
        #     ax.legend(loc='best', fontsize=sizes.legend_size)

    return freq_tmp, Flxs, fig


def SF_3d_plot(figures, title, temporal, fig=None, ax=None, depths=None, freq=None, Flxs=None, vmin=None, vmax=None,
               cmap=None):

    Flxs, freq = flx_dictionary(Flxs, depths, figures, freq, temporal)
    if fig is None:
        fig, axes = plot_initializing(figures)
    fig.suptitle(title, fontsize=sizes.title_size, y=0.995, x=0.515)
    freq_len = 362 if temporal else 725
    depths_len=129

    # data initialization
    for m in figures.keys():
        if ax is None:
            ax = fig.axes[m]
        for direction in figures[m]['directions']:
            for dat_head, flags_list in zip(figures[m]['data_headers'], figures[m]['flags']):

                param = SimulationParameter(dat_head.exp)
                flx_tot = np.zeros((depths_len, freq_len))

                for flags in flags_list:

                    if flags == ((0, 0, 0),):
                        filter_tmp = (0, None, None, 0, None, None)
                    else:
                        filter_tmp = dat_head.filter

                    flx_tmp = np.zeros((depths_len, freq_len))
                    if direction == 'Total':
                        flx_tmp += Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags]['Vertical']
                        flx_tmp += Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags]['Horizontal']
                    else:
                        flx_tmp += Flxs[dat_head.exp][filter_tmp][dat_head.parameters][flags][direction]

                    if flags != (0, 0, 0):
                        flx_tot += flx_tmp

                if temporal:
                    # freq_tmp = freq / (F_cor / 2 /np.pi)
                    freq_tmp = freq * 3600
                else:
                    freq_tmp = freq * param.L
                flx_tmp = freq_tmp * flx_tmp

                if temporal:
                    # xlabel = '$f/f_{cor}$'
                    xlabel = 'cph'
                    if normalize_work:
                        ylabel = r'\Pi/\tau_{work}\cdot cph'
                    else:
                        ylabel = r'\Pi \cdot cph'
                else:
                    xlabel = '$k_{h}L$'
                    ylabel = '$depth [km]$'
                    if normalize_work:
                        zlabel = r'\Pi/\tau_{work}\cdot k_{h}L'
                    else:
                        zlabel = r'$\Pi \cdot k_{h}L$'
                        # ylabel = '${(cm/s)}^2$'

                plot_3d_log(freq_tmp, np.linspace(0,2,129), flx_tmp, ax=ax, fig=None, title=figures[m]['title'],
                            xlabel=xlabel, ylabel=ylabel, colorbar_label=zlabel, vmin=vmin, vmax=vmax, cmap=cmap,
                            log_axes=['x'])

    return freq_tmp, flx_tmp, fig


def plotting_parameters(ci, figures, m):
    # plotting parameters
    _label = figures[m]['labels'][ci]
    color = figures[m]['color'][ci]
    linestyle = figures[m]['style'][ci]
    linewidth = figures[m]['size'][ci]
    zorder = figures[m]['zorder'][ci]
    if 'outline' in figures[m].keys():
        outline = figures[m]['outline'][ci]
    else:
        outline = False
    return _label, color, linestyle, linewidth, zorder, outline


def flx_dictionary(Flxs, depths, figures, freq, temporal, testing=False, flux_remainder=False):
    if Flxs is None or freq is None:
        Flxs = {}
        for m in figures.keys():
            calc_directions = ['Vertical', 'Horizontal'] if figures[m]['directions'] == ['Total'] else figures[m][
                'directions']
            for data_header, flag_list in zip(figures[m]['data_headers'], figures[m]['flags']):
                Flxs_tmp, freq = create_Flx_dict(flag_list, calc_directions, data_header, temporal=temporal,
                                                 x_direction=False, depths=depths, testing=testing,
                                                 flux_remainder=flux_remainder)
                Flxs = concatenate_dictionaries(Flxs_tmp, Flxs)
    return Flxs, freq


def plot_initializing(figures):
    # plotting initialization
    fig = plt.figure()
    fig.subplots_adjust(left=0.08, right=0.95)
    fig.subplots_adjust(hspace=0.1, wspace=0.2)

    if len(figures.keys()) == 1:
        axes = fig.subplots(nrows=1, ncols=1)
        fig.set_figheight(6)
        fig.set_figwidth(7)
    elif len(figures.keys()) == 2:
        axes = fig.subplots(nrows=1, ncols=2)  # sharey=True
        fig.subplots_adjust(left=0.08, right=0.95)
        fig.subplots_adjust(bottom=0.1, top=0.86, hspace=0.1, wspace=0.2)
        fig.set_figheight(5)
        fig.set_figwidth(10.5)
    elif len(figures.keys()) == 3:
        axes = fig.subplots(nrows=1, ncols=3)  # sharey=True
        fig.subplots_adjust(left=0.08, right=0.95)
        fig.subplots_adjust(bottom=0.1, top=0.86, hspace=0.1, wspace=0.2)
        fig.set_figheight(5)
        fig.set_figwidth(16)
    elif len(figures.keys()) == 4:
        axes = fig.subplots(nrows=2, ncols=2)
        fig.set_figheight(15)
        fig.set_figwidth(15)
    return fig, axes

# if axes_ratio == 'equal':
#     ylim=[0,0]
#     for i in range(len(figures.keys())):
#         ylim_min, ylim_max = fig.axes[i].get_ylim()
#         ylim[1] = np.max((ylim[1], np.max(ylim_max)))
#         ylim[0] = np.min((ylim[0], np.min(ylim_min)))
#     # factor = np.max(np.abs(ylim))*0.1
#     # ylim = [ylim[0] - factor, ylim[1] + factor]
#     for i in range(len(figures.keys())):
#         fig.axes[i].set_ylim(ylim[0], ylim[1])
#         fig.axes[i].set_aspect(1.0 / fig.axes[i].get_data_ratio(), adjustable='box')
