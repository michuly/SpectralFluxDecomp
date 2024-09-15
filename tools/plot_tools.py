import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm, Normalize
import matplotlib.patheffects as pe
import matplotlib.animation as animation
from matplotlib.lines import Line2D
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
import numpy as np

from open_data import *
from tools.save_netCDF import create_netcdf
from matplotlib.ticker import FormatStrFormatter

# Show the major grid lines with dark grey lines
# ax.grid(b=True, which='major', color='#666666', linestyle='-')
# ax.minorticks_on()
# ax.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.6)
# https://stackoverflow.com/questions/14324477/bold-font-weight-for-latex-axes-label-in-matplotlib

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    'text.latex.preamble': r"\usepackage{amsmath} \DeclareMathSizes{11}{11}{7}{7} \DeclareMathSizes{12}{12}{8}{8} "
                           r"\DeclareMathSizes{16}{16}{14}{14} \DeclareMathSizes{16}{20}{12}{12}"
})



class font_size(object):
    def __init__(self):
        self.title_size = 16
        self.label_size = 11
        self.label_pad = 7
        self.legend_size = 12
        self.tick_size = 14
        self.title_weight = 700
        # self.label_weight = 1000
        self.colorbar_size = 12
        # https://www.overleaf.com/learn/latex/Font_sizes%2C_families%2C_and_styles#Reference_guide

sizes = font_size()


exp_name = {'Stochastic': 'COMB', 'Stochastic_only': 'HFW', 'Steady': 'LFW'}
helm_name = {'div': 'Divergent', 'rot': 'Rotational', None: 'no filter'}


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def get_row_col(graphs_num):
    if graphs_num == 4 or graphs_num == 3:
        row, col = 2, 2
    elif graphs_num == 2:
        row, col = 1, 2
    elif graphs_num == 1:
        row, col = 1, 1
    return row, col


def curve_fit(func, xdata, ydata, slices):
    popt, pcov = curve_fit(func, xdata[slices], ydata[slices])
    return func(xdata, *popt)[slices]


class DataHeader(object):
    def __init__(self, exp=None, filter_width=None, freq_domain=None, helm=None, coor_shift=False, geostrophic=False,
                 dim=2, reduce_mean=None):
        self.exp = exp
        self.filter_width = filter_width
        self.freq_domain = freq_domain
        self.helm = helm
        self.coor_shift = coor_shift
        self.geostrophic = geostrophic
        self.dim = dim
        self.reduce_mean = reduce_mean
        self.header = (exp, filter_width, freq_domain, helm, coor_shift, geostrophic, dim, reduce_mean)


class TwoDataHeader(object):
    def __init__(self, exp=None, filter_width1=None, freq_domain1=None, helm1=None, filter_width2=None,
                 freq_domain2=None, helm2=None, coor_shift=False, geostrophic=False, dim=2, reduce_mean=False):
        self.data_header1 = DataHeader(exp, filter_width1, freq_domain1, helm1, coor_shift, geostrophic,
                 dim, reduce_mean)
        self.data_header2 = DataHeader(exp, filter_width2, freq_domain2, helm2, coor_shift, geostrophic,
                 dim, reduce_mean)
        self.exp = exp
        self.coor_shift = coor_shift
        self.geostrophic = geostrophic
        self.dim = dim
        self.filter = (filter_width1, freq_domain1, helm1, filter_width2, freq_domain2, helm2)
        self.parameters = (coor_shift, geostrophic, reduce_mean)


def replace_all(text):
    dic = {'(': '', ')': '', ':': '_', ' ': '_', '$': '', 'depth': 'd', 'x_ind': 'x', ',': '', '=': '', '-': '_',
           '\\': '_'}
    for i, j in dic.items():
        text = text.replace(i, j)
    return text


def save_fig(fig, title, format='png', fps=10, data_folder=False, dpi=200):
    if os.uname()[1] == 'Michals-MacBook-Pro.local' or os.uname()[1]=='Michals-MBP' or os.uname()[1]=='kaplun510-5.tau.ac.il':
        # fig_file_name = '/Users/michal/Data/Figures/{}.{}'.format(replace_all(title), png_gif)
        if data_folder:
            fig_file_name = '/Users/michal/Data/Figures/{}.{}'.format(replace_all(title), format)
        else:
            fig_file_name = '/Users/michal/Documents/master_project/Paper/Figures/{}.{}'.format(replace_all(title),
                                                                                                format)
    elif os.uname()[1] == 'atlantic.tau.ac.il':
        fig_file_name = '/atlantic/michalshaham/Data/Figures/Paper_figs/{}.{}'.format(replace_all(title), format)
    print(fig_file_name)
    if format == 'png':
        fig.savefig(fig_file_name, bbox_inches='tight', dpi=dpi)
    elif format == 'jpg':
        fig.savefig(fig_file_name, bbox_inches='tight', dpi=dpi, format='jpg')
    elif format == 'gif':
        fig.write_gif(fig_file_name, fps=fps)


def get_time_slicer(t_start=300, duration=10, fps=10):
    Nt = fps * duration
    return slice(t_start, Nt * 2 + t_start, 2)


def cretae_subplots(cmap, Title=''):

    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_figheight(10)
    fig.set_figwidth(10.75)
    fig.subplots_adjust(left=0.12, right=0.88)
    fig.subplots_adjust(hspace=0.17, wspace=0.1)
    cax = fig.add_axes([0.9, 0.11, 0.03, 0.77])
    fig.suptitle(Title, fontsize=sizes.title_size + 2, y=0.99, x=0.5)

    for i in range(4):
        fig.axes[i].set_aspect(1.0, adjustable='box')
        fig.axes[i].set_xticks([0, 100, 200])
        fig.axes[i].set_yticks([0, 100, 200])
        fig.axes[i].tick_params(labelsize=sizes.tick_size)
    fig.text(x=0.01, y=0.66, s='NO trans.', rotation=90,
             fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
    fig.text(x=0.01, y=0.24, s='Coor. trans.', rotation=90,
             fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
    fig.axes[0].set_title(r'$\mathbf{%s}$' % 'LFW'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2, y=1.1,
                          x=0.5, family='Helvetica')
    fig.axes[1].set_title(r'$\mathbf{%s}$' % 'COMB'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2,
                          y=1.1, x=0.5, family='Helvetica')
    fig.axes[2].set_xlabel('x [km]', fontsize=sizes.label_size)
    fig.axes[3].set_xlabel('x [km]', fontsize=sizes.label_size)
    fig.axes[0].set_ylabel('y [km]', fontsize=sizes.label_size)
    fig.axes[2].set_ylabel('y [km]', fontsize=sizes.label_size)
    cbar = fig.colorbar(cmap, cax=cax)
    cbar.ax.tick_params(labelsize=sizes.tick_size)
    cax.text(x=-0.3, y=1.4, s=r'$\boldsymbol{\zeta/f}$', fontsize=sizes.tick_size)


def plot_1d(x, y, ax=None, title='', xlabel='', ylabel='', xscale='log', yscale='log', xlim=None, ylim=None,
            color=None, linestyle='-', linewidth=1.5, label='', axvline=False, zorder=2, outline=False, alpha=1):
    if not ax:
        fig, ax = plt.subplots()
    if outline:
        # path_effects = [pe.Stroke(linewidth=linewidth*2, foreground='k'), pe.Normal()]
        l1 = Line2D(x, y, linestyle="--", color=color, lw=linewidth+0.3)
        l1.set_dashes([0, 1.5 + 1, 3.5, 1])
        # l1.set_dashes([0, 2, 3.5, 2])
        l2 = Line2D(x, y, linestyle="--", color='grey', lw=linewidth+0.3)
        l2.set_dashes([1.5, 1, 0, 3.5 + 1])
        # l2.set_dashes([0, 1, 5.5, 1])
        ax.add_artist(l1)
        ax.add_artist(l2)
        handle=(l1, l2)
    else:
        # path_effects = None
        ax.plot(x, y, color=color, linestyle=linestyle, linewidth=linewidth, label=label, zorder=zorder, alpha=alpha)
        handle = ax.lines[-1]

    ax.set_title(title, fontsize=sizes.title_size, fontweight=sizes.title_weight, y=1.02, x=0.5)
    ax.set_xlabel(r'$\mathbf{%s}$' % xlabel, fontsize=sizes.label_size+2, labelpad=sizes.label_pad)
    ax.set_ylabel(r'$\mathbf{%s}$' % ylabel, fontsize=sizes.label_size+2, labelpad=sizes.label_pad)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(x[1], x[-2])
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.grid(True)
    # ax.legend(loc='best', fontsize=sizes.legend_size)
    ax.tick_params(labelsize=sizes.tick_size, width=1.5)
    tx = ax.yaxis.get_offset_text()
    tx.set_fontsize(sizes.tick_size)
    if axvline:
        ax.axvline(x=axvline, color='k', linestyle='--', label='f')
    return handle



def plot_3d_log(x, y, z, ax=None, fig=None, title='', xlabel='', ylabel='', colorbar_label='', vmin=None, vmax=None,
                cmap=None, log_axes=['z']):
    if cmap is None:
        cmap = 'jet'

    x_array = np.log10(x) if 'x' in log_axes else x
    y_array = np.log10(y) if 'y' in log_axes else y
    if 'z' in log_axes:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    X, Y = np.meshgrid(x_array, y_array)
    im = ax.pcolor(X, Y, z, norm=norm, cmap=cmap)
    # im = ax.pcolor(X, Y, z, norm=LogNorm(vmin=1e3, vmax=1e4), cmap=cmap)
    if 'x' in log_axes:
        ax.xaxis.set_major_formatter(FormatStrFormatter('$10^{%.1f}$'))
    if 'y' in log_axes:
        ax.yaxis.set_major_formatter(FormatStrFormatter('$10^{%.1f}$'))
    ax.set_title(title, fontsize=sizes.title_size)
    ax.set_xlabel(r'$\mathbf{%s}$' % xlabel, fontsize=sizes.label_size+2)
    ax.set_ylabel(r'$\mathbf{%s}$' % ylabel, fontsize=sizes.label_size+2)
    ax.set_ylim(top=np.max(y_array))
    ax.set_xlim(right=np.max(x_array))
    ax.tick_params(labelsize=sizes.tick_size, width=1.5)
    # ax.text(np.max(x_array) * 1.05, np.max(y_array) * 1.05, colorbar_label, fontsize=sizes.label_size)
    ax.grid()
    if fig:
        cbar = fig.colorbar(im, fraction=0.046, pad=0.04)
        cbar.ax.get_xaxis().labelpad = 6
        cbar.ax.set_xlabel(colorbar_label, fontsize=sizes.colorbar_size)
        cbar.ax.xaxis.set_label_position('top')

    return im


def plot_3d(data, ax=None, fig=None, plot_colorbar=True, hor_ver='hor', vmin=None, vmax=None, title='', _cmap=None,
            ticks=[], grid={}, yticklabels=[], xvline=0, shape=(512,513), colorbar_label=''):

    if _cmap is None:
        _cmap = 'viridis'

    if ax is None or fig is None:
        fig = plt.figure(0)
        ax = plt.subplot(111)

    if hor_ver == 'hor':
        # # using 201,202 to make the label 200 appear
        # y = np.arange(shape[0]) * 201 / shape[1]
        # x = np.arange(shape[1]) * 202 / shape[1]
        y = np.arange(shape[1]+1) * 200 / shape[1]
        x = np.arange(shape[0]+1) * 200 / shape[0]
        X, Y = np.meshgrid(x, y)
        cmap = ax.pcolor(X, Y, data, vmin=vmin, vmax=vmax, cmap=_cmap)
        # ax.plot(256*0.39, 333*0.39, '*k')
        # ax.plot(128*0.39, 256*0.39, '*k')
        # ax.plot(256*0.39, 128*0.39, '*k')
        # ax.plot(25*0.39, 384*0.39, '*k')
        ax.set_xticks([0,100,200])
        ax.set_yticks([0,100,200])
        ax.tick_params(labelsize=sizes.tick_size)
        ax.set_xlabel('x [km]', fontsize=sizes.label_size)
        ax.set_ylabel('y [km]', fontsize=sizes.label_size)
    elif hor_ver == 'ver':
        x = np.arange(shape[0]) * 15 / 1000
        y = np.arange(shape[1]) * 390 / 1000
        X, Y = np.meshgrid(x, y)
        cmap = ax.pcolor(np.transpose(Y), np.transpose(np.flip(X)), np.flip(data, axis=0), vmin=vmin, vmax=vmax,
                         cmap=_cmap)
        ax.set_xlabel('y [km]', fontsize=sizes.label_size)
        ax.set_ylabel('z [km]', fontsize=sizes.label_size)

    if xvline:
        ax.vlines(xvline, 0, np.max(y), linestyle='dashed')
    ax.set_aspect(1.0 / ax.get_data_ratio(), adjustable='box')
    ax.grid(grid)
    ax.minorticks_on()
    ax.set_title(title, fontsize=sizes.title_size, pad=10, loc='center')
    ax.locator_params(nbins=4)
    if plot_colorbar:
        # add if ticks?
        if ticks:
            cbar = fig.colorbar(cmap, ax=ax, fraction=0.046, pad=0.04, ticks=ticks)
        else:
            cbar = fig.colorbar(cmap, ax=ax, fraction=0.046, pad=0.04)

        if yticklabels:
            cbar.ax.set_yticklabels(yticklabels, fontsize=sizes.tick_size)
        else:
            cbar.ax.tick_params(labelsize=sizes.tick_size)
        # cbar.ax.get_xaxis().labelpad = 6
        # cbar.ax.set_xlabel(colorbar_label, fontsize=sizes.colorbar_size)
        # cbar.ax.xaxis.set_label_position('top')

    return fig, cmap


def plot_3d_subplots(data, axes=None, fig=None, plot_colorbar=True, hor_ver='hor', vmin=None, vmax=None, title='',
                     _cmap=None, ticks=[], grid={}, yticklabels=[], xvline=0, shape=(512,513)):
    if vmin is None:
        vmin = np.min(np.concatenate(data))
    if vmax is None:
        vmax = np.max(np.concatenate(data))
    if axes is None or fig is None:
        if len(data) == 2:
            fig, axes = plt.subplots(nrows=1, ncols=2)
            fig.set_figheight(5)
            fig.set_figwidth(13)
        elif len(data) == 4:
            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig.set_figheight(15)
            fig.set_figwidth(15)
            # fig.subplots_adjust(right=0.75, left=0.3, hspace=None, wspace=None)
        # fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [8, 8, 1]})
        # fig.set_figheight(5)
        # fig.set_figwidth(13)

    for _data, ax, _title, _hor_ver, _vmin, _vmax in zip(data, axes.flat, title[1:], hor_ver, vmin, vmax):
        _xvline = xvline if _hor_ver == 'hor' else None
        # plot_3d(_data, ax=ax, fig=fig, plot_colorbar=False, hor_ver=_hor_ver, vmin=vmin, vmax=vmax, title=_title,
        #         _cmap=_cmap, ticks=ticks, grid=grid, yticklabels=yticklabels, xvline=_xvline, shape=shape)
        plot_3d(_data, ax=ax, fig=fig, plot_colorbar=False, hor_ver=_hor_ver, vmin=_vmin, vmax=_vmax, title=_title,
                _cmap=_cmap, ticks=ticks, grid=grid, yticklabels=yticklabels, xvline=_xvline, shape=shape)
    # plot_3d(data[0], ax=ax[0], fig=fig, plot_colorbar=False, hor_ver=hor_ver, vmin=vmin, vmax=vmax, title=title[1],
    #         _cmap=_cmap, ticks=ticks, grid=grid, yticklabels=yticklabels, xvline=xvline)
    # plot_3d(data[1], ax=ax[1], fig=fig, plot_colorbar=False, hor_ver=hor_ver, vmin=vmin, vmax=vmax, title=title[2],
    #         _cmap=_cmap, ticks=ticks, grid=grid, yticklabels=yticklabels, xvline=xvline)
        if plot_colorbar:
            if _cmap is None:
                cmap = mpl.cm.viridis
            else:
                cmap = _cmap
            norm = mpl.colors.Normalize(vmin=_vmin, vmax=_vmax)
            fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
            # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes.ravel().tolist())
        # fig.subplots_adjust(right=0.8)
        # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax)
        #
        # fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax[2])
    fig.suptitle(title[0])
    return fig


def animating(data, vmin=None, vmax=None, title='', _cmap=None, ticks=[], grid={}, yticklabels=[],
              hor_ver='hor', duration=10, fps=10, xvline=0, time_stamp=True, shape=(512,513)):
    """
    try to reduce the number of arguments, use a dictionary
    """
    if type(data) == list:
        if len(data) == 2:
            fig, axes = plt.subplots(nrows=1, ncols=2)
            fig.set_figheight(5)
            fig.set_figwidth(13)
            # fig.subplots_adjust(hspace=0.3, wspace=0.2)
        elif len(data) == 4:
            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig.set_figheight(12)
            fig.set_figwidth(13)
            # fig.subplots_adjust(right=0.75, left=0.25, hspace=0.3, wspace=0.2)
        # fig, axes = plt.subplots(1, 3, gridspec_kw={'width_ratios': [8, 8, 1]})
        # fig.set_figheight(5)
        # fig.set_figwidth(13)
        # fig.subplots_adjust(hspace=0.3, wspace=0.2)
    else:
        fig, ax = plt.subplots()
        # ax.set_aspect(aspect=1)

    def make_frame(t):
        ind = int(t * fps)
        time_ind = ind if time_stamp else 0
        plot_colorbar = True if str(t) == '0' else False
        _slices = {}
        _slices['hor'] = (ind, slice(None), slice(None))
        _slices['ver'] = (slice(None), ind, slice(None))
        if type(data) == list:
            [ax.clear() for ax in axes.flat]
            datas = []
            for _data, _hor_ver in zip(data, hor_ver):
                # _data = organize_data(_data, _hor_ver, _slices, t, 0)
                datas.append(_data)
            plot_3d_subplots(datas, axes=axes, fig=fig, plot_colorbar=plot_colorbar,
                             hor_ver=hor_ver, vmin=vmin, vmax=vmax, title=title, _cmap=_cmap, ticks=ticks, grid=grid,
                             yticklabels=yticklabels, xvline=xvline, shape=shape)
        else:
            ax.clear()
#            _data = organize_data(data, hor_ver, _slices, t, 0)
            plot_3d(data, ax=ax, fig=fig, plot_colorbar=plot_colorbar, hor_ver=hor_ver, vmin=vmin, vmax=vmax,
                    title=title, _cmap=_cmap, ticks=ticks, grid=grid, yticklabels=yticklabels, xvline=xvline, shape=shape)
            ax.text(0.85, 0.85, str(ind), fontsize=16)

        return mplfig_to_npimage(fig)

    def organize_data(_data, _hor_ver, _slices, t, U):
        data_tmp = _data[_slices[_hor_ver]]
        _data = np.zeros(np.shape(data_tmp))
        x_0 = int(np.round(U * (t - 100)))  # U [ind/t]
        # x_0=0
        print(x_0)
        if x_0 > 0:
            _data[:, 0:(512 - x_0)] = data_tmp[:, x_0:512]
            _data[:, (512 - x_0):512] = data_tmp[:, 0:x_0]
        elif x_0 < 0:
            _data[:, 0:(-x_0)] = data_tmp[:, (512 + x_0):512]
            _data[:, -x_0:512] = data_tmp[:, 0:(512 + x_0)]
        elif x_0 == 0:
            _data = data_tmp
        return _data

    # my_animation = animation.FuncAnimation(fig, func=make_frame)
    my_animation = VideoClip(make_frame, duration=duration)
    title2 = title[0] if type(title) == list else title
    name = replace_all(title2)
    # animation.write_gif('/atlantic/michalshaham/Data/Figures/Figs/{}.gif'.format(name), fps=fps)
    my_animation.write_videofile(Data_folder + 'Figures/{}.mp4'.format(name), fps=fps)


def PDF_2D(x, y, z, title, x_label='', y_label='', bins=1000, vmax=None):
    plt.hist2d(x=x.flatten(), y=y.flatten(), weights=z.flatten(), bins=bins, vmax=vmax)
    plt.colorbar()

    plt.grid(b=True, which='major', color='#999999', linestyle='-')
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    # plt.locator_params(nbins=4)


def getting_data_from_data_set(var_str, filt_width, freq_domain, helm, geostrophic, exp, depth, time_slicer):
    data_set = DataSet(filter_width=filt_width, freq_domain=freq_domain, helm=helm, geostrophic=geostrophic, exp=exp,
                       depth=depth)
    data = data_set.return_data(var_str)[time_slicer, :, :]
    if np.shape(data)[0] == 1:
        data = data[0, :, :]
    return data


def building_vertical_matrix(func, var_str, exp, filter_width, freq_domain, helm, geostrophic, time_slicer, x_ind=220,
                             Nt=1,
                             saving_file_name=''):
    print('Building Vertical Matrix for variable: ', var_str)
    if Nt == 1:
        energy = np.zeros((129, 513))
    else:
        energy = np.zeros((129, Nt, 513))
    energy.fill(np.NaN)
    # depths = range(76, 90)
    depths = range(1, 130)
    for depth in depths:
        if depth in (77, 83, 85, 87):
            continue
        if Nt == 1:
            energy[depth - 1, :] = func(var_str=var_str, filt_width=filter_width, freq_domain=freq_domain,
                                        helm=helm, geostrophic=geostrophic, exp=exp, depth=depth,
                                        time_slicer=time_slicer)[:,
                                   x_ind]
        else:
            energy[depth - 1, :, :] = func(var_str=var_str, filt_width=filter_width, freq_domain=freq_domain,
                                           helm=helm, geostrophic=geostrophic, exp=exp, depth=depth,
                                           time_slicer=time_slicer)[:, :, x_ind]

    if saving_file_name:
        create_netcdf(name=saving_file_name, description=var_str, axes_names=('z', 't', 'y'),
                      shape=(129, 722, 513), dimensions=(2000, 200000, 200000), fields_keys=(var_str,),
                      fields_values=(energy,))
    return energy

