# sys.path.extend(['/Users/michal/PythonProjects/Oceanography'])
from tools.plot_tools import *
import matplotlib.colors as mcolors
import sys


# data parameters
experiments = ['Steady', 'Stochastic', 'Steady', 'Stochastic']
freq_domains = [None] * 6
filter_widths = [0] * 6
helm_domains = [None] * 6
geostrophics = [False] * 6
reduce_means = [False] * 6
coor_shifts = [False, False, True, True]
dims = [2] * 6
_time = 500

data_headers = list(map(DataHeader, experiments, filter_widths, freq_domains, helm_domains, coor_shifts, geostrophics,
                        dims, reduce_means))

# plotting parameters
sizes.title_size = 18
sizes.label_size = 18
sizes.legend_size = 10
sizes.tick_size = 20
sizes.title_weight = 800
sizes.colorbar_size = 12

# Title='Ro = $\\frac{\zeta}{f_{cor}}$'
Title=''
figures = {}
figures[0] = {'title': '', 'filters': data_headers[0]}
figures[1] = {'title': '', 'filters': data_headers[1]}
figures[2] = {'title': '', 'filters': data_headers[2]}
figures[3] = {'title': '', 'filters': data_headers[3]}

colors = [(0,0,0.6),(0,0,1),'deepskyblue',(1,1,1),'yellow',(1,0,0),(0.6,0,0)]
_cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

fps=10
duration=int(722/fps)
# duration=4
depth = 128
datas = []

# for m in figures.keys():
#     dat_head = figures[m]['filters']
#     data_set = DataSet(filter_width=dat_head.filter_width, freq_domain=dat_head.freq_domain, helm=dat_head.helm,
#                        geostrophic=dat_head.geostrophic, dim=dat_head.dim, reduce_mean=dat_head.reduce_mean,
#                        exp=dat_head.exp, depth=depth, coor_shift=dat_head.coor_shift)
#     datas.append(data_set.return_data('Vort') / F_cor)

def plot_coor(t, min):

    fig, axes = plt.subplots(nrows=2, ncols=2)
    fig.set_figheight(10)
    fig.set_figwidth(10.5)
    fig.subplots_adjust(left=0.12, right=0.87)
    fig.subplots_adjust(hspace=0.14, wspace=0.1)
    cax = fig.add_axes([0.9, 0.11, 0.03, 0.77])

    for m in figures.keys():
        title = figures[m]['title']
        ax = fig.axes[m]

        vmin, vmax = -1.3, 1.3
        shape = (512, 513)
        ticks = [-1, 0, 1]
        yticklabel = []

        _, cmap = plot_3d(datas[m][t-min, :, :], ax=ax, fig=fig, plot_colorbar=False, hor_ver='hor', vmin=vmin,
                            vmax=vmax, title=title, _cmap=_cmap, ticks=ticks, grid={}, yticklabels=yticklabel,
                            xvline=0, shape=shape)
        fig.axes[m].tick_params(length=7, width=1.3, direction='inout', which='major')
        fig.axes[m].tick_params(length=5, width=1, direction='inout', which='minor')

    fig.axes[0].set(xlabel=None)
    fig.axes[1].set(xlabel=None)
    fig.axes[1].set(ylabel=None)
    fig.axes[3].set(ylabel=None)
    fig.axes[0].tick_params(labelbottom=False)
    fig.axes[1].tick_params(labelleft=False, labelbottom=False)
    fig.axes[3].tick_params(labelleft=False)

    cbar = fig.colorbar(cmap, cax=cax)
    cbar.ax.tick_params(labelsize=sizes.tick_size)
    cax.text(x=-0.3, y=1.4, s=r'$\boldsymbol{\zeta/f}$', fontsize=sizes.tick_size)

    fig.text(x=0.01, y=0.655, s='Before  shift', rotation=90,
             fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
    fig.text(x=0.01, y=0.24, s='After  shift', rotation=90,
             fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
    fig.axes[0].set_title(r'$\mathbf{%s}$' % 'LFW'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2,
                          y=1.1, x=0.5, family='Helvetica')
    fig.axes[1].set_title(r'$\mathbf{%s}$' % 'COMB'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2,
                          y=1.1, x=0.5, family='Helvetica')

    save_fig(fig, title='Ro_movie_2/Ro_coor_shift_movie_%d' % t, data_folder=True, dpi=400)

    return cmap


def make_frame(t):
    ind = int(t * fps)
    plot_colorbar = True if str(t) == '0' else False

    for m in figures.keys():
        fig, axes = plt.subplots(nrows=2, ncols=2)
        fig.set_figheight(10)
        fig.set_figwidth(10.5)
        fig.subplots_adjust(left=0.12, right=0.87)
        fig.subplots_adjust(hspace=0.14, wspace=0.1)
        cax = fig.add_axes([0.9, 0.11, 0.03, 0.77])

        title = figures[m]['title']
        ax = fig.axes[m]

        vmin, vmax = -1.3, 1.3
        shape = (512, 513)
        ticks = [-1, 0, 1]
        yticklabel = []

        _, cmap = plot_3d(datas[m][ind, :, :], ax=fig.axes[m], fig=fig, plot_colorbar=False, hor_ver='hor', vmin=vmin,
                            vmax=vmax, title=title, _cmap=_cmap, ticks=ticks, grid={}, yticklabels=yticklabel,
                            xvline=0, shape=shape)
    fig.axes[0].set(xlabel=None)
    fig.axes[1].set(xlabel=None)
    fig.axes[1].set(ylabel=None)
    fig.axes[3].set(ylabel=None)

    #
    if plot_colorbar:
        cbar = fig.colorbar(cmap, cax=cax)
        cbar.ax.tick_params(labelsize=sizes.tick_size)
        cax.text(x=-0.3, y=1.4, s=r'$\boldsymbol{\zeta/f}$', fontsize=sizes.tick_size)

        fig.text(x=0.01, y=0.68, s='Before  shift', rotation=90,
                 fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
        fig.text(x=0.01, y=0.24, s='After  shift', rotation=90,
                 fontsize=sizes.title_size, fontweight=sizes.title_weight, family='Helvetica')
        fig.axes[0].set_title(r'$\mathbf{%s}$' % 'LFW'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2,
                              y=1.1, x=0.5, family='Helvetica')
        fig.axes[1].set_title(r'$\mathbf{%s}$' % 'COMB'.replace('', '\hspace{0.05em}'), fontsize=sizes.title_size + 2,
                              y=1.1, x=0.5, family='Helvetica')

    return mplfig_to_npimage(fig)


# my_animation = VideoClip(make_frame, duration=duration)
# name = replace_all('Ro_coor_shift')
# my_animation.write_videofile(Data_folder + 'Figures/{}.mp4'.format(name), fps=fps)

if len(sys.argv)>1:
    min = int(sys.argv[1])
    max = int(sys.argv[2])
    print(min, max)


for m in figures.keys():
    dat_head = figures[m]['filters']
    data_set = DataSet(filter_width=dat_head.filter_width, freq_domain=dat_head.freq_domain, helm=dat_head.helm,
                       geostrophic=dat_head.geostrophic, dim=dat_head.dim, reduce_mean=dat_head.reduce_mean,
                       exp=dat_head.exp, depth=depth, coor_shift=dat_head.coor_shift)
    datas.append(data_set.return_data('Vort')[min:max,:,:] / F_cor)

for i in range(min, max):
    cmap = plot_coor(i, min)

# ffmpeg -r 12 -i Ro_coor_shift_movie_%d.png -vcodec mpeg4 -y /Users/michal/Documents/master_project/Paper/Figures/Ro_Movie.mp4
# if __name__ == '__main__':
#     main()

