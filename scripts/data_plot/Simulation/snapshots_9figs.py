from tools.plot_tools import *
import matplotlib.colors as mcolors

"""
final plot
"""
# data parameters
experiments = ['Steady', 'Steady', 'Steady', 'Stochastic', 'Stochastic', 'Stochastic', 'Stochastic', 'Stochastic', 'Stochastic']
freq_domains = [None, 'LF', 'LF'] * 2 + [None, 'LF', 'HF']
filter_widths = [0, 24, 24] * 2 + [0, 24, 16]
helm_domains = [None] * 9
geostrophics = [False] * 9
reduce_means = [False] * 9
coor_shifts = [True] * 9
dims = [2] * 9
depth_orientation = False
depths = [128] * 3 + [128] * 3 + [65] * 3
_time = 14
depth_name = {128:'Surface', 126:'Surface', 119:'150m', 65:'1000m'}

data_headers = list(map(DataHeader, experiments, filter_widths, freq_domains, helm_domains, coor_shifts, geostrophics,
                        dims, reduce_means))

# data parameters

# plotting parameters
sizes.title_size = 12
sizes.label_size = 10
sizes.legend_size = 8
sizes.tick_size = 14
sizes.title_weight = 1000
sizes.colorbar_size = 10

colors = [(0,0,0.6),(0,0,1),'deepskyblue',(1,1,1),'yellow',(1,0,0),(0.6,0,0)]
_cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

# Title='Ro = $\\frac{\zeta}{f_{cor}}$'
Title=''
figures = {}
for i in range(9):
    figures[i] = {'title': i, 'filters': data_headers[i]}

# figure
fig = plt.figure()
a=20
d=1
b=2*d
c=3*a+4*b
e1=4
e2 = 1
f=3*a+e1 +e2
gs0 = fig.add_gridspec(f, c)

ax1 = fig.add_subplot(gs0[0:a, 0:a])
ax2 = fig.add_subplot(gs0[0:a, a:2*a])
ax3 = fig.add_subplot(gs0[0:a, 2*a+2*b:3*a+2*b])
cax1 = fig.add_subplot(gs0[0:2*a+d, 2*a:2*a+b])

ax4 = fig.add_subplot(gs0[a+d:2*a+d, 0:a])
ax5 = fig.add_subplot(gs0[a+d:2*a+d, a:2*a])
ax6 = fig.add_subplot(gs0[a+d:2*a+d, 2*a+2*b:3*a+2*b])
cax2 = fig.add_subplot(gs0[0:2*a+d, 3*a+2*b:3*a+3*b])

ax7 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 0:a])
ax8 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, a:2*a])
ax9 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 2*a+2*b:3*a+2*b])
cax3 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 3*a+2*b:3*a+3*b])

# ax7 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 0:a+d])
# ax8 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, a+d:2*a+2*d])
# ax9 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 2*a+2*d:3*a+3*d])
# cax3 = fig.add_subplot(gs0[2*a+d+e1:3*a+d+e1+e2, 3*a+2*b:3*a+3*b])


fig.set_figheight(10)
fig.set_figwidth(11.5)
gs0.update(bottom=0.07, left=0.1, right=0.95, top=0.93)
gs0.update(wspace=8, hspace=8)

axes=[ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
caxes=[cax1, cax2, cax3]

# fig.subplots_adjust(left=0.08, right=0.95)
# fig.subplots_adjust(hspace=0.2, wspace=0.05)
# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.99, x=0.5)
ticks = []
yticklabel = []
shape = (512, 513)

cmaps=[]
for m in figures.keys():
    title = figures[m]['title']
    ax = axes[m]
    ci = 0
    dat_head = figures[m]['filters']
    depth = depths[m]

    data_set = DataSet(filter_width=dat_head.filter_width, freq_domain=dat_head.freq_domain, helm=dat_head.helm,
                       geostrophic=dat_head.geostrophic, dim=dat_head.dim, reduce_mean=dat_head.reduce_mean,
                       exp=dat_head.exp, depth=depth, coor_shift=dat_head.coor_shift)
    data_set.ind_min = _time
    data_set.ind_max = _time + 1

    x = np.arange(data_set.SHAPE[2]) * data_set.DX
    y = np.arange(data_set.SHAPE[1]) * data_set.DX

    if m in [0,1,3,4]:
        try:
            # data = np.empty((513,512))
            data = np.squeeze(data_set.return_data('Vort'))/F_cor
        except:
            data = (np.squeeze(data_set.return_data('dvx'))-np.squeeze(data_set.return_data('duy')))/F_cor
        vmax = 1.5
        _cmap_tmp = _cmap
    elif m in [2,5]:
        # data = np.empty((513,512))
        data = (np.squeeze(data_set.return_data('dux')) + np.squeeze(data_set.return_data('dvy'))) / F_cor
        if dat_head.exp=='Stochastic':
            vmax= 0.5
        elif dat_head.exp == 'Steady':
            vmax= 0.5
        # vmax = np.nanmax(np.abs(data)) * 0.4
        _cmap_tmp = _cmap
    else:
        # data = np.empty((513,512))
        data = np.squeeze(data_set.return_data('w'))
        vmax = 0.04
        _cmap_tmp = 'Greys'
    vmin = -vmax

    #
    # data_set = {}
    # data_set['SHAPE']=(11,10)
    # data_set['DX']=1
    # shape=data_set['SHAPE']
    #
    # x = np.arange(data_set['SHAPE'][1]) * data_set['DX']
    # y = np.arange(data_set['SHAPE'][0]) * data_set['DX']
    # data = np.random.random(shape)
    # if m in [0,1,3,4]:
    #     vmax = 1.5
    #     _cmap_tmp = _cmap
    # elif m in [2,5]:
    #     vmax= 0.5
    #     _cmap_tmp = _cmap
    # else:
    #     vmax = 0.04
    #     _cmap_tmp = 'Greys'
    # vmin = -vmax

    _, cmap = plot_3d(data, ax=ax, fig=fig, plot_colorbar=False, hor_ver='hor', vmin=vmin, vmax=vmax,
                        title=title, _cmap=_cmap_tmp, ticks=ticks, grid={}, yticklabels=yticklabel, xvline=0,
                      shape=shape)
    cmaps.append(cmap)


cbar = fig.colorbar(cmaps[0], cax=caxes[0], pad=0.07, fraction=0.1)
caxes[0].tick_params(labelsize=sizes.tick_size-2)
caxes[0].set_title(r'$\mathbf{Ro}$', y=1.02, fontsize=sizes.tick_size+3)

cbar = fig.colorbar(cmaps[5], cax=caxes[1])
caxes[1].tick_params(labelsize=sizes.tick_size-2)
caxes[1].set_title(r'$\mathbf{\nabla_h \cdot u /f}$', y=1.02, fontsize=sizes.tick_size+3)

cbar = fig.colorbar(cmaps[8], cax=caxes[2])
caxes[2].tick_params(labelsize=sizes.tick_size-2)
caxes[2].set_title(r'$\mathbf{w[m/s]}$', y=1.03, fontsize=sizes.tick_size+3)


# row titles
fig.text(x=0.01, y=0.78, s=r'$\mathbf{%s}$' % 'LFW', rotation=90, fontsize=sizes.title_size, fontweight=sizes.title_weight)
fig.text(x=0.01, y=0.5, s=r'$\mathbf{%s}$' % 'COMB', rotation=90, fontsize=sizes.title_size, fontweight=sizes.title_weight)
fig.text(x=0.01, y=0.17, s=r'$\mathbf{%s}$' % 'COMB', rotation=90, fontsize=sizes.title_size, fontweight=sizes.title_weight)
# fig.text(x=0.42, y=0.33, s=r'$\mathbf{Vertical\ Velocity}$', fontsize=sizes.title_size+3, fontweight=sizes.title_weight)
# fig.text(x=0.32, y=0.95, s=r'$\mathbf{Vorticity}$', fontsize=sizes.title_size+3, fontweight=sizes.title_weight)
# fig.text(x=0.73, y=0.95, s=r'$\mathbf{Divergence}$', fontsize=sizes.title_size+3, fontweight=sizes.title_weight)
for i,letter in enumerate('abcdefghi'):
    axes[i].text(x=8, y=177, s=r'$\mathbf{(%s)}$' % letter, fontsize=sizes.title_size+3, fontweight=sizes.title_weight)
    # axes[i].text(x=8, y=180, s=r'$(%s)$' % letter, fontsize=sizes.title_size+4, fontweight=sizes.title_weight, family='Helvetica')

# column titles
for i in range(3):
    filter_str = filter_to_text(filter_widths[i], freq_domains[i], helm_domains[i], fw=False)
    axes[i].set_title(filter_str, fontsize=sizes.title_size+2, y=1.01, x=0.5)

for i in range(3,6):
    axes[i].set_title('')

for i in range(6,9):
    filter_str = filter_to_text(filter_widths[i], freq_domains[i], helm_domains[i], fw=False)
    axes[i].set_title(filter_str, fontsize=sizes.title_size+2, y=1.0, x=0.5)

for i in [0,1,2,3,4,5]:
    axes[i].set(xlabel=None)
    axes[i].set_xticklabels([])
for i in [1,2,4,5,7,8]:
    axes[i].set(ylabel=None)
    axes[i].set_yticklabels([])

plt.show()
save_fig(fig, title='Solutions_decomposition_128_65', format='jpg')


