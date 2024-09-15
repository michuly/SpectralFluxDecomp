import matplotlib.pyplot as plt

from tools.plot_tools import *
import matplotlib.colors as mcolors

"""
final plot
"""
# data parameters
depth_name = {128:'Surface', 126:'Surface', 119:'150m', 65:'1000m'}


# data parameters

# plotting parameters
sizes.title_size = 14
sizes.label_size = 10
sizes.legend_size = 8
sizes.tick_size = 14
sizes.colorbar_size = 10

# Title='Ro = $\\frac{\zeta}{f_{cor}}$'
Title=''

# figure
fig, axes = plt.subplots(nrows=2, ncols=2)
axes = axes.flatten()
fig.set_figheight(10)
fig.set_figwidth(10.5)
fig.subplots_adjust(left=0.12, right=0.93)
fig.subplots_adjust(hspace=0.18, wspace=0.2)
# cax = fig.add_axes([0.9, 0.11, 0.03, 0.4])

ticks = []
yticklabel = []
shape = (512, 513)

cmaps=[]
for m in range(4):
    ax = axes[m]
    ci = 0

    x = np.arange(512) * 390
    y = np.arange(513) * 390
    ind_t=444

    if m==0:
        dat = Dataset('/atlantic2/michalshaham/Data/Coarse_Graining/Coor_shift/Stochastic/fw24_LF_rot_fw16_HF/Ray/CG_local_2/XY_n10_128.nc')
        data = dat.variables['SPHeeE'][ind_t,:,:]
        title = '$\Pi^{Eee}$ $E=Rot-LP$ (srf)'
        vmin, vmax = -6e-7, 6e-7

    if m==1:
        dat = Dataset('/atlantic2/michalshaham/Data/Coarse_Graining/Coor_shift/Stochastic/fw24_LF_fw16_HF/Ray/CG_local_2/XY_n16_128.nc')
        data = dat.variables['SPHeeE'][ind_t,:,:]
        data = +dat.variables['SPVeeE'][ind_t,:,:]
        title = '$\Pi^{Eee}$ $E=LP$ (srf)'
        vmin, vmax = -4e-9, 4e-9

    if m==2:
        dat = Dataset('/atlantic2/michalshaham/Data/Coarse_Graining/Coor_shift/Stochastic/fw24_LF_fw16_HF/Ray/CG_local_2/XY_n32_65.nc')
        data = dat.variables['SPHwwE'][ind_t,:,:]
        data = +dat.variables['SPVwwE'][ind_t,:,:]
        title = '$\Pi^{Eww}$ (1km)'
        vmin, vmax = -6e-8, 6e-8

    if m==3:
        dat = Dataset('/atlantic2/michalshaham/Data/Coarse_Graining/Coor_shift/Stochastic/fw24_LF_fw16_HF/Ray/CG_local_2/XY_n32_65.nc')
        data = dat.variables['SPHweW'][ind_t,:,:]
        data = dat.variables['SPVweW'][ind_t,:,:]
        data = +dat.variables['SPHewW'][ind_t,:,:]
        data = +dat.variables['SPVewW'][ind_t,:,:]
        title = '$\Pi^{scatt}$ (1km)'
        vmin, vmax = -6e-8, 6e-8

    _, cmap = plot_3d(data, ax=ax, fig=fig, plot_colorbar=True, hor_ver='hor', vmin=vmin, vmax=vmax,
                        title=title, _cmap='seismic', ticks=ticks, grid={}, yticklabels=yticklabel, xvline=0,
                      shape=shape)

for i in range(2):
    axes[i].set(xlabel=None)
    axes[i].set_xticklabels([])
for i in [1, 3]:
    axes[i].set(ylabel=None)
    axes[i].set_yticklabels([])
# cbar = fig.colorbar(cmap, cax=cax)
# cbar.ax.tick_params(labelsize=sizes.tick_size)
# cax.text(x=-0.3, y=1.06e-8, s=r'$\Pi\ [W/kg]$', fontsize=sizes.tick_size)
for i,letter in enumerate('abcd'):
    axes[i].text(x=7, y=178, s=r'$(%s)$' % letter, fontsize=sizes.title_size+3, fontweight=sizes.title_weight)

plt.show()
save_fig(fig, title='CG_snapshot_interactions', format='jpg')


