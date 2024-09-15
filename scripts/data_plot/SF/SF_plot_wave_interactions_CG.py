from scripts.plot_data.plot_SF import *

"""
final plot
"""

##
# determine data plots
depths=None

# what to plot
direction = 'Total'
filter_widths1 = [24] * 6
filter_widths2 = [16] * 6
freq_domains1 = ['LF'] * 6  # LF
freq_domains2 = ['HF'] * 6 # HF
helms1 = [None, None, 'Rot']   # rot
helms1_cg = [None, None, None, None, 'Rot', 'Rot']   # rot
helms2 = [None] * 6 # divs
coor_shifts = [True] * 6
reduce_means = [False] * 6
exps = ['Stochastic'] * 6
ab = 'EW'
temporal = False

# plotting parameters
sizes.title_size = 18
sizes.label_size = 16
sizes.tick_size = 14

def get_label(data_headers, flags_list_plot):
    labels=[]
    for dat_head, flags_list in zip(data_headers, flags_list_plot):
        for plot_flags in flags_list:
            if plot_flags == (2,2,2):
                flags_str = '$\Pi^{Www}$'
            elif plot_flags == ((2,2,1),(2,1,2)):
                flags_str = '$\Pi^{scatt}$'
            elif plot_flags == (2,1,1):
                flags_str = '$\Pi^{Wee}$'
            elif plot_flags == (2,1,2):
                flags_str = '$\Pi^{Wew}$'
            elif plot_flags == (2,2,1):
                flags_str = '$\Pi^{Wwe}$'
            elif plot_flags == (1,2,2):
                flags_str = '$\Pi^{Eww}$'
            elif plot_flags == ((2,1,2),(2,2,1)):
                flags_str = '$\Pi^{scatt}$'
            dat_head1 = dat_head.data_header1
            filter1 = filter_to_text(dat_head1.filter_width, dat_head1.freq_domain, dat_head1.helm, fw=False)
            if flags_str!='$\Pi^{Www}$':
                flags_str = '%s (\hspace{0.1em}e = $%s$\hspace{0.1em})' % (flags_str, filter1.replace('-','{\\text -}'))

            labels.append(flags_str)
    return labels

### plotting parameters


data_headers1 = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
data_headers1_cg = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1_cg, filter_widths2, freq_domains2,
                         helms2, coor_shifts))

flags_list_plot1 = ([(2,2,2)], [(2,1,2)], [(2,1,2)])
flags_list_plot1_cg = ([(2,2,2)], [(2,1,2)], [(2,2,1)], [(1,2,2)], [(2,2,1)], [(1,2,2)])
# flags_list_plot1_cg = ([((2,1,2),(2,2,1))],)
labels1=get_label(data_headers1, flags_list_plot1)
labels1_cg=get_label(data_headers1_cg, flags_list_plot1_cg)
color1=['blue', 'darkgreen', 'goldenrod']
color1_cg=['blue', 'darkgreen', 'c', 'm', 'brown', 'red']

Title = 'Wave Interactions'.replace('', '\hspace{0.05em}') + ' $\Pi^{Www}, \Pi^{scatt}, \widetilde{\Pi}^{Eww}$'
size = [2] * 6
style1 = ['-'] * 6
style1_cg = [(0, (2,2))] * 6
zorder=[1] * 6


### initializing figure

fig, axis = plt.subplots(nrows=1, ncols=1)
fig.set_figheight(6)
fig.set_figwidth(6)
# fig.subplots_adjust(left=0.08, right=0.98, top=0.83)

figures = {}
figures[0] = {'title':'', 'flags': flags_list_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels1, 'color':color1, 'style':style1, 'size':size, 'zorder':zorder}
freq_sf, Flxs_sf, fig = SF_plot(figures, '', temporal, depths, fig=fig, ax=axis, testing=False)

figures[0] = {'title':'', 'flags': flags_list_plot1_cg, 'directions': [direction], 'data_headers': data_headers1_cg,
              'labels':labels1_cg, 'color':color1_cg, 'style':style1_cg, 'size':size, 'zorder':zorder}
freq_cg, Flxs_cg, fig = SF_plot(figures, '', temporal, depths, fig=fig, ax=axis, testing=False, flux_remainder='CG_ray')

fig.axes[0].legend(fontsize=14)

### Bottom panels

fig, axis = plt.subplots(nrows=1, ncols=1)
fig.set_figheight(6)
fig.set_figwidth(6)

def process_cg(freq_tmp, flx_tmp):
    index = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 112, 144, 176, 240, 368, 496]) - 1
    flx_tmp = np.nanmean(flx_tmp, axis=0)
    flx_tmp = np.concatenate((np.zeros(1), np.flip(flx_tmp[index])))
    flx_tmp = np.concatenate((flx_tmp[0:6], savgol_filter(flx_tmp, 7, 3)[6:17], flx_tmp[17:19]))
    flx_tmp = np.concatenate((savgol_filter(flx_tmp, 5, 3)[0:6], savgol_filter(flx_tmp, 5, 2)[6:17],
                              savgol_filter(flx_tmp, 5, 3)[17:19]))
    if no_factor:
        pass
    else:
        flx_tmp = flx_tmp * freq_tmp
    return freq_tmp, flx_tmp

def process_sf(freq_tmp, flx_tmp):
    flx_tmp = np.nanmean(flx_tmp, axis=0)
    flx_tmp = np.concatenate((np.array([flx_tmp[0]]), savgol_filter(flx_tmp[1:], 15, 3)))
    if no_factor:
        pass
    else:
        flx_tmp = flx_tmp * freq_tmp
    return freq_tmp, flx_tmp


exp='Stochastic'
sf_Www = Flxs_sf[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,2,2)]['Vertical'] + \
          Flxs_sf[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,2,2)]['Horizontal']
freq_sf, sf_Www = process_sf(freq_sf, sf_Www)

cg_Wwe = Flxs_cg[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,2,1)]['Vertical'] + \
          Flxs_cg[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,2,1)]['Horizontal']
freq_cg, cg_Wwe = process_cg(freq_cg, cg_Wwe)

cg_Wwe_rot = Flxs_cg[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(2,2,1)]['Vertical'] + \
          Flxs_cg[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(2,2,1)]['Horizontal']
freq_cg, cg_Wwe_rot = process_cg(freq_cg, cg_Wwe_rot)

sf_Wew = Flxs_sf[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,1,2)]['Vertical'] + \
          Flxs_sf[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(2,1,2)]['Horizontal']
freq_sf, sf_Wew = process_sf(freq_sf, sf_Wew)
cg_Wew = np.interp(freq_cg, freq_sf, sf_Wew)

sf_Wew_rot = Flxs_sf[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(2,1,2)]['Vertical'] + \
          Flxs_sf[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(2,1,2)]['Horizontal']
freq_sf, sf_Wew_rot = process_sf(freq_sf, sf_Wew_rot)
cg_Wew_rot = np.interp(freq_cg, freq_sf, sf_Wew_rot)

cg_Eww = Flxs_cg[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(1,2,2)]['Vertical'] + \
          Flxs_cg[exp][(24, 'LF', None, 16, 'HF', None)][(True, False, False)][(1,2,2)]['Horizontal']
freq_cg, cg_Eww = process_cg(freq_cg, cg_Eww)

cg_Eww_rot = Flxs_cg[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(1,2,2)]['Vertical'] + \
          Flxs_cg[exp][(24, 'LF', 'Rot', 16, 'HF', None)][(True, False, False)][(1,2,2)]['Horizontal']
freq_cg, cg_Eww_rot = process_cg(freq_cg, cg_Eww_rot)

cg_scatt = cg_Wew + cg_Wwe
cg_scatt_rot = cg_Wew_rot + cg_Wwe_rot
xlabel = 'k_{h}\hspace{0.04em}L'
ylabel = r'\Pi \cdot k_{h}L'

to_plot=[sf_Www]
size = [2] * 3
style1 = ['-'] * 3
zorder=[2] * 3
color = ['blue']
labels = ['$\Pi^{Www}$']
for i in range(1):
    plot_1d(freq_sf, to_plot[i], ax=fig.axes[0], xscale='log', yscale='linear',  color=color[i], linestyle=style1[i],
            linewidth=size[i], xlabel=xlabel, ylabel=ylabel, label=labels[i], zorder=1, title='')

to_plot=[cg_scatt, cg_scatt_rot, cg_Eww, cg_Eww_rot]
size = [2] * 4
style1 = ['-',(0,(3,3)),'-',(0,(3,3))]
zorder=[2] * 4
color=['darkgreen', 'darkgreen', 'goldenrod', 'goldenrod']
labels = ['$\Pi^{scatt}\ (E\hspace{0.04em}=\hspace{0.04em}LP)$', '$\Pi^{scatt}\ (E\hspace{0.04em}=\hspace{0.04em}LP\\text{-}Rot)$',
          '$\widetilde{\Pi}^{Eww}\ (E\hspace{0.04em}=\hspace{0.04em}Rot)$', '$\widetilde{\Pi}^{Eww}\ (E\hspace{0.04em}=\hspace{0.04em}LP\\text{-}Rot)$']
for i in range(4):
    plot_1d(freq_cg, to_plot[i], ax=fig.axes[0], xscale='log', yscale='linear',  color=color[i], linestyle=style1[i],
            linewidth=size[i], xlabel=xlabel, ylabel=ylabel, label=labels[i], zorder=1, title='')

fig.axes[0].legend(fontsize=sizes.legend_size+1, handlelength=1.5)
# fig.suptitle(Title, fontsize=sizes.title_size-2, y=0.97, x=0.52, horizontalalignment='center')

if no_factor:
    fig.axes[0].set_ylim(-5e-11, 6e-10)
    fig.axes[0].set_xlim(None, 256)
else:
    fig.axes[0].set_ylim(-0.2e-8, 2.4e-8)
fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')

plt.show()
if no_factor:
    save_fig(fig, title='SF_24_16_wave_CG_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_24_16_wave_CG', format='jpg')


