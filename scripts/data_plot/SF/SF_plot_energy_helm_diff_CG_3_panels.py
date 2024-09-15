from scripts.plot_data.plot_SF import *

"""
final plot
"""

Title = 'Role of Rot. and Div. Components'.replace('', '\hspace{0.05em}')
## plotting upper fig

# what to plot

direction = 'Total'
filter_widths1 = [24] * 4
filter_widths2 = [16] * 4
helms1 = ['rot', 'rot', 'div', 'div']  # rot
helms2 = ['rot', 'div', 'rot', 'div'] # div
freq_domains1 = ['LF'] * 4  # rot
freq_domains2 = ['HF'] * 4  # div
coor_shifts = [True] * 4
reduce_means = [False] * 4
exps1 = ['Stochastic'] * 4
ab = 'EW'
temporal = 0


def get_label(data_headers1, flags_list_plot1):
    labels = []
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        for plot_flags in flag_list:
            flags_str = ''
            # flags_str = '%s - ' % exp_name[dat_head.exp]
            if len(np.shape(plot_flags)) == 1:
                plot_flags = (plot_flags,)
            dat_head1 = dat_head.data_header1
            dat_head2 = dat_head.data_header2
            filter1 = filter_to_text(dat_head1.filter_width, dat_head1.freq_domain, dat_head1.helm, fw=False).replace('-','{\\text -}')
            filter2 = filter_to_text(dat_head2.filter_width, dat_head2.freq_domain, dat_head2.helm, fw=False).replace('-','{\\text -}')
            flags_str += '$%s$\ \&\ $%s$' % (filter1, filter2)
            labels.append(flags_str)
    return labels

## plotting parameters

data_headers1 = list(map(TwoDataHeader, exps1, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))

flags_lists_plot1 = ([(2, 1, 2)], [(2, 1, 2)], [(2, 1, 2)], [(2, 1, 2)])  # figure # filter # line # flags
flags_lists_plot1_cg = ([(2, 2, 1), (1,2,2)], [(2, 2, 1), (1,2,2)], [(2, 2, 1), (1,2,2)], [(2, 2, 1), (1,2,2)] )  # figure # filter # line # flags
# flags_lists_plot1 = ([((2, 1, 2), (2, 2, 1))],)  # figure # filter # line # flags
labels1 = get_label(data_headers1, flags_lists_plot1)
labels1_cg = get_label(data_headers1, flags_lists_plot1_cg)
color1 = ['darkred', 'goldenrod', 'darkred', 'orangered', 'goldenrod', 'olive', 'darkgreen', 'blue'] * 2
line_s = [2] * 2 + [1.7] * 6
s1 = [(0,(3,3))] * 2 + ['-'] * 6
zorder = [2] * 8

### figure


figures = {}
figures[0] = {'title': '$(a)\ \ \Pi^{scatt} SF$', 'flags': flags_lists_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels': labels1, 'color': color1, 'style': s1, 'size': line_s, 'zorder': zorder}
freq_sf, Flxs_sf, fig1 = SF_plot(figures, '', temporal, depths=None, smoothing=1, fig=None, ax=None)
fig1.axes[0].legend()
figures[0] = {'title': '$(a)\ \ \widetilde{\Pi}^{scatt} CG$', 'flags': flags_lists_plot1_cg, 'directions': [direction], 'data_headers': data_headers1,
              'labels': labels1_cg, 'color': color1, 'style': s1, 'size': line_s, 'zorder': zorder}
freq_cg, Flxs_cg, fig2 = SF_plot(figures, '', temporal, depths=None, smoothing=1, flux_remainder='CG_ray', fig=None, ax=None)
fig2.axes[0].legend()


######

### initializing figure

fig, axes = plt.subplots(nrows=1, ncols=3)
fig.set_figheight(4.5)
fig.set_figwidth(12.5)
# fig.subplots_adjust(left=0.07, right=0.98, top=0.92, bottom=0.03)
fig.subplots_adjust(left=0.07, right=0.98, bottom=0.03)
fig.subplots_adjust(hspace=0.13, wspace=0.2)

def process_cg(freq_tmp, flx_tmp):
    index = np.array([2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 80, 112, 144, 176, 240, 368, 496]) - 1
    flx_tmp = np.nanmean(flx_tmp, axis=0)
    flx_tmp = np.concatenate((np.zeros(1), np.flip(flx_tmp[index])))
    # flx_tmp = np.concatenate((flx_tmp[0:6], savgol_filter(flx_tmp, 7, 3)[6:17], flx_tmp[17:19]))
    flx_tmp = np.concatenate((savgol_filter(flx_tmp, 5, 3)[0:6], savgol_filter(flx_tmp, 5, 2)[6:16],
                              savgol_filter(flx_tmp, 5, 3)[16:19]))
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


######### plotting left fig

colors=['darkred', 'tab:red', 'tab:orange', 'gold']
exp='Stochastic'
for i, header in enumerate([(24, 'LF', 'rot', 16, 'HF', 'rot'), (24, 'LF', 'rot', 16, 'HF', 'div'),
                            (24, 'LF', 'div', 16, 'HF', 'rot'), (24, 'LF', 'div', 16, 'HF', 'div')]):

    sf_Wew = Flxs_sf[exp][header][(True, False, False)][(2, 1, 2)]['Vertical'] + \
             Flxs_sf[exp][header][(True, False, False)][(2, 1, 2)]['Horizontal']
    freq_sf, sf_Wew = process_sf(freq_sf, sf_Wew)
    cg_Wew = np.interp(freq_cg, freq_sf, sf_Wew)

    cg_Wwe = Flxs_cg[exp][header][(True, False, False)][(2, 2, 1)]['Vertical'] + \
              Flxs_cg[exp][header][(True, False, False)][(2, 2, 1)]['Horizontal']
    freq_cg, cg_Wwe = process_cg(freq_cg, cg_Wwe)

    cg_Eww = Flxs_cg[exp][header][(True, False, False)][(1, 2, 2)]['Vertical'] + \
              Flxs_cg[exp][header][(True, False, False)][(1, 2, 2)]['Horizontal']
    freq_cg, cg_Eww = process_cg(freq_cg, cg_Eww)

    cg_scatt = cg_Wew + cg_Wwe
    xlabel = 'k_{h}\hspace{0.04em}L'
    ylabel = r'\Pi \cdot k_{h}L'
    filter1 = filter_to_text(header[0], header[1], header[2], fw=False).replace('-','\\text{-}')
    filter2 = filter_to_text(header[3], header[4], header[5], fw=False).replace('-','\\text{-}')
    label = '$%s\ \&\ %s$' % (filter1, filter2)
    plot_1d(freq_cg, cg_scatt, ax=fig.axes[0], xscale='log', yscale='linear',  color=colors[i], linestyle='-',
            linewidth=2, xlabel=xlabel, ylabel=ylabel, label=label, zorder=1, title='')

fig.axes[0].legend(fontsize=sizes.legend_size-2, handlelength=1.2)
fig.axes[0].set_title('$(a)\ \ \ \Pi^{scatt}$', fontsize=sizes.title_size-1, fontweight=sizes.title_weight, y=1.02, x=0.5)
fig.axes[0].set_ylim(-0.3e-8, 1.3e-8)
fig.axes[0].set_xlim(0.5, None)
fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')


######### plotting middle fig

colors=['darkred', 'tab:red', 'tab:orange', 'gold']
exp='Stochastic'
for i, header in enumerate([(24, 'LF', 'rot', 16, 'HF', 'rot'), (24, 'LF', 'rot', 16, 'HF', 'div'),
                            (24, 'LF', 'div', 16, 'HF', 'rot'), (24, 'LF', 'div', 16, 'HF', 'div')]):

    cg_Eww = Flxs_cg[exp][header][(True, False, False)][(1, 2, 2)]['Vertical'] + \
              Flxs_cg[exp][header][(True, False, False)][(1, 2, 2)]['Horizontal']
    freq_cg, cg_Eww = process_cg(freq_cg, cg_Eww)

    xlabel = 'k_{h}\hspace{0.04em}L'
    ylabel = r'\Pi \cdot k_{h}L'
    filter1 = filter_to_text(header[0], header[1], header[2], fw=False).replace('-','\\text{-}')
    filter2 = filter_to_text(header[3], header[4], header[5], fw=False).replace('-','\\text{-}')
    label = '$%s\ \&\ %s$' % (filter1, filter2)
    plot_1d(freq_cg, cg_Eww, ax=fig.axes[1], xscale='log', yscale='linear',  color=colors[i], linestyle='-',
            linewidth=2, xlabel=xlabel, ylabel=ylabel, label=label, zorder=1, title='')

fig.axes[1].legend(fontsize=sizes.legend_size-2, handlelength=1.2)
fig.axes[1].set_title('$(b)\ \ \ \widetilde{\Pi}^{Eww}$', fontsize=sizes.title_size-1, fontweight=sizes.title_weight, y=1.02, x=0.5)
fig.axes[1].set_ylim(-1e-8, 1.3e-8)
fig.axes[1].set_xlim(0.5, None)
fig.axes[1].set_ylabel(None)
fig.axes[1].set_aspect(1.0 / fig.axes[1].get_data_ratio(), adjustable='box')


## plotting right fig

# what to plot
direction = 'Total'
filter_widths1 = [24] * 6
filter_widths2 = [16] * 6
freq_domains1 = ['LF'] * 6
freq_domains2 = ['HF'] * 6
helms1 = [None] * 6
helms2 = [None, 'rot', 'div'] * 2
coor_shifts = [True] * 6
reduce_means = [False] * 6
exps1 = ['Stochastic'] * 3
ab = 'EW'
temporal = 0

def get_label(data_headers1, flags_list_plot1):
    labels=[]
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        dat_head2 = dat_head.data_header2
        filter2 = filter_to_text(dat_head2.filter_width, dat_head2.freq_domain, dat_head2.helm, fw=False).replace('-','{\\text -}')
        flags_str = '$%s$' % filter2
        labels.append(flags_str)
    return labels

## plotting parameters

data_headers1 = list(map(TwoDataHeader, exps1, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))

flags_lists_plot1 = ([(2,2,2)],[(2,2,2)],[(2,2,2)])
# flags_lists_plot1 = ([(2,2,2)], )
labels1=get_label(data_headers1, flags_lists_plot1)
color1=['darkred', 'orangered', 'goldenrod']
color1=['darkred', 'tab:orange', 'gold']

size = [2] * 3
style1 = ['-'] * 3
zorder = [2] * 3

## figure

figures = {}
figures[0] = {'title': '', 'flags': flags_lists_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels': labels1, 'color': color1, 'style': style1, 'size': size, 'zorder': zorder}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=None, smoothing=1, fig=fig, ax=axes[2])
axes[2].set_ylim(-0.5e-8, 2.5e-8)
axes[2].set_aspect(1.0 / axes[2].get_data_ratio(), adjustable='box')
fig.axes[2].set_ylabel(None)
fig.axes[2].legend(fontsize=sizes.legend_size-2, handlelength=1.2)
fig.axes[2].set_title('$(c)\ \ \ \Pi^{Www}$', fontsize=sizes.title_size-1, fontweight=sizes.title_weight, y=1.02, x=0.5)

if no_factor:
    fig.axes[0].set_ylim(-1e-10, 4e-10)
    fig.axes[0].set_aspect(1 / fig.axes[0].get_data_ratio(), adjustable='box')
    fig.axes[1].set_ylim(-6.5e-10, 12e-10)
    fig.axes[1].set_aspect(1 / fig.axes[1].get_data_ratio(), adjustable='box')
    fig.axes[2].set_ylim(-1e-10, 4e-10)
    fig.axes[2].set_aspect(1 / fig.axes[2].get_data_ratio(), adjustable='box')


#############

# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.97, x=0.52, horizontalalignment='center')
if no_factor:
    save_fig(fig, title='SF_helm_diff_CG_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_helm_diff_CG', format='jpg')

plt.show()
