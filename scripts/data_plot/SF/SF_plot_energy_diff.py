from scripts.plot_data.plot_SF import *

"""
final plot
"""

##
# determine data plots

# what to plot
direction = 'Total'
filter_widths1 = [0, 24, 24] * 2
filter_widths2 = [0, 16, 16] * 2
freq_domains1 = [None, 'LF', 'LF'] * 2 # rot
freq_domains2 = [None, 'HF', 'HF'] * 2 # div
helms1 = [None, None, 'rot'] * 2 # rot
helms2 = [None, None, None] * 2 # div
coor_shifts = [True] * 6
reduce_means = [False] * 6
exps = ['Steady'] * 3 + ['Stochastic'] * 3
ab = 'EW'
temporal = 0

# plotting parameters
sizes.title_size = 18
sizes.label_size = 16
sizes.tick_size = 14

def get_label(data_headers1, flags_list_plot1):
    labels=[]
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        for plot_flags in flag_list:
            flags_str = ''
            # flags_str = '%s - ' % exp_name[dat_head.exp]
            dat_head1 = dat_head.data_header1
            dat_head2 = dat_head.data_header2
            if len(np.shape(plot_flags)) == 1:
                plot_flags = (plot_flags,)
            if plot_flags == ((3, 3, 3),):
                filter1 = filter_to_text(dat_head1.filter_width, dat_head1.freq_domain, dat_head1.helm, fw=True)
                filter2 = filter_to_text(dat_head2.filter_width, dat_head2.freq_domain, dat_head2.helm, fw=True)
                flags_str += '$%s\hspace{0.25em}\&\hspace{0.25em}%s$' % (filter1.replace('$',''), filter2.replace('$',''))
            elif plot_flags == ((0, 0, 0),):
                flags_str += '$\Pi$'

            labels.append(flags_str)
    return labels


### plotting parameters


data_headers1 = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_lists_plot1 = ([(0,0,0)],[(3,3,3)],[(3,3,3)],[(0,0,0)],[(3,3,3)],[(3,3,3)])
labels1=get_label(data_headers1, flags_lists_plot1)
color1=['k', 'darkred', 'orangered'] * 2

Title = 'Differences in'.replace('', '\hspace{0.05em}') + ' $\Pi$'
line_s = [2] * 3 + [1.7] * 3
s = [(0,(3,3))] * 3 + ['-'] * 3
zorder=[2] * 3 + [1] * 3

### initializing figure

fig, axis = plt.subplots(nrows=1, ncols=1)
fig.set_figheight(6)
fig.set_figwidth(6)
# fig.subplots_adjust(left=0.08, right=0.98, top=0.83)


figures = {}
figures[0] = {'title':'', 'flags': flags_lists_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels1, 'color':color1, 'style':s, 'size':line_s, 'zorder':zorder}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=None, smoothing=True, fig=fig, ax=axis, testing=False)

# fig.suptitle(Title, fontsize=sizes.title_size, y=0.97, x=0.52, horizontalalignment='center')
for i in range(1):
    if no_factor:
        fig.axes[i].set_ylim(-1e-9, 2.5e-9)
    else:
        fig.axes[i].set_ylim(-0.5e-8, 6.5e-8)
    fig.axes[i].set_aspect(1.0 / fig.axes[i].get_data_ratio(), adjustable='box')
lines = fig.axes[0].get_lines()
if no_factor:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                                 fontsize=sizes.legend_size - 1, handlelength=1.5, loc='upper right', bbox_to_anchor=(0.995, 0.995),
                                 bbox_transform=fig.axes[0].transAxes)
    legend2 = fig.axes[0].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3, 6)], title='COMB',
                                 fontsize=sizes.legend_size - 1, handlelength=1.5, loc='lower right', bbox_to_anchor=(0.995, 0),
                                 bbox_transform=fig.axes[0].transAxes)
else:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                                 fontsize=sizes.legend_size, loc=(0.015, 0.78))
    legend2 = fig.axes[0].legend([lines[i] for i in range(3,6)], [labels1[i] for i in range(3,6)], title='COMB',
                                 fontsize=sizes.legend_size, loc=(0.015, 0.554))
fig.axes[0].add_artist(legend1)
fig.axes[0].add_artist(legend2)

plt.show()
if no_factor:
    save_fig(fig, title='SF_energy_sum_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_energy_sum', format='jpg')
