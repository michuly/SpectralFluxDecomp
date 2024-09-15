from scripts.plot_data.plot_SF import *

"""
final plot
"""
### figure
test = 0

sizes.title_size = sizes.title_size-2
# sizes.legend_size_size = 12

### initializing figure

fig, axes = plt.subplots(nrows=1, ncols=3)
fig.set_figheight(4.5)
fig.set_figwidth(12.5)
fig.subplots_adjust(left=0.06, right=0.98, bottom=0.05)
# fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.05)
fig.subplots_adjust(hspace=0.13, wspace=0.18)
# Title = 'Sensitivity to '.replace('', '\hspace{0.05em}') + '$\sigma_{L}$,  $\sigma_{H}$'
Title = ''

## plotting upper fig

# what to plot

direction = 'Total'
filter_widths1 = [48, 24, 16] * 2
filter_widths2 = [16] * 6
freq_domains1 = ['LF'] * 6  # rot
freq_domains2 = ['HF'] * 6  # div
helms1_1 = ['rot'] * 6  # rot
helms1_2 = [None] * 6  # rot
helms2 = [None] * 6  # div
coor_shifts = [True] * 6
reduce_means = [False] * 6
exps_A = ['Steady'] * 3 + ['Stochastic'] * 3
ab = 'EW'
temporal = 0


def get_label(data_headers1, flags_list_plot1):
    labels = []
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        dat_head1 = dat_head.data_header1
        flags_str = '$\sigma_{L}\hspace{0.1em}$=$\hspace{0.1em}%dhr$' % dat_head1.filter_width
        labels.append(flags_str)
    return labels


data_headers1 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1_2, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
data_headers2 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1_1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_lists_plot1 = [[(1, 1, 1)], [(1, 1, 1)], [(1, 1, 1)], [(1, 1, 1)], [(1, 1, 1)], [(1, 1, 1)]]  # # lines
if test:
    flags_lists_plot1 = [[(1, 1, 1)]]
labels1 = get_label(data_headers1, flags_lists_plot1)
color1 = ['orangered', 'darkred', 'goldenrod'] * 3 # olive

# plotting parameters

line_s = [2] * 3 + [1.7] * 3
s = [(0, (3, 3))] * 3 + ['-'] * 3
zorder = [2] * 3 + [1] * 3


### plotting

figures = {}
figures[0] = {'title': '$(a)\ \ \Pi^{Eee}\ (E\ =\ LP)$', 'flags': flags_lists_plot1, 'directions': [direction],
              'data_headers': data_headers1,
              'labels': labels1, 'color': color1, 'style': s, 'size': line_s, 'zorder': zorder}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=None, smoothing=1, fig=fig, ax=fig.axes[0])

figures[0] = {'title': '$(b)\ \ \Pi^{Eee}\ (E\ =\ LP\\text{-}Rot)$', 'flags': flags_lists_plot1, 'directions': [direction],
              'data_headers': data_headers2,
              'labels': labels1, 'color': color1, 'style': s, 'size': line_s, 'zorder': zorder}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=None, smoothing=1, fig=fig, ax=fig.axes[1])

if no_factor:
    fig.axes[0].set_ylim(-7e-10, 4e-10)
    fig.axes[1].set_ylim(-8e-10, 2e-10)
else:
    fig.axes[0].set_ylim(-2.6e-9, 7.3e-9)
    fig.axes[1].set_ylim(-2.5e-9, 2e-9)
fig.axes[0].set_aspect(1 / fig.axes[0].get_data_ratio(), adjustable='box')
fig.axes[1].set_aspect(1 / fig.axes[1].get_data_ratio(), adjustable='box')
fig.axes[1].set_ylabel(None)

lines = fig.axes[0].get_lines()
if no_factor:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW0',
                                 fontsize=sizes.legend_size - 1, loc='lower right', bbox_to_anchor=(0.995, 0),
                                 bbox_transform=fig.axes[0].transAxes, handlelength=1.6)
    bbox = legend1.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox_data_coords = bbox.transformed(fig.axes[0].transAxes.inverted())
    legend_width = bbox_data_coords.width
    legend_height = bbox_data_coords.height
    legend2 = fig.axes[0].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3,6)], title='COMB0',
                                 fontsize=sizes.legend_size - 1, loc='lower right',
                                 bbox_to_anchor=(0.995, legend_height + 0.04),
                                 bbox_transform=fig.axes[0].transAxes, handlelength=1.6)
else:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                             fontsize=sizes.legend_size-1, loc=(0.01, 0.73), handlelength=1.6)
    legend2 = fig.axes[0].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3, 6)], title='COMB',
                             fontsize=sizes.legend_size-1, loc=(0.01, 0.46), handlelength=1.6)
fig.axes[0].add_artist(legend1)
fig.axes[0].add_artist(legend2)


lines = fig.axes[1].get_lines()
if test:
    lines*=6
    labels1*=6
if no_factor:
    legend1 = fig.axes[1].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW1',
                                 fontsize=sizes.legend_size - 1, loc='lower right', bbox_to_anchor=(0.995, 0),
                                 bbox_transform=fig.axes[1].transAxes, handlelength=1.6)
    bbox = legend1.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox_data_coords = bbox.transformed(fig.axes[1].transAxes.inverted())
    legend_width = bbox_data_coords.width
    legend_height = bbox_data_coords.height
    legend2 = fig.axes[1].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3,6)], title='COMB1',
                                 fontsize=sizes.legend_size - 1, loc='lower right',
                                 bbox_to_anchor=(0.995, legend_height + 0.04),
                                 bbox_transform=fig.axes[1].transAxes, handlelength=1.6)
else:
    legend1 = fig.axes[1].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                             fontsize=sizes.legend_size-1, loc=(0.01, 0.73), handlelength=1.6)
    legend2 = fig.axes[1].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3, 6)], title='COMB',
                             fontsize=sizes.legend_size-1, loc=(0.01, 0.46), handlelength=1.6)
fig.axes[1].add_artist(legend1)
fig.axes[1].add_artist(legend2)

## plotting bottom fig

# what to plot
direction = 'Total'

filter_widths1 = [24] * 6
filter_widths2 = [24, 16, 12] * 2
helms1 = [None] * 6 # rot
helms2 = [None] * 6 # div
freq_domains1 = ['LF'] * 6 # rot
freq_domains2 = ['HF'] * 6 # div
coor_shifts = [True] * 6
reduce_means = [False] * 6
exps = ['Stochastic'] * 3
ab = 'EW'
temporal = 0

def get_label(data_headers1, flags_list_plot1):
    labels=[]
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        for plot_flags in flag_list:
            flags_str = ''
            dat_head2 = dat_head.data_header2
            flags_str += '$\sigma_{H}\hspace{0.1em}$=$\hspace{0.1em}%dhr$' % dat_head2.filter_width
            labels.append(flags_str)
    return labels


## plotting parameters

data_headers1 = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_lists_plot1 = ([(2,2,2)],[(2,2,2)],[(2,2,2)])
if test:
    flags_lists_plot1 = ([(2,2,2)],)
labels=get_label(data_headers1, flags_lists_plot1)
color1=['darkred', 'orangered', 'goldenrod'] * 2

size = [2] * 3 + [1.7] * 3
style = ['-'] * 3
zorder = [2] * 3 + [1] * 3

## figure

figures = {}
figures[0] = {'title':'', 'flags': flags_lists_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels, 'color':color1, 'style':style, 'size':size, 'zorder':zorder}

freq, Flxs, fig = SF_plot(figures, Title, temporal, depths=None, smoothing=1, fig=fig, ax=fig.axes[2])


fig.axes[2].set_title('$(c)\ \ \Pi^{Www}\ (W\ =\ HP)$', fontsize=sizes.title_size, fontweight=sizes.title_weight, y=1.02, x=0.5)
if no_factor:
    fig.axes[2].set_ylim(-1e-10, 6e-10)
else:
    fig.axes[2].set_ylim(-0.5e-8, 3.5e-8)
fig.axes[2].set_aspect(1 / fig.axes[2].get_data_ratio(), adjustable='box')
fig.axes[2].set_ylabel(None)


lines = fig.axes[2].get_lines()
if test:
    lines*=6
    labels*=6
if no_factor:
    legend1 = fig.axes[2].legend([lines[i] for i in range(3)], [labels[i] for i in range(3)], title='COMB2',
                                 fontsize=sizes.legend_size - 1, loc='upper right', bbox_to_anchor=(0.998, 0.998),
                                 bbox_transform=fig.axes[2].transAxes, handlelength=1.6)
else:
    legend1 = fig.axes[2].legend([lines[i] for i in range(3)], [labels[i] for i in range(3)], title='COMB',
                             fontsize=sizes.legend_size-1, loc=(0.01, 0.73), handlelength=1.6)
fig.axes[2].add_artist(legend1)

# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.99, x=0.5)

if no_factor:
    save_fig(fig, title='SF_filter_diff_CG_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_filter_diff_CG', format='jpg')
plt.show()

