from scripts.plot_data.plot_SF import *

"""
final plot
"""

##
# determine data plots
depths=None

# what to plot
direction = 'Total'
filter_widths1 = [24]
filter_widths2 = [16]
freq_domains1 = ['LF']  # rot
freq_domains2 = ['HF']  # div
helms1_1 = ['rot']  # rot
helms1_2 = [None]  # rot
helms2 = [None]  # div
coor_shifts = [True]
reduce_means = [False]
exps_A = ['Stochastic']
exps_B = ['Steady']
ab = 'EW'
temporal = False

# plotting parameters
Title=None
# from matplotlib.mathtext import _mathtext as mathtext
# print(mathtext.FontConstantsBase.sup1)
# print(mathtext.SHRINK_FACTOR)
# mathtext.FontConstantsBase.sup1 = 0.5
# mathtext.SHRINK_FACTOR = 0.3
# print(mathtext.FontConstantsBase.sup1)
# print(mathtext.SHRINK_FACTOR)

def get_label(data_headers, flags_list_plot):
    labels=[]
    for dat_head, flag_list in zip(data_headers, flags_list_plot):
        for plot_flags in flag_list:
            if plot_flags == (1, 1, 1):
                flags_str = '$\Pi^{Eee}$'
            elif plot_flags == (2,2,2):
                flags_str = '$\Pi^{Www}$'
            if plot_flags == (1, 2, 1):
                flags_str = '$\Pi^{Ewe}$'
            elif plot_flags == (2,1,2):
                flags_str = '$\Pi^{Wew}$'
            elif plot_flags == ((2,2,1),(2,1,2)):
                flags_str = '$\Pi^{scatt}$'
            elif plot_flags == ((1,2,1),(1,1,2)):
                flags_str = '$\Pi^{Ewe}+\Pi^{Eew}$'
            elif plot_flags == (1,2,2):
                flags_str = '$\Pi^{Eww}$'
            elif plot_flags == (2,1,1):
                flags_str = '$\Pi^{Wee}$'
            elif plot_flags == (2,2,1):
                flags_str = '$\Pi^{Wwe}$'
            elif plot_flags == (1,1,2):
                flags_str = '$\Pi^{Eew}$'
            elif plot_flags == ((1,2,1),(1,1,2),(2,1,1)):
                flags_str = '$\Pi^{Ewe}+\Pi^{Eew}+\Pi^{Wee}$'
            elif plot_flags == ((1,2,2),(1,2,1),(1,1,2),(2,1,1)):
                flags_str = 'res.'
            elif plot_flags == (3, 3, 3):
                flags_str = 'sum'
            elif plot_flags == (0,0,0):
                flags_str = 'total'
            labels.append(flags_str)
    return labels

### plotting parameters
SF_CG_TEST=True

data_headers3 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1_2, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_lists_plot3 = ([(2,2,1), (1,2,2), (1,1,2), (2,1,1)],)
flags_lists_plot3_cg = ([(2,2,1), (1,2,2), (1,1,2), (2,1,1)],)
if SF_CG_TEST:
    flags_lists_plot3 = ([(2,2,2), (1,2,1), (2,1,2), (1,1,1)],)
    flags_lists_plot3_cg = ([(2,2,2), (1,2,1), (2,1,2), (1,1,1)],)

labels3=get_label(data_headers3, flags_lists_plot3)
labels3_cg=list(map(lambda s: s.replace('\Pi','\widetilde{\Pi}'), get_label(data_headers3, flags_lists_plot3)))
color3=['royalblue', 'green', 'maroon', 'orange']
color3_cg=['royalblue', 'green', 'maroon', 'orange']
size3 = [1.7] *4
style3 = ['-']*4
size3_cg = [2] *4
style3_cg = [(0, (3,3))]*4

data_headers4 = list(map(TwoDataHeader, exps_B, filter_widths1, freq_domains1, helms1_2, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_lists_plot4 = ([(2,2,1), (1,2,2), (1,1,2), (2,1,1)],)
flags_lists_plot4_cg = ([(2,2,1), (1,2,2), (1,1,2), (2,1,1)],)
if SF_CG_TEST:
    flags_lists_plot4 = ([(2,2,2), (1,2,1), (2,1,2), (1,1,1)],)
    flags_lists_plot4_cg = ([(2,2,2), (1,2,1), (2,1,2), (1,1,1)],)

labels4=get_label(data_headers4, flags_lists_plot4)
labels4_cg=list(map(lambda s: s.replace('\Pi','\\widetilde{\Pi}'), get_label(data_headers4, flags_lists_plot4)))
color4=['royalblue', 'green', 'maroon', 'orange']
color4_cg=['royalblue', 'green', 'maroon', 'orange']
size4 = [1.7] *4
style4 = ['-']*4
size4_cg = [2] *4
style4_cg = [(0, (3,3))]*4
# style4_cg = [(0,(4,4))]*4


### plotting titles

Title = '$\Pi$' + ' Decompositions'.replace('', '\hspace{0.053em}')
filter3 = filter_to_text(filter_widths1[0], freq_domains1[0], helms1_2[0], fw=0)
filter4 = filter_to_text(filter_widths2[0], freq_domains2[0], helms2[0], fw=0)
title3 = '$E = %s\  \ \&\  \ W = %s$' % (filter3, filter4)
Title = '$\Pi$' + ' Decomp.  '.replace('', '\hspace{0.053em}') + '(%s)' % title3

### initializing figure

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_figheight(5.5)
fig.set_figwidth(10)
fig.subplots_adjust(left=0.08, right=0.98, hspace=0.1)
# fig.subplots_adjust(left=0.08, right=0.98, top=0.83, hspace=0.1)


figures = {}
figures[0] = {'title': '$(a)$\ \ LFW', 'flags': flags_lists_plot4, 'directions': [direction], 'data_headers': data_headers4,
              'labels': labels4, 'color': color4, 'style':style4, 'size':size4, 'zorder': [1]*4}
freq, Flxs4, fig = SF_plot(figures, '', temporal, depths, smoothing=0, fig=fig, ax=axes[0])

figures[0] = {'title': '$(a)$\ \ LFW', 'flags': flags_lists_plot4_cg, 'directions': [direction], 'data_headers': data_headers4,
              'labels': labels4_cg, 'color': color4_cg, 'style':style4_cg, 'size':size4_cg, 'zorder': [1]*4}
freq_cg, Flxs4_cg, fig = SF_plot(figures, '', temporal, depths, smoothing=0, fig=fig, ax=axes[0], flux_remainder='CG_ray')

figures[0] = {'title': '$(b)$\ \ COMB', 'flags': flags_lists_plot3, 'directions': [direction], 'data_headers': data_headers3,
              'labels': labels3, 'color': color3, 'style':style3, 'size':size3, 'zorder': [1]*4}
freq, Flxs3, fig = SF_plot(figures, '', temporal, depths, smoothing=0, fig=fig, ax=axes[1])

figures[0] = {'title': '$(b)$\ \ COMB', 'flags': flags_lists_plot3_cg, 'directions': [direction], 'data_headers': data_headers3,
              'labels': labels3_cg, 'color': color3_cg, 'style':style3_cg, 'size':size3_cg, 'zorder': [1]*4}
freq_cg, Flxs3_cg, fig = SF_plot(figures, '', temporal, depths, smoothing=0, fig=fig, ax=axes[1], flux_remainder='CG_ray')

fig.axes[0].legend(fontsize=sizes.legend_size+1, loc=2, handlelength=1.5)
# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.97, x=0.52, horizontalalignment='center')

lines = fig.axes[0].get_lines()
if no_factor:
    legend1 = fig.axes[0].legend([lines[i] for i in range(4)], [labels3[i] for i in range(4)],
                                 fontsize=sizes.legend_size+1.1, handlelength=1.5,
                                 loc='upper right', bbox_to_anchor=(0.997, 0.997),
                                 bbox_transform=fig.axes[0].transAxes)
    bbox = legend1.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox_data_coords = bbox.transformed(fig.axes[0].transAxes.inverted())
    legend_width = bbox_data_coords.width
    legend_height = bbox_data_coords.height
    legend2 = fig.axes[0].legend([lines[i] for i in range(4,8)], [labels3_cg[i] for i in range(4)],
                                 fontsize=sizes.legend_size+1, handlelength=1.5,
                                 loc='upper right', bbox_to_anchor=(0.997 - (legend_width + 0.02), 0.997),
                                 bbox_transform=fig.axes[0].transAxes)
else:
    legend1 = fig.axes[0].legend([lines[i] for i in range(4)], [labels3[i] for i in range(4)],
                                 fontsize=sizes.legend_size+1.1, loc=(0.72, 0.71), handlelength=1.5)
    legend2 = fig.axes[0].legend([lines[i] for i in range(4,8)], [labels3_cg[i] for i in range(4)],
                                 fontsize=sizes.legend_size+1, loc=(0.46, 0.71), handlelength=1.5)
fig.axes[0].add_artist(legend1)
fig.axes[0].add_artist(legend2)

if no_factor:
    fig.axes[1].set_ylim(-6e-10, 9e-10)
    fig.axes[0].set_ylim(-4e-11, 4e-11 * 9 / 6)
else:
    fig.axes[1].set_ylim(-0.2e-8, 1.2e-8)
    fig.axes[0].set_ylim(-0.3e-9, 1.5e-9)


fig.axes[1].set_xlim(freq[1], np.max(freq))
fig.axes[0].set_xlim(freq[1], np.max(freq))
fig.axes[1].set_aspect(1.0 / fig.axes[1].get_data_ratio(), adjustable='box')
fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')


plt.show()

if no_factor:
    save_fig(fig, title='SF_16_24_All_CG_err_no_factor', format='jpg')
else:
    # save_fig(fig1, title='SF_16_24_rot_All')
    save_fig(fig, title='SF_16_24_All_CG_err', format='jpg')
