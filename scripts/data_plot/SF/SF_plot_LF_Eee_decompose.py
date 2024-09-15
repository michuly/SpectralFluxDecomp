from scripts.plot_data.plot_SF import *

##
# determine data plots
depths=None
# depths=range(90,100)

# what to plot
direction = 'Total'
filter_widths1 = [24] * 3
filter_widths2 = [24] * 3
freq_domains1 = ['LF'] * 3  # rot
freq_domains2 = ['LF'] * 3  # div
helms1 = ['rot']  # rot
helms2 = ['div'] * 3  # div
coor_shifts = [1] * 3
reduce_means = [False] * 3
exps_B = ['Stochastic']
exps_A = ['Steady']
ab = 'RD'
temporal = False


sizes.label_size = 12


def get_label(data_headers, flags_list_plot):
    labels=[]
    for dat_head, flags_list in zip(data_headers, flags_list_plot):
        for plot_flags in flags_list:
            if len(np.shape(plot_flags)) == 1:
                plot_flags = (plot_flags,)
            flags_str = '$\Pi^{%s}$' %flags_to_str(plot_flags[0], ab=ab)
            for flags in plot_flags[1:]:
                flags_str = flags_str + '+' + '$\Pi^{%s}$' % flags_to_str(flags, ab=ab)
            # flags_str = flags_str.replace('d','$\delta$')
            # flags_str = flags_str.replace('D','$\Delta$')
            # flags_str = flags_str.replace('r','$\omega$')
            # flags_str = flags_str.replace('R','$\Omega$')
            labels.append(flags_str)
    return labels

### plotting parameters

data_headers1 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
# flags_list_plot1 = ([(1,1,1),((1,2,1),(1,1,2)), (1,2,2), (2,2,2), ((2,2,1),(2,1,2)), (2,1,1)],)
# # flags_list_plot1 = ([(1,1,1)],)
# color1=['orangered', 'darkred', 'orange', 'green', 'royalblue', 'darkturquoise', ]
flags_list_plot1 = ([(1,1,1),((1,2,1),(1,1,2),(2,1,1)), (2,2,2), ((2,2,1),(2,1,2),(1,2,2))],)
color1=['orange', 'darkred', 'green', 'royalblue']
labels1=get_label(data_headers1, flags_list_plot1)

data_headers2 = list(map(TwoDataHeader, exps_B, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
# flags_list_plot2 = ([(1,1,1), ((1,2,1),(1,1,2)), (1,2,2), (2,2,2), ((2,2,1),(2,1,2)), (2,1,1)],)
# # flags_list_plot2 = ([(1,1,1)],)
# color2=['orangered', 'darkred', 'orange', 'green', 'royalblue', 'darkturquoise', ]
flags_list_plot2 = ([(1,1,1),((1,2,1),(1,1,2),(2,1,1)), (2,2,2), ((2,2,1),(2,1,2),(1,2,2))],)
color2=['orange', 'darkred', 'green', 'royalblue']
labels2=get_label(data_headers2, flags_list_plot2)

filter3 = filter_to_text(filter_widths1[0], freq_domains1[0], helms1[0], fw=0)
filter4 = filter_to_text(filter_widths2[0], freq_domains2[0], helms2[0], fw=0)
title3 = '$R = %s\  \ \&\  \ D = %s$' % (filter3.replace('-','{\\text -}'), filter4.replace('-','{\\text -}'))
Title = '$\Pi^{Eee}$' + ' Decomp.  '.replace('', '\hspace{0.053em}') + '(%s)' % title3
line_s = [1.7] *6
s = [(0, (3, 3))] * 3 + ['-'] * 3
s = ['-'] * 4
zorder=[3,3,3,2,2,1]
outline=[False]*6

### initializing figure

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_figheight(5.5)
fig.set_figwidth(10)
# fig.subplots_adjust(left=0.08, right=0.98, top=0.83)
fig.subplots_adjust(left=0.08, right=0.98)
fig.subplots_adjust(hspace=0.1)


### plotting

figures = {}
figures[0] = {'title':'$(a)$\ \ LFW', 'flags': flags_list_plot1, 'directions': [direction], 'data_headers': data_headers1,
               'labels':labels1, 'color':color1, 'style':s, 'size':line_s, 'zorder':zorder, 'outline':outline}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=depths, smoothing=1, fig=fig, ax=axes[0])

figures[0] = {'title':'$(b)$\ \ COMB', 'flags': flags_list_plot2, 'directions': [direction], 'data_headers': data_headers2,
              'labels':labels2, 'color':color2, 'style':s, 'size':line_s, 'zorder':zorder, 'outline':outline}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=depths, smoothing=1, fig=fig, ax=axes[1])

# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.97, x=0.52, horizontalalignment='center')


for i in range(2):
    if no_factor:
        fig.axes[i].set_ylim(-8e-10, 4e-10)
        fig.axes[1].legend(fontsize=sizes.legend_size+3, loc=4)
    else:
        fig.axes[i].set_ylim(-2.6e-9, 4.9e-9)
        fig.axes[0].legend(fontsize=sizes.legend_size+3, loc=2)
    fig.axes[i].set_aspect(1.0 / fig.axes[i].get_data_ratio(), adjustable='box')
plt.show()
if no_factor:
    save_fig(fig, title='SF_24_LF_Eee_decompose_CG_no_factor', format='jpg')
else:
    # save_fig(fig, title='SF_24_LF_Eee_decompose')
    save_fig(fig, title='SF_24_LF_Eee_decompose_CG', format='jpg')
