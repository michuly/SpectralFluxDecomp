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
helms1 = [None]  # rot
helms2 = [None]  # div
coor_shifts = [True]
reduce_means = [False]
exps_A = ['Stochastic']
exps_A = ['Steady']
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
            elif plot_flags == (2,1,2):
                flags_str = '$\Pi^{Wew}$'
            elif plot_flags == (1,2,1):
                flags_str = '$\Pi^{Ewe}$'
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


data_headers3 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
# flags_lists_plot3 = ([(1,1,1)],)
flags_lists_plot3 = ([(1,1,1), (2,2,2), (2,1,2), (1,2,1), (2,1,1), (1,2,2), (1,1,2), (2,2,1)],)

labels3=get_label(data_headers3, flags_lists_plot3)
color3=['red', 'royalblue', 'green', 'orange', 'red', 'royalblue', 'green', 'orange']
size3 = [1.7] *8
style3 = ['-']*8


### plotting titles

Title = '$\Pi$' + ' Decompositions'.replace('', '\hspace{0.053em}')
filter3 = filter_to_text(filter_widths1[0], freq_domains1[0], helms1[0], fw=0)
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
figures[0] = {'title': '$(b)$\ \ COMB', 'flags': flags_lists_plot3, 'directions': [direction], 'data_headers': data_headers3,
              'labels': labels3, 'color': color3, 'style':style3, 'size':size3, 'zorder': [1]*8}
freq, Flxs3, fig = SF_plot(figures, '', temporal, depths, smoothing=1, fig=fig, ax=axes[1])

plt.show()
save_fig(fig, title='SF_16_24_All_Test', format='jpg')

