from scripts.plot_data.plot_SF import *

"""
final plot
"""

# data parameters
# what to plot
direction = 'Total'
filter_widths1 = [0] * 4
filter_widths2 = [0] * 4
freq_domains1 = ['LF'] * 4  # rot
freq_domains2 = ['HF'] * 4  # div
helms1 = [None] * 4  # rot
helms2 = [None] * 4  # div
coor_shifts = [True] * 2
reduce_means = [False] * 4
exps = ['Steady', 'Stochastic']
ab = 'EW'

# plotting parameters
# sizes.title_size = 18
# sizes.label_size = 16
# sizes.tick_size = 14
# sizes.legend_size = 12
# sizes.title_weight = 600

def get_label(data_headers1, flags_list_plot1):
    labels=[]
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        for plot_flags in flag_list:
            flags_str = exp_name[dat_head.exp]
            labels.append(flags_str)
    return labels

data_headers1 = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_list_plot1 = [[(0,0,0)], [(0,0,0)]]
labels1=get_label(data_headers1, flags_list_plot1)
color1 = ['royalblue', 'darkred']

# plotting parameters
Title = ''
line_s = [1.7] * 2
s = ['-'] * 2
zorder=[2]*4

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_figheight(5.5)
fig.set_figwidth(10)
fig.subplots_adjust(left=0.08, right=0.98, hspace=0.1)

####  ax 1  ####
##
# determine data plots
figures = {}
figures[0] = {'title':'$(a)$\ \ $\Pi (k_h)$', 'flags': flags_list_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels1, 'color':color1, 'style':s, 'size':line_s, 'zorder':zorder}

kh, Flxs, fig = SF_plot(figures, Title, temporal=0, smoothing=1, fig=fig, ax=axes[0])

####  ax 2  ####

direction = 'Total'
filter_widths1 = [0] * 4
filter_widths2 = [0] * 4
freq_domains1 = ['LF'] * 4  # rot
freq_domains2 = ['HF'] * 4  # div
helms1 = [None] * 4  # rot
helms2 = [None] * 4  # div
coor_shifts = [True] * 2
reduce_means = [False] * 4
exps = ['Steady', 'Stochastic'] * 2
ab = 'EW'

# plotting parameters
# sizes.title_size = 18
# sizes.label_size = 16
# sizes.tick_size = 14
# sizes.title_weight = 600

def get_label(data_headers1, flags_list_plot1):
    labels=[]
    for dat_head, flag_list in zip(data_headers1, flags_list_plot1):
        for plot_flags in flag_list:
            flags_str = exp_name[dat_head.exp]
            if not dat_head.coor_shift:
                flags_str += '\hspace{0.5em}[before coor shift]'
            labels.append(flags_str)
    return labels

data_headers1 = list(map(TwoDataHeader, exps, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_list_plot1 = [[(0,0,0)],[(0,0,0)],[(0,0,0)],[(0,0,0)]]
labels1=get_label(data_headers1, flags_list_plot1)
color1 = ['royalblue', 'royalblue', 'darkred', 'darkred', ]
line_s = [1.7, 2] * 4
s = ['-', (0, (3, 2))] * 2
color1 = ['royalblue', 'darkred']
line_s = [2] * 4
s = ['-'] * 2

figures = {}
figures[0] = {'title':r'$(b)$\ \ ${\Pi} (\omega)$', 'flags': flags_list_plot1, 'directions': [direction],
              'data_headers': data_headers1, 'labels':labels1, 'color':color1, 'style':s, 'size':line_s, 'zorder':zorder}

freq, Flxs, fig = SF_plot(figures, Title, temporal=1, smoothing=1, fig=fig, ax=axes[1])

fig.axes[1].legend(loc=0, fontsize=sizes.legend_size, handlelength=1.6)

fig.axes[1].axvline(x=1, color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2,
           label='$f\hspace{0.05em},  (16hr)^{-1}\hspace{0.05em},\\newline(24hr)^{-1}\hspace{0.05em},  (48hr)^{-1}$')
fig.axes[1].axvline(x=1 / 16 / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
fig.axes[1].axvline(x=1 / 24 / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
fig.axes[1].axvline(x=1 / 48 / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)
# fig.axes[1].axvline(x=1 / (24*7) / 3600 / (F_cor / 2 / np.pi), color='0.25', linestyle=(0, (4, 3)), linewidth=1.7, zorder=2)

# fig.axes[1].grid(visible=True, which='major') # color='b', linestyle='-'
# fig.axes[1].grid(visible=True, which='minor', axis='y')

if no_factor:
    fig.axes[0].set_ylim(-1.1e-9, 2.2e-9)
    fig.axes[1].set_ylim(-1.1e-9, 2.2e-9)
    fig.axes[1].set_xlim(1 / (24*6) / 3600 / (F_cor / 2 / np.pi), None)
    fig.axes[1].set_aspect(1.0 / fig.axes[1].get_data_ratio(), adjustable='box')
    fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')
    
fig.axes[1].set_ylabel('')

plt.show()
if no_factor:
    save_fig(fig, title='SF_All_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_full', format='jpg')

