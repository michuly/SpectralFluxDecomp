from scripts.plot_data.plot_SF import *


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
sizes.title_size = 18
sizes.label_size = 16
sizes.tick_size = 14
sizes.legend_size = 12
sizes.title_weight = 600

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

fig, ax = plt.subplots(nrows=1, ncols=1)
fig.subplots_adjust(left=0.1, right=0.98)
fig.set_figheight(6)
fig.set_figwidth(6.2)

####  ax 1  ####
##
# determine data plots
figures = {}
figures[0] = {'title':'', 'flags': flags_list_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels1, 'color':color1, 'style':s, 'size':line_s, 'zorder':zorder}

kh, Flxs, fig = SF_plot(figures, Title, temporal=0, smoothing=1, fig=fig, ax=ax)
fig.axes[0].legend()

# fig.axes[1].grid(visible=True, which='major') # color='b', linestyle='-'
if no_factor:
    fig.axes[0].set_ylim(-1.1e-9, 2.2e-9)
    fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')

plt.show()
if no_factor:
    save_fig(fig, title='SF_All_spatial_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_full_spatial', format='jpg')

