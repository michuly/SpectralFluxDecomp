from scripts.plot_data.plot_SF import *

"""
final plot
"""
##
# determine data plots
depths=None
sigma=24
# depths=range(90,100)

# what to plot
direction = 'Total'
filter_widths1 = [sigma] * 6

filter_widths2 = [0] * 6
freq_domains1 = ['LF'] * 6  # rot
freq_domains2 = [None] * 6  # div
helms1 = [None, 'rot', 'div'] * 2 # rot
helms2 = [None] * 6  # div
coor_shifts = [1] * 6
reduce_means = [False] * 6
exps_A = ['Steady'] * 3 + ['Stochastic'] * 3
ab = 'EE'
temporal = False


# sizes.label_size = 12
# Title=None


def get_label(data_headers, flags_list_plot):
    labels=[]
    for dat_head in data_headers:
        dat_head = dat_head.data_header1
        flags_str = '$%s$' % filter_to_text(
            dat_head.filter_width, dat_head.freq_domain, dat_head.helm, fw=False).replace('-','{\\text -}')
        labels.append(flags_str)
    return labels


data_headers1 = list(map(TwoDataHeader, exps_A, filter_widths1, freq_domains1, helms1, filter_widths2, freq_domains2,
                         helms2, coor_shifts))
flags_list_plot1 = [[(1,1,1)]]*6
labels1=get_label(data_headers1, flags_list_plot1)
color1=['royalblue', 'orangered', 'green'] * 2


###  plotting parameters

if sigma == 24:
    Title = 'Eddy-Eddy interactions'.replace('', '\hspace{0.05em}') + ' $\Pi^{Eee}$'

Title = 'Eddy-Eddy interactions'.replace('', '\hspace{0.05em}') + r' $\Pi^{Eee}\ \ (\hspace{0.05em}\sigma_L = %dhr$\hspace{0.05em})' % sigma
size = [2] * 3 + [1.7] * 3
style1 = [(0, (3, 3))] * 3 + ['-'] * 3
zorder=[2] * 3 + [1] * 3
vmin=-1e-8
vmax=1e-8

sizes.legend_size -=1

###  initializing figure

fig, axes = plt.subplots(nrows=1, ncols=2)
fig.set_figheight(5.2)
fig.set_figwidth(10.5)
fig.subplots_adjust(left=0.08, right=0.98)
fig.subplots_adjust(hspace=0.1, wspace=0.2)

###  Upper panels

figures = {}

figures[0] = {'title':'', 'flags': flags_list_plot1, 'directions': [direction], 'data_headers': data_headers1,
              'labels':labels1, 'color':color1, 'size':size, 'style':style1, 'zorder':zorder}
freq, Flxs, fig = SF_plot(figures, '', temporal, depths=depths, smoothing=1, fig=fig, ax=axes[0])

# fig.suptitle(Title, fontsize=sizes.title_size+2, y=0.97, x=0.48)
# _title = '$(a)$'
_title =''
axes[0].set_title(_title, fontsize=sizes.title_size - 2, fontweight=sizes.title_weight, y=1.02, x=0.5)
axes[0].set_ylim(-0.26e-8, 0.49e-8)
axes[0].set_aspect(1.0 / axes[0].get_data_ratio(), adjustable='box')

lines = fig.axes[0].get_lines()

if no_factor:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                                 fontsize=sizes.legend_size, loc='lower right',
                                 bbox_to_anchor=(0.995, 0),
                                 bbox_transform=fig.axes[0].transAxes, handlelength=1.6)
    bbox = legend1.get_window_extent(renderer=fig.canvas.get_renderer())
    bbox_data_coords = bbox.transformed(fig.axes[0].transAxes.inverted())
    legend_width = bbox_data_coords.width
    legend_height = bbox_data_coords.height
    legend2 = fig.axes[0].legend([lines[i] for i in range(3, 6)], [labels1[i] for i in range(3, 6)], title='COMB',
                                 fontsize=sizes.legend_size, loc='lower right',
                                 bbox_to_anchor=(0.995, legend_height + 0.02),
                                 bbox_transform=fig.axes[0].transAxes, handlelength=1.6)
else:
    legend1 = fig.axes[0].legend([lines[i] for i in range(3)], [labels1[i] for i in range(3)], title='LFW',
                                fontsize=sizes.legend_size, loc=(0.015, 0.758))
    legend2 = fig.axes[0].legend([lines[i] for i in range(3,6)], [labels1[i] for i in range(3,6)], title='COMB',
                                 fontsize=sizes.legend_size, loc=(0.015, 0.51))


fig.axes[0].add_artist(legend1)
fig.axes[0].add_artist(legend2)

### Bottom panels

exp='Stochastic'
sf_comb = Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)]['Vertical'] + \
          Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)]['Horizontal']
# sf_comb = Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)][direction]

exp='Steady'
sf_lfw = Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)]['Vertical'] + \
         Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)]['Horizontal']
# sf_lfw = Flxs[exp][(24, 'LF', None, 0, None, None)][(True, False, False)][(1,1,1)][direction]

sf_lfw = np.concatenate((sf_lfw[:, 0, np.newaxis], savgol_filter(sf_lfw[:, 1:], 15, 3, axis=1)), axis=1)
sf_comb = np.concatenate((sf_comb[:, 0, np.newaxis], savgol_filter(sf_comb[:, 1:], 15, 3, axis=1)), axis=1)
if no_factor:
    pass
else:
    sf_comb = freq * sf_comb
    sf_lfw = freq * sf_lfw

top=300
bot=1000
top_ind=int(128*(2000-top)/2000)
bot_ind=int(128*bot/2000)
print(bot_ind, top_ind)
comb_upper = np.mean(sf_comb[top_ind:-1,:], axis=0)
comb_lower = np.mean(sf_comb[:bot_ind,:], axis=0)
lfw_upper = np.mean(sf_lfw[top_ind:-1,:], axis=0)
lfw_lower = np.mean(sf_lfw[:bot_ind,:], axis=0)
to_plot=[lfw_upper, lfw_lower, comb_upper, comb_lower]

xlabel = 'k_{h}\hspace{0.04em}L'
if no_factor:
    ylabel = r'\Pi'
else:
    ylabel = r'\Pi \cdot k_{h}L'

size = [2] * 2 + [1.7] * 2
style1 = [(0, (3, 3))] * 2 + ['-'] * 2
zorder=[2] * 3 + [1] * 3
color = ['maroon', 'blue', 'maroon', 'blue']
labels = list(map(lambda s: '$%s$' % s, ['LFW\hspace{0.35em}\mathbf{top}', 'LFW\hspace{0.35em}\mathbf{bot}', 'COMB\hspace{0.35em}\mathbf{top}', 'COMB\hspace{0.35em}\mathbf{bot}']))

for i in range(4):
    plot_1d(freq, to_plot[i], ax=fig.axes[1], xscale='log', yscale='linear',  color=color[i], linestyle=style1[i],
            linewidth=size[i], xlabel=xlabel, ylabel=ylabel, label=labels[i], zorder=1, title='')

# _title = '$(b)$'
_title =''
axes[1].set_title(_title, fontsize=sizes.title_size-2, fontweight=sizes.title_weight, y=1.02, x=0.5)
# tx = fig.axes[0].yaxis.get_offset_text()
# tx.set_fontsize(0)
# tx = fig.axes[1].yaxis.get_offset_text()
# tx.set_fontsize(0)

if no_factor:
    axes[0].set_ylim(-8e-10, 3e-10)
    axes[1].set_ylim(-4e-9, 1.5e-9)
    axes[1].legend(loc=4)
    axes[0].text(.98, .97, s=r'$\mathbf{(a)}$', fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight, ha='right', va='top', transform=axes[0].transAxes)

    axes[1].text(.98, .97, s=r'$\mathbf{(b)}$', fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight, ha='right', va='top', transform=axes[1].transAxes)
    # txt1 = fig.axes[0].text(0.5, 3.1e-10, r'x$10^{-10}$', color='black', fontsize=15,
    #                         bbox=dict(facecolor='None', edgecolor='grey', boxstyle='circle,pad=0.5',
    #                                   mutation_aspect=0.6))
    # txt2 = fig.axes[1].text(0.5, 1.6e-9, r'x$10^{-9}$', color='black', fontsize=15,
    #                         bbox=dict(facecolor='None', edgecolor='grey', boxstyle='circle,pad=0.5',
    #                                   mutation_aspect=0.6))
else:
    axes[1].set_ylim(-0.9e-8, 1.15e-8)
    axes[0].text(x=190, y=4.3e-9, s=r'$\mathbf{(%s)}$' % 'a', fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight)
    axes[1].text(x=190, y=9.9e-9, s=r'$\mathbf{(%s)}$' % 'b', fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight)
    # txt = fig.axes[0].text(0.5, 5.1e-9, r'x$10^{-9}$', color='black', fontsize=15,
    #                        bbox=dict(facecolor='None', edgecolor='grey', boxstyle='circle,pad=0.5',
    #                                  mutation_aspect=0.7))
    # txt = fig.axes[1].text(0.5, 1.2e-8, r'x$10^{-8}$', color='black', fontsize=15,
    #                        bbox=dict(facecolor='None', edgecolor='grey', boxstyle='circle,pad=0.5',
    #                                  mutation_aspect=0.7))
    axes[1].legend()
fig.axes[0].set_xlim(0.5,256)
fig.axes[1].set_xlim(0.5,256)
axes[0].set_aspect(1.0 / axes[0].get_data_ratio(), adjustable='box')
axes[1].set_aspect(1.0 / axes[1].get_data_ratio(), adjustable='box')



plt.show()
if no_factor:
    save_fig(fig, title='SF_%d_LF_Eee' % sigma +'_no_factor', format='jpg')
else:
    save_fig(fig, title='SF_%d_LF_Eee' % sigma, format='jpg')