from scripts.plot_data.plot_PSD import *

"""
final plot
"""
###
fig, axes = plt.subplots(nrows=1, ncols=3)
fig.set_figheight(5)
fig.set_figwidth(16)
fig.subplots_adjust(left=0.08, right=0.98)
fig.subplots_adjust(hspace=0.13, wspace=0.1)


# plotting parameters
sizes.title_size = 16
sizes.label_size = 16
sizes.legend_size = 13
sizes.tick_size = 14
sizes.title_weight = 600


### ax 1

# data parameters
experiments = ['Steady'] * 2 + ['Stochastic'] * 2
freq_domains = [None] *4
filter_widths = [0] * 4
helm_domains = [None] *4
geostrophics = [False] * 4
psd_types = ['KE'] * 4
dims = [1] * 4
reduce_means = [False] * 4
coor_shifts = [True, False] * 2
hours_cut = 24
data_headers = map(DataHeader, experiments, filter_widths, freq_domains, helm_domains, coor_shifts, geostrophics, dims,
                   reduce_means)

depth_orientation = False
c = ['skyblue', 'royalblue', 'lightcoral', 'darkred', ]
c = ['royalblue', 'royalblue', 'darkred', 'darkred', ]
line_s = [1.7, 2] * 4
s = ['-', (0, (3, 2))] * 2
Title = ''
labels=[]
for exp, coor_shift in zip(experiments, coor_shifts):
    _str = exp_name[exp]
    if not coor_shift:
        _str += '\hspace{0.5em}[before coor shift]'
    labels.append(_str)

print(fig.axes[0])
print(fig.axes[1])

figures = {}
# _title='$(a)$   Freq. Spectra'
_title=''
figures[0] = {'title': _title, 'filters': data_headers, 'labels':labels}
psd_dict, freq_dict, fig = PSD_plot(figures, Title, hours_cut, depth_orientation, c, line_s, s, smoothing=1,
                                    fig=fig, ax=fig.axes[0])

fig.axes[0].legend(loc=3, fontsize=sizes.legend_size, handlelength=1.6)

# fig.axes[0].set_xlim(0.01, None)
fig.axes[0].set_xlim(0.05, 5.15)
fig.axes[0].set_ylim(10**(-1.6), 10**(4.4))
fig.axes[0].set_aspect(1.0 / fig.axes[0].get_data_ratio(), adjustable='box')


####  ax 2  ####
# data parameters
experiments = ['Stochastic'] * 3
freq_domains = [None] * 4
filter_widths = [0] * 4
helm_domains = [None, 'rot', 'div']
coor_shifts = [True] * 4
dims = [1] * 4
hours_cut = 0
depth_orientation = False
data_headers = map(DataHeader, experiments, filter_widths, freq_domains, helm_domains, coor_shifts, geostrophics, dims,
                   reduce_means)

# c = ['cornflowerblue', 'coral', 'tab:red', 'maroon', 'gold']
c = ['maroon', 'tab:red', 'tab:orange']
line_s = [1.6, 2, 2.5] * 4
s = ['-', (0, (5, 4)), (0, (1, 1)), (0, (4, 4))]

Title = ''
labels=[]
for freq_domain, filter_width, helm, coor_shift, exp in zip(freq_domains, filter_widths, helm_domains, coor_shifts,
                                                            experiments):
    labels.append('%s\ \ (%s)' % (filter_to_text(filter_width, freq_domain, helm, fw=True), exp_name[exp]))

figures = {}
# _title='$(b)$   Decomposed Freq. Spectra'
_title=''
figures[0] = {'title': _title, 'filters': data_headers, 'labels':labels}

psd_dict, freq_dict, fig = PSD_plot(figures, Title, hours_cut, depth_orientation, c, line_s, s, smoothing=1,
                                    plot_f=True, fig=fig, ax=fig.axes[1])

fig.axes[1].set_xlim(0.05, 5.15)
fig.axes[1].set_ylim(10**(-1.6), 10**(4.4))
fig.axes[1].set_aspect(1.0 / fig.axes[1].get_data_ratio(), adjustable='box')
fig.axes[1].legend(loc=3, fontsize=sizes.legend_size, handlelength=1.6)
# fig.axes[0].set_title(title, fontsize=title_size, fontweight=title_weight)


####  ax 2  ####

# data parameters
experiments = ['Stochastic', 'Steady'] * 2
freq_domains = ['LF'] * 4
filter_widths = [24] * 4
helm_domains = ['rot'] * 2 + ['div'] *2
dims = [2] * 4
coor_shifts = [True] * 4
hours_cut = 0
data_headers = map(DataHeader, experiments, filter_widths, freq_domains, helm_domains, coor_shifts, geostrophics, dims,
                   reduce_means)


depth_orientation = False
c = ['maroon', 'royalblue', 'maroon', 'royalblue', 'green']
line_s = [2, 2, 2, 2]
s = [(0, (4, 4)), (0, (4, 4)), (0, (1, 1)), (0, (1, 1))]
Title = ''
labels=[]
for freq_domain, filter_width, helm, exp in zip(freq_domains, filter_widths, helm_domains, experiments):
    labels.append('%s\ \ (%s)' % (filter_to_text(filter_width, freq_domain, helm, fw=False), exp_name[exp]))

figures = {}
_title=''
# _title='$(c)$   Decomposed Wavenumber Spectra'
figures[0] = {'title': _title, 'filters': data_headers, 'labels':labels}

psd_dict, freq_dict, fig = PSD_plot(figures, Title, hours_cut, depth_orientation, c, line_s, s, smoothing=1, fig=fig,
                                    ax=fig.axes[2])

fig.axes[2].set_xlim(0.7, 230)
fig.axes[2].set_ylim(1e-10, 20)
fig.axes[2].set_aspect(1.0 / fig.axes[2].get_data_ratio(), adjustable='box')
fig.axes[2].legend(loc=3, fontsize=sizes.legend_size, handlelength=1.6)

for i, letter in enumerate('ab'):
    axes[i].text(x=3, y=7000, s=r'$\mathbf{(%s)}$' % letter, fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight)

axes[2].text(x=120, y=2, s=r'$\mathbf{(%s)}$' % 'c', fontsize=sizes.title_size + 3,
                 fontweight=sizes.title_weight)

plt.show()
save_fig(fig, title='PSD_coordinate_shift', format='jpg')
