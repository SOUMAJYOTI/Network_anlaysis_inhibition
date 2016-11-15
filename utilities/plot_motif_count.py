from pylab import *
import datetime
import pickle
import numpy as np
import matplotlib.image as image
import matplotlib.pyplot as plt
import os
from matplotlib import rc

path = '../data/motifs/motifs_count/5/08/v1/inhib'

cnt_1 = 0
cnt_2 = 0
motif_count_interval = {}
titles = []

# number_intervals = len([name for name in os.listdir(path)])
number_intervals = 21
for filename in os.listdir(path):
    # print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_count = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval

    # print(int_reverse)
    if int_reverse < 0:
        continue
    for m in motif_count:
        if m not in motif_count_interval:
            motif_count_interval[m] = [[] for i in range(21)]
        list_filtered = []
        for v in motif_count[m]:
            if v == inf:
                continue
            if v == 0:
                continue
            # v = math.log(v)
            list_filtered.append(v)
        motif_count_interval[m][int_reverse-20].extend(list_filtered)
    # titles.append('I' + str(interval))

print('Saving plots...')

limits_y_inhib = {'M19': 0.00094817215452424166, 'M20': 0.05881844380403458, 'M2': 0.045579084845375248, 'M0': 0.41850331489814846, 'M16': 0.53090054538198905, 'M12': 0.10348765024435345, 'M8': 0.0879971181556196, 'M10': 0.091385099685204613, 'M5': 0.12023234829686442, 'M11': 0.066681715575620762, 'M1': 0.065612535612535605, 'M9': 0.44399181301993929, 'M14': 0.11685871090934383, 'M18': 0.10481751824817517, 'M7': 0.08319148936170212, 'M3': 0.098428750522356875, 'M13': 0.10423963903743316, 'M15': 0.094497716894977157, 'M4': 0.084761904761904761, 'M17': 0.47151162790697676, 'M6': 0.41162041001359995}
limits_y_steep = {'M19': 0.00083403417818740406, 'M6': 0.46843888454521759, 'M18': 0.11928571428571429, 'M2': 0.043956521739130436, 'M5': 0.40271084337349394, 'M12': 0.095994236311239198, 'M4': 0.11875739644970414, 'M1': 0.10183908045977011, 'M15': 0.10524752475247524, 'M14': 0.4173777101505417, 'M20': 0.11181818181818182, 'M0': 0.52477172544080597, 'M9': 0.5161979526109961, 'M17': 0.08747126436781609, 'M3': 0.10467372568332924, 'M7': 0.099047619047619051, 'M16': 0.597428549883116, 'M10': 0.12266227657572906, 'M8': 0.11538681948424069, 'M11': 0.10111111111111111, 'M13': 0.11090257023311416}

# limits_y_steep = {'M13': 0.10090257023311418, 'M5': 0.30271084337349397, 'M12': 0.091902017291066285,
#                   'M20': 0.12000000000000001, 'M6': 0.36843888454521762, 'M4': 0.10875739644970414,
#                   'M0': 0.42477172544080605, 'M11': 0.10000000000000001, 'M10': 0.11266227657572907,
#                   'M9': 0.41619795261099612, 'M17': 0.077471264367816095, 'M3': 0.094673725683329227,
#                   'M19': 0.00075203533026113665, 'M14': 0.31737771015054173, 'M16': 0.49742854988311597,
#                   'M8': 0.1053868194842407, 'M18': 0.10928571428571429, 'M1': 0.091839080459770114,
#                   'M7': 0.089047619047619056, 'M15': 0.095247524752475249, 'M2': 0.034391304347826092}

#
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

mpl.rcParams['text.latex.preamble'] = [
    r'\usepackage{siunitx}',  # i need upright \micro symbols, but you need...
    r'\sisetup{detect-all}',  # ...this to force siunitx to actually use your fonts
    r'\usepackage{helvet}',  # set the normal font here
    r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
    r'\sansmath'  # <- tricky! -- gotta actually tell tex to use!
]

for m in motif_count_interval:
    print(m)
    data_to_plot = []
    max_value_interval = []
    for idx in range(len(motif_count_interval[m])):
        if len(motif_count_interval[m][idx]) == 0:
            continue
        max_value_interval.append(max(motif_count_interval[m][idx]))

    max_val = max(max_value_interval)
    for idx in range(len(motif_count_interval[m])):
        if len(motif_count_interval[m][idx]) == 0:
            data_to_plot.append([])
            continue
        for idx_v in range(len(motif_count_interval[m][idx])):
            motif_count_interval[m][idx][idx_v] = motif_count_interval[m][idx][idx_v] / max_val
        data_to_plot.append(motif_count_interval[m][idx])
        titles.append(str(idx+1))

    # Create the box_plots
    fig = plt.figure(1, figsize=(10, 8))

    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    for box in bp['boxes']:
        # change outline color
        box.set(color='#0000FF', linewidth=2)
        # change fill color
        box.set(facecolor='#FFFFFF')

    ## change color and linewidth of the whiskers
    # for whisker in bp['whiskers']:
    #     whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    # for cap in bp['caps']:
    #     cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#FF0000', linewidth=4)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    #ax.set_ylim([0, 0.05])

    third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
    third_quartile = max(third_quartile)

    dir_save = '../plots/motif_count_plots/11_14/v1/inhib'
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    file_save = dir_save + '/' + 'mc_inhib_' + str(m) + '.png'
    # plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
    # plt.ylim([0, limits_y_steep[m]])
    # major_ticks = np.arange(0, 21, 5)
    # plt.xticks(major_ticks)
    plt.tick_params('y', labelsize=20)
    plt.tick_params('x', labelsize=20)
    plt.xlabel(r"\textbf{Intervals before inhibition region}", fontsize=25)
    plt.ylabel(r"\textbf{Normalized Motif counts}", fontsize=25)
    # limits_y_steep[m] = third_quartile + 0.3*math.pow(10, int(math.log10(third_quartile)))
    plt.ylim([0, limits_y_inhib[m]])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()
    # plt.show()

print(limits_y_steep)