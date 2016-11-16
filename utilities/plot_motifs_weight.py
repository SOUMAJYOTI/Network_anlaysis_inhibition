import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle
import numpy as np
import matplotlib.image as image

path = '../data/motifs/motifs_weights/generative/5/08/v2/steep'
cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
motif_count_interval = {}
titles = []

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

    if interval >= 21 or interval == 0:
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
            list_filtered.append(log(v))
        motif_count_interval[m][int_reverse-20].extend(list_filtered)
    # titles.append('I' + str(interval))

print('Saving plots...')
limits_y_steep = {'M20': 0.86106343302304023, 'M7': 0.64438863366874344, 'M0': 0.52041106724605668, 'M14': 0.82361752469746308, 'M1': 0.68151792239231535, 'M16': 0.61215238307654141, 'M11': 0.71258912196837776, 'M18': 0.75033299173532964, 'M15': 0.63653643105590496, 'M12': 0.57797418489392127, 'M10': 0.85805878338268959, 'M9': 0.64632714590710494, 'M17': 0.82255026770471429, 'M13': 0.7994586223008413, 'M6': 0.73689229271002721, 'M5': 0.91587583606713596, 'M19': 0.56613519121015421, 'M4': 0.76217142193461984, 'M8': 0.74789042769091307, 'M2': 0.67967976102171346, 'M3': 0.82927241552399344}

limits_y_inhib = {'M4': 0.6119413589629128, 'M17': 0.86972305093836666, 'M20': 0.97369384484410837, 'M13': 0.41529521258120816, 'M11': 0.66327600104024076, 'M0': 0.76435835449444878, 'M8': 0.55645140340386767, 'M9': 0.59648458010303251, 'M1': 0.88633079899268052, 'M2': 0.70509815145500632, 'M12': 0.71837095061478551, 'M19': 0.95843817663788955, 'M3': 0.7034706865702518, 'M5': 0.94900075671069617, 'M6': 0.79633182558263638, 'M14': 0.65208384262460695, 'M16': 0.87655907654473819, 'M7': 0.6152942772733162, 'M10': 1.0384931884158897, 'M18': 0.6588883159225376, 'M15': 0.51821996407136606}

# limits_y_steep = {'M0': 0.37004071345856793, 'M12': 0.46211606187440929, 'M18': 0.05304876869887281,
#                   'M5': 0.078728494978215019, 'M4': 0.33054649586121532, 'M8': 0.43699005161852106,
#                   'M1': 0.59111832024447453, 'M11': 0.46045738089362498, 'M20': 0.38708133964946784,
#                   'M13': 0.083492146369315443, 'M7': 0.62113487902589126, 'M15': 0.40014426075444243,
#                   'M9': 0.30937759738947024, 'M3': 0.10585969672779087, 'M16': 0.42398788129952036,
#                   'M2': 0.51231722060716445, 'M17': 0.49101915684395109, 'M19': 0.3370070481280063,
#                   'M6': 0.36705330453365875, 'M10': 0.35369554616907251, 'M14': 0.32429041419793164}

# limits_y_steep = {}
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

for m in motif_count_interval:
    print(m)
    data_to_plot = []
    max_value_interval = []
    for idx in range(len(motif_count_interval[m])-1):
        if len(motif_count_interval[m][idx]) == 0:
            continue
        max_value_interval.append(max(motif_count_interval[m][idx]))

    max_val = max(max_value_interval)
    for idx in range(len(motif_count_interval[m])-1):
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

    dir_save = '../plots/motif_weights_plots/generative/11_14/v2/steep'
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    file_save = dir_save + '/' + 'mw_steep_' + str(m) + '.png'
    # plt.ylim([0, third_quartile + 0.3*math.pow(10, int(math.log10(third_quartile)))])
    # plt.ylim([0, limits_y_steep[m]])
    # major_ticks = np.arange(0, 21, 5)
    # plt.xticks(major_ticks)
    plt.tick_params('y', labelsize=25)
    plt.tick_params('x', labelsize=25)
    plt.xlabel(r'\textbf{Intervals before steep region}', fontsize=25)
    plt.ylabel(r'\textbf{Normalized motif weights}', fontsize=25)
    # limits_y_steep[m] = third_quartile + 0.2*math.pow(10, int(math.log10(third_quartile)))
    plt.ylim([0, limits_y_steep[m]])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()

print(limits_y_steep)