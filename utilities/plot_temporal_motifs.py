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

path = '../data/motifs/frontier_motifs_count/5/frontiers_lt/inhib'
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

    if interval >= 21 or interval <= 0:
        continue
    print(20 - int_reverse)
    for m in motif_count:
        if m not in motif_count_interval:
            motif_count_interval[m] = [[] for i in range(20)]
        list_filtered = []
        for v in motif_count[m]:
            if v == inf:
                continue
            if v == 0:
                continue
            # v = math.log(v)
            list_filtered.append(v)
        motif_count_interval[m][20-int_reverse].extend(list_filtered)
    # titles.append('I' + str(interval))

print('Saving plots...')

limits_y_inhib_lt = {'M15': 41.0, 'M13': 545.0, 'M2': 36.0, 'M12': 2765.5, 'M10': 2883.0, 'M1': 29.0, 'M20': 10.0, 'M6': 53.0, 'M19': 22.0, 'M11': 229.0, 'M8': 2707.5, 'M7': 2931.5, 'M4': 2580.0, 'M17': 23.0, 'M18': 301.0, 'M9': 263.25, 'M5': 222.0, 'M14': 22.0, 'M3': 328.5, 'M16': 450.0, 'M0': 711.5}
limits_y_steep_lt = {'M3': 325.0, 'M19': 20.0, 'M4': 2759.0, 'M20': 10.0, 'M10': 2815.25, 'M6': 49.0, 'M14': 6.0, 'M17': 22.5, 'M15': 41.0, 'M12': 2922.0, 'M16': 454.0, 'M7': 2997.5, 'M11': 234.5, 'M18': 341.0, 'M9': 249.25, 'M2': 36.0, 'M5': 233.0, 'M0': 673.5, 'M8': 2900.0, 'M1': 34.0, 'M13': 737.0}

limits_y_inhib_gt = {'M9': 107.0, 'M0': 387.0, 'M19': 5.0, 'M20': 5.0, 'M8': 2551.0, 'M3': 226.0, 'M11': 94.5, 'M1': 50.25, 'M5': 57.0, 'M15': 24.0, 'M17': 5.5, 'M13': 372.0, 'M18': 254.25, 'M14': 4.25, 'M7': 2736.0, 'M16': 370.5, 'M2': 29.0, 'M10': 2543.0, 'M12': 2769.0, 'M6': 30.0, 'M4': 2576.5}
limits_y_steep_gt = {'M9': 66.0, 'M13': 220.0, 'M19': 4.0, 'M18': 278.25, 'M3': 237.0, 'M14': 8.75, 'M12': 2245.75, 'M7': 2595.0, 'M20': 10.5, 'M5': 59.0, 'M1': 375.5, 'M11': 66.25, 'M2': 25.5, 'M6': 39.0, 'M17': 5.0, 'M10': 2336.0, 'M15': 22.5, 'M4': 1079.25, 'M8': 923.0, 'M0': 313.25, 'M16': 310.5}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

limits_y_steep = {}
for m in motif_count_interval:
    print(m)
    data_to_plot = []
    max_value_interval = []
    for idx in range(len(motif_count_interval[m])):
        if len(motif_count_interval[m][idx]) == 0:
            continue
        max_value_interval.append(max(motif_count_interval[m][idx]))

    if max_value_interval == []:
        continue
    max_val = max(max_value_interval)
    if max_val == 0:
        continue
    for idx in range(len(motif_count_interval[m])-1, -1, -1):
        if len(motif_count_interval[m][idx]) == 0:
            data_to_plot.append([])
            continue
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

    dir_save = '../plots/temporal_motif_count_plots/12_08/frontiers_lt/v2/inhib'
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    file_save = dir_save + '/' + 'tmcl_inhib_' + str(m) + '.png'
    # try:
    #     plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
    # except:
    #     pass
    # plt.ylim([0, limits_y_steep[m]])
    # major_ticks = np.arange(0, 51, 5)
    # plt.xticks(major_ticks)
    plt.tick_params('y', labelsize=25)
    plt.tick_params('x', labelsize=25)
    plt.xlabel(r'\textbf{Network subsequences leading to $N_{inhib}$}', fontsize=25)
    plt.ylabel(r'\textbf{Motif counts}', fontsize=25)
    try:
        # limits_y_steep[m] = third_quartile + math.pow(10, int(math.log10(third_quartile)))
        plt.ylim([0, max(limits_y_steep_lt[m], limits_y_inhib_lt[m])])
    except:
        pass
    # plt.ylim([0, limits_y_steep[m]])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()

print(limits_y_steep)