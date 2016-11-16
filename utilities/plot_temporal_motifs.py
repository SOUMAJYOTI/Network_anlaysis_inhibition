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

number_intervals = 35
for filename in os.listdir(path):
    # print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_count = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval

    if interval >= 35 or interval <= 14:
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

limits_y_inhib_gt = {'M4': 1.1528002281734042, 'M10': 1.1658064813847313, 'M17': 0.58487546239074018, 'M11': 0.83837992037728637, 'M16': 1.003631400372563, 'M2': 0.68715548440355212, 'M19': 0.58685280723454158, 'M9': 0.86078993385998959, 'M13': 0.9007331403204657, 'M8': 1.1414792873978443, 'M7': 1.1783216855843752, 'M1': 0.68318044580940862, 'M3': 0.89093642185655453, 'M5': 0.82862802139849645, 'M20': 0.50658097730535623, 'M15': 0.66063050799552014, 'M14': 0.59073345008644307, 'M18': 0.86282224781086114, 'M0': 0.94453133338449247, 'M12': 1.1729819450718328, 'M6': 0.78197872519301725}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

limits_y_steep = {}
for m in motif_count_interval:
    print(m)
    data_to_plot = []
    max_value_interval = []
    for idx in range(len(motif_count_interval[m])-1):
        if len(motif_count_interval[m][idx]) == 0:
            continue
        max_value_interval.append(max(motif_count_interval[m][idx]))

    if max_value_interval == []:
        continue
    max_val = max(max_value_interval)
    if max_val == 0:
        continue
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

    dir_save = '../plots/temporal_motif_count_plots/11_16/frontiers_lt/v1/inhib'
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
    plt.xlabel(r'\textbf{Intervals before inhib region}', fontsize=25)
    plt.ylabel(r'\textbf{Normalized motif counts}', fontsize=25)
    try:
        limits_y_steep[m] = third_quartile + 0.3*math.pow(10, int(math.log10(third_quartile)))
        plt.ylim([0, limits_y_steep[m]])
    except:
        pass
    # plt.ylim([0, limits_y_steep[m]])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()

print(limits_y_steep)