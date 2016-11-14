import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import seaborn

import datetime
import pickle
import numpy as np
import matplotlib.image as image

path = 'motifs_weights/generative/5/08/steep'
cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
motif_count_interval = {}
titles = []

number_intervals = len([name for name in os.listdir(path)])
for filename in os.listdir(path):
    # print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_count = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval

    if interval >= 51:
        continue
    for m in motif_count:
        if m not in motif_count_interval:
            motif_count_interval[m] = [[] for i in range(50)]
        list_filtered = []
        for v in motif_count[m]:
            if v == inf:
                continue
            if v == 0:
                continue
            # v = math.log(v)
            list_filtered.append(log(v))
        motif_count_interval[m][int_reverse-50].extend(list_filtered)
    # titles.append('I' + str(interval))

print('Saving plots...')
limits_y_steep = {'M16': 0.8682747128042938, 'M18': 0.94771363847784729, 'M11': 0.86327783579165485, 'M1': 0.85165038303016249, 'M5': 1.0334808251552958, 'M2': 0.81972941671463451, 'M8': 0.74308752071433237, 'M3': 0.77510075753449281, 'M12': 0.85856400583588344, 'M7': 0.88446815076679552, 'M0': 0.77690135385164738, 'M15': 0.53695488569176419, 'M20': 0.92229001406757116, 'M19': 1.0288987355433072, 'M17': 0.88731435954200522, 'M10': 1.0770096915739313, 'M14': 0.99241453959280257, 'M4': 0.69631949993558395, 'M6': 0.80042017367611429, 'M9': 0.72430682229544563, 'M13': 0.90886885259377936}

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
        box.set(color='#7570b3', linewidth=2)
        # change fill color
        box.set(facecolor='#1b9e77')

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    #ax.set_ylim([0, 0.05])

    third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
    third_quartile = max(third_quartile)

    dir_save = 'motif_weights_plots/generative/v4/5/steep'
    if not os.path.exists(dir_save):
        os.makedirs(dir_save)
    file_save = dir_save + '/' + 'mw_steep_' + str(m) + '.png'
    plt.ylim([0, third_quartile + 0.2*math.pow(10, int(math.log10(third_quartile)))])
    plt.ylim([0, limits_y_steep[m]])
    major_ticks = np.arange(0, 51, 5)
    plt.xticks(major_ticks)
    plt.tick_params('y', labelsize=15)
    plt.tick_params('x', labelsize=15)
    plt.xlabel(r'\textbf{Intervals before steep region}', fontsize=15)
    plt.ylabel(r'\textbf{Normalized motif weights}', fontsize=15)
    # limits_y_steep[m] = third_quartile + 0.2*math.pow(10, int(math.log10(third_quartile)))
    # plt.ylim([0, y])
    plt.savefig(file_save)
    plt.close()

# print(limits_y_steep)