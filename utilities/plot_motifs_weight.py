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

path = '../data/motifs/motifs_weights/generative/5/08/steep'
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

    if interval >= 21:# or interval <= 14:
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
limits_y_steep = {'M15': 4.2987502900671553, 'M3': 7.8972724550828648, 'M9': 7.147837355023527, 'M0': 8.0530439432653047, 'M4': 5.9830612244455352, 'M1': 9.4275811043785698, 'M8': 5.7100822887254381, 'M12': 7.0283429888893556, 'M16': 10.050544629696276, 'M13': 2.7189886804605292, 'M14': 5.0850547841377196, 'M20': 9.2170853417721545, 'M19': 20.426930939454287, 'M11': 7.1841775316819589, 'M7': 7.9072539792257421, 'M5': 10.715593272122231, 'M6': 9.0364294988516463, 'M17': 9.8198160389498685, 'M2': 8.0488419817160093, 'M10': 20.855825913361898, 'M18': 6.2422368124944878}

limits_y_inhib = {'M1': 10.252633802156769, 'M11': 9.9390502958172302, 'M3': 20.694111073943727, 'M5': 20.059685878677747, 'M9': 10.012212002035568, 'M12': 8.0273688046929319, 'M19': 20.737676430293448, 'M4': 7.1708109037725203, 'M2': 9.1369680065207231, 'M7': 9.1700362028079656, 'M0': 10.107869935993531, 'M14': 7.2909256631621684, 'M8': 8.8893947942588731, 'M10': 21.142869122808531, 'M13': 5.7262114967542956, 'M17': 10.739823384170785, 'M18': 9.2533418209668632, 'M6': 9.4987740902157078, 'M20': 9.7709782769547147, 'M15': 6.0611239954656417, 'M16': 10.755161821126348}

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

max_val_list = {}
for m in motif_count_interval:
    print(m)
    data_to_plot = []
    max_value_interval = []
    for idx in range(len(motif_count_interval[m])-1):
        if len(motif_count_interval[m][idx]) == 0:
            continue
        max_value_interval.append(max(motif_count_interval[m][idx]))

    max_val = max(max_value_interval)
    max_val_list[m] = max_val
    for idx in range(len(motif_count_interval[m])-1):
        if len(motif_count_interval[m][idx]) == 0:
            data_to_plot.append([])
            continue
        for idx_v in range(len(motif_count_interval[m][idx])):
            motif_count_interval[m][idx][idx_v] = motif_count_interval[m][idx][idx_v]
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

    dir_save = '../plots/motif_weights_plots/generative/11_14/v1/steep'
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
    plt.ylabel(r'\textbf{Motif weights}', fontsize=25)
    # limits_y_steep[m] = third_quartile + math.pow(10, int(math.log10(third_quartile)))
    plt.ylim([0, max(limits_y_steep[m], limits_y_inhib[m])])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()

print(limits_y_steep)
