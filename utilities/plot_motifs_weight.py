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

number_intervals = 20
for filename in os.listdir(path):
    # print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_count = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval

    # if interval > 30 or interval <= 10:
    #     continue
    if interval <=0 or interval > 20:
        continue
    # print(int_reverse, interval)
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
        motif_count_interval[m][20-interval].extend(list_filtered)
    # titles.append('I' + str(interval))

print('Saving plots...')
# limits_y_steep = {'M15': 4.2987502900671553, 'M3': 7.8972724550828648, 'M9': 7.147837355023527, 'M0': 8.0530439432653047, 'M4': 5.9830612244455352, 'M1': 9.4275811043785698, 'M8': 5.7100822887254381, 'M12': 7.0283429888893556, 'M16': 10.050544629696276, 'M13': 2.7189886804605292, 'M14': 5.0850547841377196, 'M20': 9.2170853417721545, 'M19': 20.426930939454287, 'M11': 7.1841775316819589, 'M7': 7.9072539792257421, 'M5': 10.715593272122231, 'M6': 9.0364294988516463, 'M17': 9.8198160389498685, 'M2': 8.0488419817160093, 'M10': 20.855825913361898, 'M18': 6.2422368124944878}

# limits_y_inhib = {'M1': 10.252633802156769, 'M11': 9.9390502958172302, 'M3': 20.694111073943727, 'M5': 20.059685878677747, 'M9': 10.012212002035568, 'M12': 8.0273688046929319, 'M19': 20.737676430293448, 'M4': 7.1708109037725203, 'M2': 9.1369680065207231, 'M7': 9.1700362028079656, 'M0': 10.107869935993531, 'M14': 7.2909256631621684, 'M8': 8.8893947942588731, 'M10': 21.142869122808531, 'M13': 5.7262114967542956, 'M17': 10.739823384170785, 'M18': 9.2533418209668632, 'M6': 9.4987740902157078, 'M20': 9.7709782769547147, 'M15': 6.0611239954656417, 'M16': 10.755161821126348}

limits_y_steep = {'M2': 26.552467589070819, 'M0': 26.011674647737504, 'M1': 31.128221181054286, 'M9': 27.172730776770639, 'M7': 23.844885259585723, 'M14': 65.705507083203727, 'M18': 49.613637514464152, 'M13': 42.000487614752714, 'M17': 49.102141715335868, 'M15': 34.403804541829231, 'M5': 72.631696943965181, 'M11': 42.072715434211084, 'M3': 21.931672987458967, 'M4': 41.139613123484565, 'M19': 23.463582241225414, 'M8': 30.316590172798094, 'M6': 23.744302643009199, 'M20': 49.836840975673674, 'M10': 10, 'M16': 23.338865690931755, 'M12': 29.724444181917896}
limits_y_inhib = {'M16': 58.271285477583653, 'M0': 34.861449909515116, 'M19': 54.605892338662798, 'M5': 86.695668182952375, 'M15': 64.01536233716827, 'M10': 10.32291567858247, 'M4': 53.869426642356174, 'M1': 53.091588220813264, 'M11': 47.826803985518303, 'M2': 50.349709928717111, 'M13': 53.152852333779066, 'M8': 56.73824231042822, 'M12': 61.442760394701757, 'M18': 60.568553725368012, 'M20': 51.345523668172191, 'M3': 211.00124343278313, 'M7': 56.917226318544408, 'M17': 61.594067750623331, 'M6': 47.640103963820508, 'M14': 71.828181551760736, 'M9': 39.895771953649088}

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
    for idx in range(len(motif_count_interval[m])):
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

    dir_save = '../plots/motif_weights_plots/generative/12_08/v1/steep'
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
    # try:
    #     limits_y_steep[m] = third_quartile + math.pow(10, int(math.log10(third_quartile)))
    # except:
    #     pass
        # limits_y_steep[m] = 10

    try:
        plt.ylim([0, max(limits_y_steep[m], limits_y_inhib[m])])
    except:
        continue
    # plt.ylim([0, 5000])
    plt.grid(True)
    plt.savefig(file_save)
    plt.close()

print(limits_y_steep)
