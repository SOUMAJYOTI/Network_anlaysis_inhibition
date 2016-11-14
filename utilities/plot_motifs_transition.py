# Plot the transitions of motif patterns form 4 sized motifs to 5 sized motifs in
# next interval. Observe which transitions occur in which of the intervals.

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
import math

path = 'motifs_transition/10_05/5/08/inhib'
number_intervals = len([name for name in os.listdir(path)])
cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
motif_trans_interval = {}
titles = []

for filename in os.listdir(path):
    print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_transitions = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval
    if int_reverse < 50:
        continue

    # print(motif_transitions)
    try:
        for m4 in motif_transitions:
            if m4 == {}:
                continue
            # print(m4)
            if m4 not in motif_trans_interval:
                motif_trans_interval[m4] = {}
            for m5 in motif_transitions[m4]:
                # print(m5)
                if m5 not in motif_trans_interval[m4]:
                    motif_trans_interval[m4][m5] = [[] for i in range(50)]
                motif_trans_interval[m4][m5][int_reverse-50].extend(motif_transitions[m4][m5])
        #titles.append('I' + str(interval))
    except:
        continue

print('Saving plots...')
# limits_y_steep  = {'M4_4': {'M5_15': 1.0281842195934305, 'M5_10': 0.79024605017994332, 'M5_8': 0.37469370131155211, 'M5_13': 0.1122609241390131, 'M5_14': 0.83144279893940376, 'M5_17': 0.83556149732620333, 'M5_16': 0.93055353828774301, 'M5_9': 0.50008256274768825, 'M5_12': 0.45782431206531599, 'M5_11': 0.97050826331825313, 'M5_5': 1.2, 'M5_18': 1.2, 'M5_7': 0.59448580093741388}, 'M4_3': {'M5_15': 1.0071915685058896, 'M5_10': 0.48117578456318916, 'M5_12': 0.56319857036761978, 'M5_4': 0.39123376962530426, 'M5_13': 0.68369943660188648, 'M5_6': 0.60438905949933597, 'M5_17': 1.2, 'M5_14': 1.2, 'M5_16': 1.2, 'M5_1': 0.95211718564710091, 'M5_9': 0.70718967362225782, 'M5_8': 0.55323881396031571, 'M5_11': 0.47055665442215316, 'M5_5': 1.2, 'M5_2': 1.2, 'M5_3': 0.46503430210802049, 'M5_7': 0.88104057822995641}, 'M4_1': {'M5_15': 0.96243330478888089, 'M5_10': 0.3369285230097161, 'M5_13': 0.75646761392729145, 'M5_8': 0.39584956402400989, 'M5_4': 0.34537198742577713, 'M5_2': 0.11274762260925313, 'M5_6': 0.3342159543339569, 'M5_9': 0.37072664359861596, 'M5_14': 0.79843546284224254, 'M5_16': 1.2, 'M5_1': 0.11817981891771752, 'M5_12': 0.43847748399425629, 'M5_11': 0.57495393007199191, 'M5_5': 0.31272000434129843, 'M5_18': 1.2, 'M5_19': 1.2, 'M5_3': 0.3077790771781973, 'M5_7': 0.31575614641553323}, 'M4_2': {'M5_15': 0.88655684104627763, 'M5_10': 0.48058295714206101, 'M5_12': 0.47037904599659286, 'M5_4': 0.42885134433717431, 'M5_13': 0.33903617521775065, 'M5_6': 0.39473248259991972, 'M5_17': 0.55817202225303175, 'M5_14': 1.2, 'M5_16': 0.95946555992427562, 'M5_1': 0.847448060641669, 'M5_9': 0.42311527393368409, 'M5_8': 0.34053429820931247, 'M5_11': 0.53369616446349621, 'M5_5': 0.46514377401511636, 'M5_18': 1.2, 'M5_2': 0.71477759328847323, 'M5_3': 0.39577537351880476, 'M5_7': 0.4159959088161912}, 'M4_0': {'M5_15': 0.96218729955099414, 'M5_10': 0.8827609687868081, 'M5_8': 0.45069762250251144, 'M5_0': 0.10698710071336379, 'M5_4': 0.10785392466689443, 'M5_13': 0.97349593822890701, 'M5_6': 0.30247573526211458, 'M5_9': 0.58179864027194561, 'M5_14': 0.49674518396786271, 'M5_16': 1.2, 'M5_1': 0.11003896323903099, 'M5_7': 0.31622543925939922, 'M5_12': 0.54057044659136155, 'M5_11': 0.30321253554359007, 'M5_5': 0.06552048295040637, 'M5_18': 1.0534177439348327, 'M5_2': 0.4320547036218727, 'M5_3': 0.1150042023924318, 'M5_19': 1.2}}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

for m4 in motif_trans_interval:
    for m5 in motif_trans_interval[m4]:
        max_value_interval = []
        for idx in range(len(motif_trans_interval[m4][m5])):
            if len(motif_trans_interval[m4][m5][idx]) == 0:
                continue
            max_value_interval.append(max(motif_trans_interval[m4][m5][idx]))

        try:
            max_val = max(max_value_interval)
        except:
            continue
        data_to_plot = []
        for idx in range(len(motif_trans_interval[m4][m5])):
            for idx_v in range(len(motif_trans_interval[m4][m5][idx])):
                motif_trans_interval[m4][m5][idx][idx_v] = motif_trans_interval[m4][m5][idx][idx_v]
            data_to_plot.append(motif_trans_interval[m4][m5][idx])
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

        #ax.set_title('Motif transition:' + str(m4) + '-->' + str(m5))
        ax.set_ylabel('Interval')
        #ax.set_ylim([0, 0.05])
        ax.set_xticklabels(titles)

        third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
        third_quartile = max(third_quartile)

        dir_save = 'motif_transition_plots/v2/inhib'
        if not os.path.exists(dir_save):
            os.makedirs(dir_save)
        file_save = dir_save + '/' + 'mt_steep_' + str(m4) + '_' + str(m5) + '.png'

        try:
            plt.ylim([0, third_quartile + 0.2 * math.pow(10, int(math.log10(third_quartile)))])
        except:
            print(third_quartile)
        # plt.ylim([0, limits_y_steep[m4][m5]])
        major_ticks = np.arange(0, 41, 5)
        plt.xticks(major_ticks)
        plt.tick_params('y', labelsize=15)
        plt.tick_params('x', labelsize=15)
        plt.xlabel(r'\textbf{Intervals before inhib region}', fontsize=15)
        plt.ylabel(r'\textbf{Normalized motif counts}', fontsize=15)
        # if m4 not in limits_y_steep:
        #     limits_y_steep[m4] = {}
        # limits_y_steep[m4][m5] = third_quartile + 0.2 * math.pow(10, int(math.log10(third_quartile)))
        # plt.ylim([0, 0.3])
        plt.savefig(file_save)
        plt.close()

# print(limits_y_steep)

