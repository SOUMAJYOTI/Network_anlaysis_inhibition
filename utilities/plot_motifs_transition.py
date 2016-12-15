# Plot the transitions of motif patterns form 4 sized motifs to 5 sized motifs in
# next interval. Observe which transitions occur in which of the intervals.

import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle
import math

path = '../data/motifs/motifs_transition/10_05/5/08/inhib'
number_intervals = 21
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
    if interval> 20 or interval <=0:
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
                    motif_trans_interval[m4][m5] = [[] for i in range(20)]
                motif_trans_interval[m4][m5][20-int_reverse].extend(motif_transitions[m4][m5])
        #titles.append('I' + str(interval))
    except:
        continue

print('Saving plots...')
limits_y_steep = {'M4_1': {'M5_8': 6352.25, 'M5_18': 3076.0, 'M5_9': 4025.0, 'M5_4': 6997.5, 'M5_10': 3592.75, 'M5_14': 2662.0, 'M5_7': 6523.0, 'M5_5': 2832.5, 'M5_19': 76.0, 'M5_1': 5675.0, 'M5_20': 952.0, 'M5_16': 3342.0, 'M5_12': 3848.5, 'M5_17': 4734.0, 'M5_3': 6538.75, 'M5_2': 3497.5, 'M5_15': 7950.0, 'M5_11': 2318.0, 'M5_13': 3038.0, 'M5_6': 5214.75}, 'M4_3': {'M5_8': 6474.0, 'M5_18': 2696.0, 'M5_14': 9411.0, 'M5_4': 30431.75, 'M5_10': 4106.75, 'M5_9': 4134.5, 'M5_7': 3552.0, 'M5_5': 32209.0, 'M5_19': 660.25, 'M5_1': 29817.0, 'M5_20': 0, 'M5_16': 2833.0, 'M5_12': 3384.0, 'M5_17': 2210.5, 'M5_3': 9467.0, 'M5_2': 22940.0, 'M5_15': 2513.0, 'M5_11': 5125.5, 'M5_13': 4432.0, 'M5_6': 6019.5}, 'M4_5': {'M5_8': 326.5, 'M5_18': 2280.0, 'M5_14': 2780.0, 'M5_10': 540.0, 'M5_9': 223.25, 'M5_7': 4898.0, 'M5_6': 515.5, 'M5_19': 2415.0, 'M5_5': 337.0, 'M5_16': 2827.0, 'M5_12': 3811.5, 'M5_17': 548.0, 'M5_15': 821.0, 'M5_11': 2688.0, 'M5_13': 2132.0}, 'M4_2': {'M5_8': 28831.0, 'M5_18': 3942.0, 'M5_9': 5612.5, 'M5_4': 8291.5, 'M5_10': 4519.5, 'M5_14': 5650.5, 'M5_7': 7715.0, 'M5_5': 3139.75, 'M5_19': 3297.0, 'M5_1': 8448.5, 'M5_16': 4377.75, 'M5_12': 5915.5, 'M5_17': 5104.0, 'M5_3': 7662.0, 'M5_2': 25174.0, 'M5_15': 4664.0, 'M5_11': 4790.0, 'M5_13': 6423.5, 'M5_6': 7243.25}, 'M4_0': {'M5_8': 5265.0, 'M5_18': 2828.0, 'M5_9': 3343.5, 'M5_4': 7294.25, 'M5_0': 7035.25, 'M5_10': 2377.5, 'M5_14': 3164.0, 'M5_7': 5415.5, 'M5_5': 2940.0, 'M5_19': 67.0, 'M5_6': 4970.0, 'M5_16': 3841.0, 'M5_12': 1039.25, 'M5_17': 680.5, 'M5_15': 20752.0, 'M5_3': 8434.0, 'M5_2': 54371.25, 'M5_1': 6794.5, 'M5_11': 874.75, 'M5_13': 2885.5, 'M5_20': 326.0}, 'M4_4': {'M5_8': 4150.0, 'M5_18': 2062.0, 'M5_14': 3289.0, 'M5_10': 3606.25, 'M5_9': 3285.0, 'M5_7': 5674.25, 'M5_6': 2436.0, 'M5_5': 2995.25, 'M5_16': 5043.75, 'M5_12': 5207.0, 'M5_17': 3022.75, 'M5_15': 2596.75, 'M5_11': 4539.0, 'M5_13': 4036.5}}

limits_y_inhib = {'M4_0': {'M5_0': 8405.5, 'M5_9': 3060.75, 'M5_2': 22222.5, 'M5_1': 8026.75, 'M5_7': 5801.25, 'M5_11': 2052.75, 'M5_18': 4414.25, 'M5_6': 5230.75, 'M5_16': 5028.0, 'M5_13': 2421.0, 'M5_20': 0, 'M5_19': 0, 'M5_5': 2707.75, 'M5_14': 3883.5, 'M5_10': 2205.0, 'M5_8': 4116.0, 'M5_12': 2361.0, 'M5_3': 9176.5, 'M5_4': 8359.0, 'M5_15': 20752.0, 'M5_17': 3298.0}, 'M4_1': {'M5_9': 5453.5, 'M5_2': 3636.5, 'M5_1': 6236.25, 'M5_7': 6924.75, 'M5_11': 4959.75, 'M5_18': 4448.5, 'M5_6': 5892.0, 'M5_16': 2969.0, 'M5_13': 5260.75, 'M5_19': 0, 'M5_20': 0, 'M5_5': 3420.75, 'M5_14': 2972.5, 'M5_10': 3096.25, 'M5_8': 6774.0, 'M5_12': 3449.25, 'M5_3': 7310.0, 'M5_4': 7359.0, 'M5_15': 7950.0, 'M5_17': 5083.0}, 'M4_3': {'M5_2': 8915.0, 'M5_16': 4284.25, 'M5_7': 4331.75, 'M5_11': 2605.75, 'M5_18': 2851.5, 'M5_6': 5033.0, 'M5_3': 7968.5, 'M5_1': 20421.0, 'M5_13': 4930.5, 'M5_20': 328.0, 'M5_19': 0, 'M5_5': 10799.75, 'M5_14': 3393.25, 'M5_10': 3038.5, 'M5_8': 3190.75, 'M5_12': 2723.25, 'M5_9': 3245.0, 'M5_4': 21848.5, 'M5_15': 4098.0, 'M5_17': 2636.0}, 'M4_2': {'M5_2': 28400.5, 'M5_16': 21687.0, 'M5_7': 8935.0, 'M5_11': 5288.75, 'M5_18': 4713.25, 'M5_6': 8492.0, 'M5_3': 9291.0, 'M5_1': 33021.5, 'M5_13': 6625.0, 'M5_19': 0, 'M5_5': 3766.0, 'M5_14': 5144.0, 'M5_10': 4587.75, 'M5_8': 7851.0, 'M5_12': 6400.5, 'M5_9': 6364.0, 'M5_4': 9305.75, 'M5_15': 4277.75, 'M5_17': 5579.0}, 'M4_4': {'M5_2': 0, 'M5_7': 6277.0, 'M5_11': 4011.75, 'M5_18': 2025.75, 'M5_6': 2602.5, 'M5_16': 2734.0, 'M5_13': 3277.0, 'M5_20': 0, 'M5_19': 0, 'M5_5': 21091.5, 'M5_14': 2353.25, 'M5_10': 4247.0, 'M5_8': 3298.0, 'M5_12': 2862.5, 'M5_9': 4366.25, 'M5_15': 2740.25, 'M5_17': 2220.25}, 'M4_5': {'M5_7': 6822.5, 'M5_8': 897.25, 'M5_18': 922.0, 'M5_6': 2673.75, 'M5_16': 4098.5, 'M5_13': 6723.0, 'M5_19': 2085.75, 'M5_5': 423.0, 'M5_14': 10904.5, 'M5_10': 456.5, 'M5_11': 3933.5, 'M5_12': 5242.0, 'M5_9': 2222.0, 'M5_15': 2371.25, 'M5_17': 2839.0}}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

for m4 in motif_trans_interval:
    for m5 in motif_trans_interval[m4]:

        data_to_plot = []
        for idx in range(len(motif_trans_interval[m4][m5])-1, -1, -1):
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

        #ax.set_title('Motif transition:' + str(m4) + '-->' + str(m5))
        ax.set_ylabel('Interval')
        #ax.set_ylim([0, 0.05])
        ax.set_xticklabels(titles)

        third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
        third_quartile = max(third_quartile)

        dir_save = '../plots/motif_transition_plots/12_08/inhib'
        if not os.path.exists(dir_save):
            os.makedirs(dir_save)
        file_save = dir_save + '/' + 'mt_inhib_' + str(m4) + '_' + str(m5) + '.png'

        # try:
        #     plt.ylim([0, third_quartile + 0.3 * math.pow(10, int(math.log10(third_quartile)))])
        # except:
        #     print(third_quartile)
        # plt.ylim([0, limits_y_steep[m4][m5]])
        # major_ticks = np.arange(0, 21, 5)
        # plt.xticks(major_ticks)
        plt.tick_params('y', labelsize=25)
        plt.tick_params('x', labelsize=25)
        plt.xlabel(r'\textbf{Intervals before inhibition region}', fontsize=25)
        plt.ylabel(r'\textbf{Motif transition count}', fontsize=25)

        # if m4 not in limits_y_steep:
        #     limits_y_steep[m4] = {}
        # try:
        #     limits_y_steep[m4][m5] = third_quartile + math.pow(10, int(math.log10(third_quartile)))
        # except:
        #     limits_y_steep[m4][m5] = 0
        try:
            plt.ylim([0, max(limits_y_steep[m4][m5], limits_y_inhib[m4][m5])])
        except:
            pass
        plt.grid(True)
        plt.savefig(file_save)
        plt.close()

print(limits_y_steep)

