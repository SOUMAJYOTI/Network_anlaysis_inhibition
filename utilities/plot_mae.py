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


titles = []
mae_model = pickle.load(open('..//..//data_files//results_granger//11_15//mae_OLS.pickle', 'rb'))

mae_new = {}
for sub in mae_model:
    for ms in mae_model[sub]:
        if ms not in mae_new:
            mae_new[ms] = {}
        if sub not in mae_new[ms]:
            mae_new[ms][sub] = []
        mae_new[ms][sub] = mae_model[sub][ms]

for ms in mae_new:
    if ms < 5:
        continue
    print(ms)
    data_to_plot = [[] for idx in range(5)]
    titles = []
    cnt_sub = 0
    for sub in mae_new[ms]:
        num_plus = len(sub.split('+'))
        # if num_plus != 1:
        #     continue
        for idx in range(len(mae_new[ms][sub])):
            val = mae_new[ms][sub][idx].tolist()[0]
            data_to_plot[cnt_sub].append(val)
        titles.append(sub)
        cnt_sub += 1

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

    ax.set_title('Absolute errors - '+ str(ms) +  '- step ahead', fontsize=30)
    # ax.set_xlabel('Features', fontsize=20)
    #ax.set_ylim([0, 0.05])
    ax.set_xticklabels(titles, size=20)

    # for tick in ax.get_xticklabels():
    #     tick.set_rotation(45)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(30)

    third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
    third_quartile = max(third_quartile)

    # # dir_save = 'motif_count_plots/5/inhib'
    # if not os.path.exists(dir_save):
    #     os.makedirs(dir_save)
    # file_save = dir_save + '/' + 'count_motif_' + str(m) + '.png'
    plt.ylim([0, third_quartile + 4*math.pow(10, int(math.log10(third_quartile)))])
    plt.subplots_adjust(bottom=0.15)
    # plt.savefig(file_save)
    plt.show()
    plt.close()

