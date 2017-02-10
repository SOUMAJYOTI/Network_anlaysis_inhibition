import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import math

if __name__ == "__main__":
    statistic_values = pickle.load(open('..//data//causality//11_15//statistics.pickle', 'rb'))
    crit_values = pickle.load(
        open('..//data//causality//11_15//critical.pickle', 'rb'))
    p_values = pickle.load(
        open('..//data//causality//11_15//p_values.pickle',  'rb'))
    cause_count = pickle.load(
        open('..//data//causality//11_15//cause_count.pickle', 'rb'))

    width = 0.35  # the width of the bars

    map_cent = {'cc': 'Clustering coeff.', 'entropy': 'Degree entropy', 'nbr_deg':'Nodal degree', 'bw': 'Betweenness',
                'pr': 'Pagerank'}
    data_mean = []
    # data_std = []
    titles = []
    for m in map_cent:
        sub = m
        temp = []
        num_feat = len(sub.split('+'))
        if num_feat != 1:
            continue
        for idx in range(len(p_values[sub])):
            if statistic_values[sub][idx] <= 0:
                continue
            temp.append(p_values[sub][idx])
        data_mean.append(np.mean(temp))
        # data_std.append(np.std(temp))
        titles.append(map_cent[m])

    print(data_mean)
    cause_percentage = []
    idx = 0
    cause_idx = []
    t = []
    for sub in statistic_values:
        num_feat = len(sub.split('+'))
        if num_feat != 1:
                continue

        cause_percentage.append(cause_count[sub]/len(statistic_values[sub])*100 )
        # print(sub, cause_count[sub]/len(statistic_values[sub])*100)
        # cause_idx.append(idx)
        print(sub, cause_count[sub]/len(statistic_values[sub])*100 )
        t.append(map_cent[sub])
        idx += 1

    # print(data_mean)
    print(t)
    hfont = {'fontname': 'Arial'}
    ind = np.arange(len(cause_percentage))
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)

    rects1 = ax.bar(ind, cause_percentage, width,
                    color='#0000ff')    # axes and labels
    ax.set_xlim(-width,len(ind)+width)
    ax.set_ylim(87, 95)
    ax.set_ylabel('Percentage', size=40, **hfont)
    ax.set_xticks(ind)
    xtickNames = ax.set_xticklabels(t, **hfont)
    plt.setp(xtickNames, rotation=45, fontsize=5)
    plt.grid(True)
    plt.xticks(size=25)
    plt.yticks(size=25)
    plt.subplots_adjust(left=0.13, bottom=0.30, top=0.9)
    ## add a legend
    # ax.legend( (rects1[0], ('Men', 'Women') )

    plt.show()
    plt.close()

    ind = np.arange(len(data_mean))  # the x locations for the groups
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    ## the bars
    # rects1 = ax.bar(ind, data_mean, width,
    #                 color='#C0C0C0')
    rects1 = ax.bar(ind, data_mean, width,
                    color='#0000ff')  # axes and labels
    ax.set_xlim(-width, len(ind) + width)
    ax.set_ylim(0, 1)
    ax.set_ylabel('Average P-values', size=40, **hfont)
    # ax.set_xlabel('Node-centric measures', size=30)
    # ax.set_title('Scores by group and gender')
    xTickMarks = titles
    ax.set_xticks(ind)
    xtickNames = ax.set_xticklabels(xTickMarks, **hfont)
    plt.setp(xtickNames, rotation=45, fontsize=5)
    plt.grid(True)
    plt.xticks(size=25)
    plt.yticks(size=25)
    plt.subplots_adjust(left=0.13, bottom=0.30, top=0.9)
    ## add a legend
    # ax.legend( (rects1[0], ('Men', 'Women') )

    plt.show()