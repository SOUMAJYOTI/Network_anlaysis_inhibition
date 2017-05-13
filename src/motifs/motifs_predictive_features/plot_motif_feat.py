from graph_tool.all import *
import graph_tool as gt
import graph_tool.stats as gts
import graph_tool.util as gtu
import graph_tool.draw as gtd
from pylab import *
from math import *
from numpy.random import *
import pickle
import pandas as pd
import os
import glob
import csv
import statistics as st
import random
import graph_tool.topology as gtt
import matplotlib.pyplot as plt
# import seaborn
import operator


def plot_line(x, y, x_title, legends):
    fig = plt.figure()
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    plt.xticks(x, x_title, size=30, rotation=50)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)
    for i in range(len(y)):
        ax.plot(x, y[i],  marker='o', linestyle='-', label=legends[i], linewidth=3)

    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    # plt.title('Network_cover vs Motif edge density ', color='#000000', weight="bold", size=50)
    plt.ylabel('% of edges covered by motif', size=40)
    plt.xlabel('Motif edge density (increasing -->)', size=50)
    plt.ylim([0, 100])
    ax.legend(loc='upper right', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)

    plt.show()


def plot_box(data_to_plot, titles):
    fig = plt.figure(1, figsize=(12, 8))

    # Create an axes instance
    ax = fig.add_subplot(111)

    # label = ['Database', 'Grade Changes', 'Phone']
    # Create the boxplot
    bp = ax.boxplot(data_to_plot, patch_artist=True)

    # colors_face = ['white', 'white', '#DCDCDC', '#DCDCDC', '#696969', '#696969']
    # hatch_pattern = ['|X', '|X', '', '', '', '']
    idx = 0
    for box in bp['boxes']:
        # change outline color
        box.set(color='#FFFFFF', linewidth=4)
        # change fill color
        # box.set(facecolor=colors_face[idx])
        # box.set(hatch=hatch_pattern[idx])
        idx += 1

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#000000', linewidth=4)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#000000', linewidth=4)

    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set(color='#000000', linewidth=6)

    # change the style of fliers    and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)

    third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
    third_quartile = max(third_quartile)

    first_quartile = [item.get_ydata()[1] for item in bp['whiskers']]
    first_quartile = max(first_quartile)

    hfont = {'fontname': 'Arial'}
    ax.set_ylabel('Amount (USD)', fontsize=50, **hfont)
    # plt.ylim([0, 2500])
    # plt.xlim([0, 5])
    plt.tick_params('y', labelsize=50)
    ax.set_xticklabels(titles, size=40, rotation=45, ha='right', **hfont)

    plt.show()


def plot_hist(data, label):
    plt.figure(figsize=(12, 8))
    hfont = {'fontname': 'Arial'}
    n, bins, patches = plt.hist(data, 30, lw=3, facecolor='b')
    # plt.yscale('log', nonposy='clip', basey=2)
    plt.xlabel(label, size=40, **hfont)
    plt.ylabel('Frequency', size=40, **hfont)
    # plt.title('Histogram of')
    # plt.xlim([0, 700])
    # plt.ylim([0, 2 ** (12)])
    plt.grid(True)
    plt.xticks(size=25)
    plt.yticks(size=25)
    # file_save = dir_save + '/' + 'count_motif_' + str(m) + '.png'
    # plt.savefig(file_save)
    plt.subplots_adjust(left=0.16, bottom=0.16)
    plt.savefig('cover_dist' + str(label) + '.png')
    plt.close()


def plot_bars(x, y, titles=[]):
    width = 1
    plt.bar(x, y, width, lw=2, color="blue")
    if len(titles) > 0:
        major_ticks = np.arange(0, len(titles), 3)
        labels = []
        for i in major_ticks:
            labels.append(str(titles[i])[:10])

        plt.xticks(major_ticks, labels, rotation=45, size=20, ha='center')
    else:
        plt.xticks(size=20)
    plt.yticks(size=20)
    plt.xlabel('Month-Year (Time)', size=25)
    plt.ylabel('Count of posts', size=25)
    # plt.title('Month-wise post counts', size=20)

    plt.subplots_adjust(left=0.13, bottom=0.25, top=0.95)
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    # Load the interim computed features
    X_counts = pickle.load(open('interim_feat/X_count.pickle', 'rb'))
    X_weights = pickle.load(open('interim_feat/X_weights.pickle', 'rb'))
    X_temporal = pickle.load(open('interim_feat/X_temporal.pickle', 'rb'))

    Y_cover = [pickle.load(open('interim_feat/Y_cover_min_1.pickle', 'rb')),
               pickle.load(open('interim_feat/Y_cover_min_4.pickle', 'rb'))]

    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))

    pattern_density = {}
    for idx_mp in dict_patterns:
        pattern_density[dict_patterns[idx_mp]] = len(list(idx_mp.edges())) / 10

    pattern_density_sort = sorted(pattern_density.items(), key=operator.itemgetter(1))
    mpatterns = [] # sorted by edge density patterns
    for p, d in pattern_density_sort:
       mpatterns.append(p)

    inv_dict_patterns = {v: k for k, v in dict_patterns.items()}
    # data_to_plot = []
    data_line = [[] for _ in range(len(Y_cover))]
    legends = ['Free nodes = 1', 'Free nodes > 1']
    x_title = []

    for idx_c in range(len(Y_cover)):
        print(len(Y_cover[idx_c]))
        data_mean = {}
        den = []
        for idx in range(len(mpatterns)):
            temp_val = []
            for mid in Y_cover[idx_c][mpatterns[idx]]:
                temp_val.append(Y_cover[idx_c][mpatterns[idx]][mid])
            # data_to_plot.append(temp_val)
            # if pattern_density[mpatterns[idx]] not in data_mean:
            #     data_mean[pattern_density[mpatterns[idx]]] = np.array(np.mean(temp_val))
            # else:
            #     if np.array(np.mean(temp_val)) >  data_mean[pattern_density[mpatterns[idx]]]:
            data_line[idx_c].append(np.array(np.mean(temp_val)))
            if idx_c == 0:
                x_title.append(mpatterns[idx])
            # den.append(pattern_density[mpatterns[idx]])

        # data_mean_sort = sorted(data_mean.items(), key=operator.itemgetter(0))
        # for d, v in data_mean_sort:
        #     den.append(d)
        #     data_line[idx_c].append(v)

    plot_line(list(range(len(mpatterns))), data_line, x_title, legends)
    # plot_box(data_to_plot, mpatterns)

    # plot the coverage
    # plot_hist(Y_cover, label)
