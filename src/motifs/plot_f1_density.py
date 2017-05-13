from graph_tool.all import *
import graph_tool as gt
import graph_tool.stats as gts
import graph_tool.util as gtu
import graph_tool.draw as gtd
from pylab import *
import pickle
import graph_tool.topology as gtt
import matplotlib.pyplot as plt
import numpy as np
from sklearn.model_selection import KFold
import sklearn
from sklearn import linear_model
import random
from matplotlib.patches import Circle
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
import matplotlib.patches as patches
import operator

def plot_line(x, Y, title, X_title):
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    plt.xticks(x, X_title, size=30, rotation=50)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)
    for i in range(len(Y)):
        ax.plot(x, Y[i],  marker='o', linestyle='-', label='k=' + str(i+1), linewidth=3)

    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('F1 score', size=30)
    plt.xlabel('Motif ID', size=30)

    ax.legend(loc='lower left', fancybox=True, shadow=True, fontsize=20)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    # plt.ylim([0.05, 1.05])
    # plt.xlim([1.5, 10.5])
    # plt.savefig('v1/f1_plots_RF/Motif_' + m + '.png')
    # plt.savefig('v1/f1_plots_RF/Motif_' + m + '.pdf')
    # plt.close()
    plt.show()


if __name__ == "__main__":
    X_string = ['Motif counts', 'Motif transitions', 'Motif weights', 'Temporal motifs', 'All']
    # X_string = ['Motif counts', 'Motif transitions', 'Temporal motifs', 'All']

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']
    motif_imgs_dir = "motif_patterns/"
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    avg_precision = pickle.load(open('../05/v1/min_4/avg_precision_RF.pickle', 'rb'))
    avg_recall = pickle.load(open('../05/v1/min_4/avg_recall_RF.pickle', 'rb'))
    avg_f1 = pickle.load(open('../05/v1/min_4/avg_f1_RF.pickle', 'rb'))
    # avg_random = pickle.load(open('../05/v1/avg_random.pickle', 'rb'))

    pattern_density = {}
    for idx_mp in dict_patterns:
        # max_deg = -1
        # for v in idx_mp.vertices():
        #     if len(list(v.out_edges())) > max_deg:
        #         max_deg = len(list(v.out_edges()))
        pattern_density[dict_patterns[idx_mp]] = (len(list(idx_mp.edges())) / 10 )

    pattern_density_sort = sorted(pattern_density.items(), key=operator.itemgetter(1))
    mpatterns = []  # sorted by edge density patterns
    for p, d in pattern_density_sort:
        mpatterns.append(p)

    intervals = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    X = [[] for _ in range(9)]
    X_title = []
    for idx_int in range(len(intervals)):
        for idx_m in range(len(mpatterns)):
            m = mpatterns[idx_m]
            X[idx_int].append(avg_f1[m][X_string[4]][intervals[idx_int]])
            X_title.append(m)

    plot_line(list(range(21)), X, '', X_title)



