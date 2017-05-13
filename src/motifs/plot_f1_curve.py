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


def plot_line(x, y1, y2, y3, y4, y5, l1, l2, l3, l4, l5, title, m, subim):
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    ax.imshow(subim, aspect='auto', extent=(5.8, 8.3, 0.8, 1), zorder=-1)
    # ax.add_patch(
    #     patches.Rectangle(
    #         (5.2, 0.77),
    #         2.5,
    #         0.28,
    #         fill=False,  # remove background
    #         linewidth=3
    #     )
    # )
    plt.xticks(x, size=30)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)
    ax.plot(x, y1,  marker='o', linestyle='-', label=l1, linewidth=3)
    ax.plot(x, y2, marker='o', linestyle='-', label=l2, linewidth=3)
    ax.plot(x, y3, marker='o', linestyle='-', label=l3, linewidth=3)
    ax.plot(x, y4, marker='o', linestyle='-', label=l4, linewidth=3)
    ax.plot(x, y5, marker='o', linestyle='-', label=l5, linewidth=3)

    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('F1 score', size=30)
    plt.xlabel('Number of intervals before inhibition (k)', size=30)

    ax.legend(loc='upper right', fancybox=True, shadow=True, fontsize=15)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    plt.ylim([0.05, 1.05])
    plt.xlim([1.5, 10.5])
    plt.savefig('v1/min_4/f1_plots_RF/Motif_' + m + '.png')
    plt.savefig('v1/min_4/f1_plots_RF/Motif_' + m + '.pdf')
    plt.close()
    # plt.show()


if __name__ == "__main__":
    X_string = ['Motif counts', 'Motif transitions', 'Motif weights', 'Temporal motifs', 'All']
    # X_string = ['Motif counts', 'Motif transitions', 'Temporal motifs', 'All']

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']
    motif_imgs_dir = "motif_patterns/"

    intervals = [2, 3, 4, 5, 6, 7, 8, 9, 10]

    avg_precision = pickle.load(open('../05/v1/avg_precision_RF.pickle', 'rb'))
    avg_recall = pickle.load(open('../05/v1/avg_recall_RF.pickle', 'rb'))
    avg_f1 = pickle.load(open('../05/v1/avg_f1_RF.pickle', 'rb'))
    # avg_random = pickle.load(open('../05/v1/avg_random.pickle', 'rb'))

    # Load the motif images
    imgs_dict = {}
    for m in motif_pattern:
        imgs_dict[m] = plt.imread(motif_imgs_dir + str(m) + ".png", format="png")

    for m in motif_pattern:
        Y = [[] for _ in range(5)]
        for idx_feat in range(len(X_string)):
            feat = X_string[idx_feat]
            for i in range(len(intervals)):
                interval = intervals[i]
                Y[idx_feat].append(avg_f1[m][X_string[idx_feat]][interval])
        # for i in range(len(intervals)):
        #     interval = intervals[i]
        #     Y[4].append(avg_f1[m][X_string[idx_feat]][interval])

        title = 'Motif: ' + m
        plot_line(intervals, Y[0], Y[1], Y[2], Y[3], Y[4], X_string[0], X_string[1],
                  X_string[2], X_string[3], X_string[4], title, m, imgs_dict[m])
        # exit()


