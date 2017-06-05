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


def plot_line(x, y1, y2, y3, y4, y5, y6, y7, l1, l2, l3, l4, l5, l6, l7, title, m, subim):
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    max_val = np.max(y1)

    # ax.imshow(subim, aspect='auto', extent=(4.8, 7.3, max_val+5.8, max_val + 9), zorder=-1)
    # ax.add_patch(
    #     patches.Rectangle(
    #         (5.2, 0.77),
    #         2.5,
    #         0.28,
    #         fill=False,  # remove background
    #         linewidth=3
    #     )
    # )
    plt.xticks(arange(len(x)), x, size=30)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)

    x = np.array(range(len(x)-2)) + 1
    # print(x.shape, len(y1))
    ax.plot(x, y1,  marker='o', linestyle='-', label=l1, linewidth=3)
    ax.plot(x, y2, marker='o', linestyle='-', label=l2, linewidth=3)
    ax.plot(x, y3, marker='o', linestyle='-', label=l3, linewidth=3)
    ax.plot(x, y4, marker='o', linestyle='-', label=l4, linewidth=3)
    ax.plot(x, y5, marker='o', linestyle='-', label=l5, linewidth=3)
    ax.plot(x, y6, marker='o', linestyle='-', label=l6, linewidth=3)
    ax.plot(x, y7, marker='o', linestyle='-', label=l7, linewidth=3)
    # ax.plot(x, y8, marker='o', linestyle='-', label=l8, linewidth=3)


    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    # plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('F1', size=40)
    plt.xlabel('Interval ranges before $N_{inhib}$', size=40)

    ax.legend(loc='upper right', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    # plt.ylim([0, 1])
    plt.xlim([0, 6])
    plt.savefig('outputs/v5/f1_comp/Motif_' + m + '.png')
    plt.savefig('outputs/v5/f1_comp/Motif_' + m + '.pdf')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.pdf')
    plt.close()
    # plt.show()


if __name__ == "__main__":
    X_string = ['$NC_{thresh}$=10', '$NC_{thresh}$=35', '$NC_{thresh}$=40', '$NC_{thresh}$=45', '$NC_{thresh}$=50', '$NC_{thresh}$=55',
                '$NC_{thresh}$=60',]

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']
    motif_imgs_dir = "motif_patterns/"

    intervals = [0, 2, 4, 6, 8]

    # avg_precision = pickle.load(open('../05/v2/avg_precision_RF_motifs.pickle', 'rb'))
    # avg_recall = pickle.load(open('../05/v2/avg_recall_RF_motifs.pickle', 'rb'))
    avg_f1_30 = pickle.load(open('../05/outputs/v5/avg_f1_RF_10.pickle', 'rb'))
    avg_f1_35 = pickle.load(open('../05/outputs/v5/avg_f1_RF_35.pickle', 'rb'))
    avg_f1_40 = pickle.load(open('../05/outputs/v5/avg_f1_RF_40.pickle', 'rb'))
    avg_f1_45 = pickle.load(open('../05/outputs/v5/avg_f1_RF_45.pickle', 'rb'))
    avg_f1_50 = pickle.load(open('../05/outputs/v5/avg_f1_RF_50.pickle', 'rb'))
    avg_f1_55 = pickle.load(open('../05/outputs/v5/avg_f1_RF_55.pickle', 'rb'))
    avg_f1_60 = pickle.load(open('../05/outputs/v5/avg_f1_RF_60.pickle', 'rb'))
    # avg_f1_65 = pickle.load(open('../05/outputs/v5/avg_f1_RF_65.pickle', 'rb'))

    avg_f1 = [avg_f1_30, avg_f1_35, avg_f1_40, avg_f1_45, avg_f1_50, avg_f1_55, avg_f1_60, ]
    # avg_random = pickle.load(open('../05/v1/avg_random.pickle', 'rb'))

    # Load the motif images
    intervals_tit = ['']
    for idx in range(5):
        intervals_tit.append(str('[') + str(intervals[idx]+1) + str(',') + str(intervals[idx] + 3) + str(']'))
        # intervals_tit = np.array(intervals) + 1
    intervals_tit.append('')
    imgs_dict = {}
    for m in motif_pattern:
        imgs_dict[m] = plt.imread(motif_imgs_dir + str(m) + ".png", format="png")

    for m in motif_pattern:
        Y = [[] for _ in range(8)]
        for idx_f1 in range(len(avg_f1)):
            # feat = X_string['All']
            for i in range(len(intervals)):
                interval = intervals[i]
                Y[idx_f1].append(avg_f1[idx_f1][m]['All'][interval])
        # for i in range(len(intervals)):
        #     interval = intervals[i]
        #     Y[4].append(avg_f1[m][X_string[idx_feat]][interval])

        # intervals_tit.append()
        title = 'Motif: ' + m
        plot_line(intervals_tit, Y[0], Y[1], Y[2], Y[3], Y[4], Y[5], Y[6],  X_string[0], X_string[1],
                  X_string[2], X_string[3], X_string[4], X_string[5], X_string[6], title, m, imgs_dict[m])
        # exit()


