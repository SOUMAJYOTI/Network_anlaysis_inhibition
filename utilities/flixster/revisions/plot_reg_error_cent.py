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


def plot_line(x, y, l, x_labels):
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)


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
    plt.xticks(arange(len(x)), x_labels, size=30)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)

    markersList = ['o', 'v', '*', 'p', 's', '>', '^', 'h']
    colorsList = ['black', 'red', 'green', '#00BFFF', 'pink', '#D2691E', '#DAA520', '#BC8F8F']
    for idx in range(len(y)):  # lat one for random case
        ax.plot(x, y[idx], marker=markersList[idx], markersize=13, linestyle='--', linewidth=4,
                color=colorsList[idx], label=l[idx])


    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    # plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('Mean Absolute Error', size=40)
    plt.xlabel('Interval ranges before $N_{inhib}$', size=40)

    ax.legend(loc='lower right', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    # plt.ylim([0, 1])
    plt.xlim([-1, 5])
    # plt.savefig('outputs/v4/comp_reg/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/comp_reg/Motif_' + m + '.pdf')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.pdf')
    plt.show()
    plt.close()


if __name__ == "__main__":

    X_string = ['Pagerank', 'Degree Entropy', 'Clustering', 'Nodal Degree', 'Betweenness', 'All']

    motif_imgs_dir = "motif_patterns/"

    intervals = [0, 2,  4, 6, 8]

    avg_reg_edges = pickle.load(open('prediction/avg_reg_error_cent_edges.pickle', 'rb'))

    # avg_random = pickle.load(open('../05/v1/avg_random.pickle', 'rb'))

    # Load the motif images
    intervals_tit = []
    for idx in range(5):
        intervals_tit.append(str('[') + str(intervals[idx]+1) + str(',') + str(intervals[idx] + 3) + str(']'))

        # intervals_tit = np.array(intervals) + 1
    # intervals_tit.append('')

    labels = []
    Y_edges = [[] for _ in range(6)]
    for idx_feat in range(len(X_string)):
        if X_string[idx_feat] == 'Clustering':
            labels.append('Pagerank')
        elif X_string[idx_feat] == 'Pagerank':
            labels.append('Clustering')
        else:
            labels.append(X_string[idx_feat])
        feat = X_string[idx_feat]
        for i in range(len(intervals)):
            interval = intervals[i]
            Y_edges[idx_feat].append(avg_reg_edges[X_string[idx_feat]][interval])

    # intervals_tit.append()
    plot_line(range(len(Y_edges[0])), Y_edges, labels, intervals_tit)
    # exit()


