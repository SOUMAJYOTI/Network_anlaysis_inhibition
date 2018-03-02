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



def checkIsomorphism(graph_list, g):
    for gr in graph_list:
        if gtt.isomorphism(gr, g):
            return gr
    return False


def plot_line(x, y1, y2,  l1, l2, title, m, subim):
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


    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    # plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('Mean Absolute Error', size=40)
    plt.xlabel('Interval ranges before $N_{inhib}$', size=40)

    ax.legend(loc='upper right', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    # plt.ylim([0, 1])
    plt.xlim([0, 6])
    # plt.savefig('outputs/v4/comp_reg/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/comp_reg/Motif_' + m + '.pdf')

    plt.show()
    plt.close()
    # plt.show()


if __name__ == "__main__":
    X_string = [ 'All']

    # X_string = ['Motif counts', 'Motif transitions', 'Temporal motifs', 'All']

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']
    motif_imgs_dir = "motif_patterns/"


    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    for m in dict_patterns:
        print(m, dict_patterns[m])

