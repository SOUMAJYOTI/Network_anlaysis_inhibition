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

global_pat = {}


def checkIsomorphism(graph_list, g):
    for gr in range(len(graph_list)):
        if gtt.isomorphism(graph_list[gr], g):
            return graph_list[gr]
    return False


def cent_feat(cent_data, k):
    print('Preparing motif counts features....')
    X_feat = {}
    cnt_mids = 0
    mid_list = []

    # Data is stored in [mid][interval][motif] format
    # output is in the format [motif_pattern ID][mid] -->  k intervals of data stacked
    for mid in cent_data:
        mid_list.append(mid)
        for interval in range(1, k):
            if mid not in X_feat:
                X_feat[mid] = []

            X_feat[mid].append(cent_data[mid][interval])
        cnt_mids += 1
        # if cnt_mids > 50:
        #     break

    return X_feat, mid_list


if __name__ == "__main__":
    pr_feat = pickle.load(open('centralities/en.pickle', 'rb'))
    # """ Number of intervals to consider"""
    interval = 20

    # motif_weights_feat(motif_weights, dict_patterns, int)
    X_feat, midList = cent_feat(pr_feat, interval)
    # for mid in midList:
    #     print(X_feat[mid])

    pickle.dump(X_feat, open('centralities/feat_intervals/en.pickle', 'wb'))



