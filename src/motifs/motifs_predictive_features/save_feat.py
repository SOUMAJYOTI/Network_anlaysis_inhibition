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


def motif_counts_feat(motif_data, motifs_patterns, k):
    print('Preparing motif counts features....')
    X_feat = {}
    cnt_mids = 0
    mpatterns = list(motifs_patterns.keys())
    mid_list = []

    # Data is stored in [mid][interval][motif] format
    # output is in the format [motif_pattern ID][mid] -->  k intervals of data stacked
    for mid in motif_data:
        mid_list.append(mid)
        for interval in range(1, k):
            for m in motif_data[mid][interval]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                X_feat[motifs_patterns[pat]][mid].append(motif_data[mid][interval][m])
        cnt_mids += 1
        # if cnt_mids > 50:
        #     break

    return X_feat, mid_list


def motif_weights_feat(motif_data, motifs_patterns, k):
    print('Preparing motif weights features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0
    mid_list = []

    # Data is stored in [mid][interval][motif] format
    # output is in the format [motif_pattern ID][mid] -->  k intervals of data stacked
    for mid in motif_data:
        mid_list.append(mid)
        for int in range(1, k):
            for m in motif_data[mid][int]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                X_feat[motifs_patterns[pat]][mid].append(motif_data[mid][int][m])

        cnt_mids += 1
        # if cnt_mids > 5:
        #     break
    return X_feat, mid_list


def motif_trans_feat(motif_data, motifs_patterns, k):
    print('Preparing motif transition features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0
    mid_list = []

    motif_data_inv = {} # store data in format [mid][interval][m_5][m_4] --> from [m_4][m_5]
    for mid in motif_data:
        motif_data_inv[mid] = {}
        for interval in range(1, k):
            motif_data_inv[mid][interval] = {}
            for m_4 in motif_data[mid][interval]:
                for m_5 in motif_data[mid][interval][m_4]:
                    if m_5 not in motif_data_inv[mid][interval]:
                        motif_data_inv[mid][interval][m_5] = {}

                    motif_data_inv[mid][interval][m_5][m_4] = motif_data[mid][interval][m_4][m_5] # reverse the indices

    # Data is stored in [mid][interval][motif_inp][motif_transitioned] format
    # output is in the format [motif_pattern ID][mid] -->  k intervals of data - mean of transitions
    for mid in motif_data_inv:
        mid_list.append(mid)
        for int in range(1, k):
            for m in motif_data_inv[mid][int]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                # Stack all the motifs transition for next interval in one list
                # as a feature - can take the mean of the transitions
                temp = []
                for m1 in motif_data_inv[mid][int][m]:
                    temp.append(motif_data_inv[mid][int][m][m1])

                X_feat[motifs_patterns[pat]][mid].append(np.mean(np.array(temp)))
        cnt_mids += 1
        # if cnt_mids > 5:
        #     break

    return X_feat, mid_list


def motif_temporal_feat(motif_data, motifs_patterns, k):
    print('Preparing temporal motif features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0
    mid_list = []

    # Data is stored in [mid][interval][motif_inp] format --> data is for lambda-lt motifs
    # output is in the format [motif_pattern ID][mid] -->  k intervals of data
    for mid in motif_data:
        mid_list.append(mid)
        for int in range(1, k):
            for m in motif_data[mid][int]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                X_feat[motifs_patterns[pat]][mid].append(motif_data[mid][int][m])
        cnt_mids += 1
        # if cnt_mids > 5:
        #     break

    return X_feat, mid_list


def motif_coverage(motif_data, motifs_patterns):
    print('Preparing labels......')
    mpatterns = list(motifs_patterns.keys())

    # Data is stored in [mid][motif_pattern]  ---> only for the inhibition interval
    # output is in the format [motif_pattern ID][mid] --->  returns the coverage percentage
    Y_cover = {}
    for mid in motif_data:
        for m in motif_data[mid]:
            pat = checkIsomorphism(mpatterns, m)
            if motifs_patterns[pat] not in Y_cover:
                Y_cover[motifs_patterns[pat]] = {}

            Y_cover[motifs_patterns[pat]][mid] = motif_data[mid][m]

    return Y_cover


if __name__ == "__main__":
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    motif_counts = pickle.load(open('output_data/motifs_inhib_count.pickle', 'rb'))
    motif_weights = pickle.load(open('output_data/motifs_inhib_weights.pickle', 'rb'))
    motif_transition = pickle.load(open('output_data/motif_inhib_transition.pickle', 'rb'))
    motif_temporal = pickle.load(open('output_data/frontiers_inhib_count.pickle', 'rb'))
    motif_covers = pickle.load(open('output_data/motifs_coverage_min_4.pickle', 'rb'))

    # # Save the figures of the motif patterns
    # cnt_motif = 0
    # for g in dict_patterns:
    #     pos = gtd.arf_layout(g)
    #     gtd.graph_draw(g, pos=pos, output="motif_patterns/pdf/" + str(dict_patterns[g]) + ".pdf")
    #     cnt_motif += 1
    #
    # # print(dict_patterns)
    # exit()
    #
    # """ Number of intervals to consider"""
    # interval = 10
    # X_count, mlist_count = motif_counts_feat(motif_counts, dict_patterns, interval)
    # pickle.dump(X_count, open('interim_feat/X_count.pickle', 'wb'))
    # X_weights, mlist_weights = motif_weights_feat(motif_weights, dict_patterns, interval)
    # pickle.dump(X_weights, open('interim_feat/X_weights.pickle', 'wb'))
    # X_trans, mlist_trans = motif_trans_feat(motif_transition, dict_patterns, interval)
    # pickle.dump(X_trans, open('interim_feat/X_trans.pickle', 'wb'))
    # X_temporal, mlist_temp = motif_temporal_feat(motif_temporal, dict_patterns, interval)
    # pickle.dump(X_temporal, open('interim_feat/X_temporal.pickle', 'wb'))

    # motif_weights_feat(motif_weights, dict_patterns, int)
    Y_cover = motif_coverage(motif_covers, dict_patterns)
    pickle.dump(Y_cover, open('interim_feat/Y_cover_min_4.pickle', 'wb'))


