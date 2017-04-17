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
        if cnt_mids > 50:
            break

    return X_feat, mid_list


def motif_weights_feat(motif_data, motifs_patterns, k):
    print('Preparing motif weights features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0

    for mid in motif_data:
        for int in range(1, k):
            for m in motif_data[mid][int]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                X_feat[motifs_patterns[pat]][mid].append(motif_data[mid][int][m])

        cnt_mids += 1
        if cnt_mids > 50:
            break
    return X_feat


def motif_trans_feat(motif_data, motifs_patterns, k):
    print('Preparing motif transition features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0

    for mid in motif_data:
        for int in range(1, k):
            for m in motif_data[mid][int]:
                pat = checkIsomorphism(mpatterns, m)


                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                # Stack all the motifs transition for next interval in one list
                # as a feature
                temp = []
                for m1 in motif_data[mid][int][m]:
                    temp.append(motif_data[mid][int][m][m1])

                X_feat[motifs_patterns[pat]][mid].append(temp)
        cnt_mids += 1
        if cnt_mids > 50:
            break

    return X_feat


def motif_temporal_feat(motif_data, motifs_patterns, k):
    print('Preparing temporal motif features....')
    X_feat = {}
    mpatterns = list(motifs_patterns.keys())
    cnt_mids = 0

    for mid in motif_data:
        for int in range(1, k):
            for m in motif_data[mid][int]:
                pat = checkIsomorphism(mpatterns, m)

                if motifs_patterns[pat] not in X_feat:
                    X_feat[motifs_patterns[pat]] = {}

                if mid not in X_feat[motifs_patterns[pat]]:
                    X_feat[motifs_patterns[pat]][mid] = []

                X_feat[motifs_patterns[pat]][mid].append(motif_data[mid][int][m])
        cnt_mids += 1
        if cnt_mids > 50:
            break

    return X_feat


def motif_coverage(motif_data, motifs_patterns):
    print('Preparing labels......')
    mpatterns = list(motifs_patterns.keys())

    Y_cover = {}
    for mid in motif_data:
        for m in motif_data[mid]:
            pat = checkIsomorphism(mpatterns, m)
            if motifs_patterns[pat] not in Y_cover:
                Y_cover[motifs_patterns[pat]] = {}

            Y_cover[motifs_patterns[pat]][mid] = motif_data[mid][m]

    return Y_cover


def getFolds(Y):
    train_folds = []
    test_folds = []
    kf = KFold(n_splits=10)
    for train_index, test_index in kf.split(Y):
        train_folds.append(train_index)
        test_folds.append(test_index)

    return train_folds, test_folds


def motif_X_Y(X_feat, Y_label, motif_label, mlist):
    print('Preparing the data.....')
    X = []
    Y = []
    for idx_f in range(len(X_feat)):
        X_motif = X_feat[idx_f][motif_label]
        Y_motif = Y_label[motif_label]
        # print(Y_motif)
        for mid in range(len(mlist)):

            if mid not in X_motif:
                X_motif[mid] += [0] * (5 - len(X_motif[mid]))

            # Zero pad those which do not have k intervals
            X_motif[mid] += [0] * (5 - len(X_motif[mid]))
            X.append(X_motif[mid])

            # This should be executed only once - not for all features repeatedly
            if idx_f == 0:
                if mid not in Y_motif:
                    Y.append(0.)
                else:
                    if Y_motif[mid] > 50:
                        Y.append(1.)
                    else:
                        Y.append(0.)


    return X, Y


if __name__ == "__main__":
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    motif_counts = pickle.load(open('output_data/motifs_inhib_count.pickle', 'rb'))
    motif_weights = pickle.load(open('output_data/motifs_inhib_weights.pickle', 'rb'))
    motif_transition = pickle.load(open('output_data/motif_inhib_transition.pickle', 'rb'))
    motif_temporal = pickle.load(open('output_data/frontiers_gt_inhib_count.pickle', 'rb'))
    motif_covers = pickle.load(open('output_data/motifs_coverage_min_4.pickle', 'rb'))

    """ Number of intervals to consider"""
    int = 5
    X_count, mlist = motif_counts_feat(motif_counts, dict_patterns, int)
    X_weights = motif_weights_feat(motif_weights, dict_patterns, int)
    X_trans = motif_trans_feat(motif_transition, dict_patterns, int)
    X_temporal = motif_temporal_feat(motif_temporal, dict_patterns, int)

    # motif_weights_feat(motif_weights, dict_patterns, int)
    Y_cover = motif_coverage(motif_covers, dict_patterns)

    exit()

    X_all = [X_count]
    avg_precision = {}
    avg_recall = {}
    avg_f1 = {}

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20', 'M21']

    for mp in motif_pattern:
        avg_precision[mp] = 0.
        avg_recall[mp] = 0.
        avg_f1[mp]= 0.
        X, Y = motif_X_Y(X_all, Y_cover, mp, list)
        X = np.array(X)
        Y = np.array(Y)

        train_fold, test_fold = getFolds(Y)

        """ Logistic Regression """
        clf = linear_model.LogisticRegression()

        for idx_fold in range(len(train_fold)):
            X_train = X[train_fold[idx_fold]]
            Y_train = Y[train_fold[idx_fold]]

            X_test = X[test_fold[idx_fold]]
            Y_test = Y[test_fold[idx_fold]]

            clf.fit(X_train, Y_train)
            Y_predict = clf.predict(X_test)

            avg_precision[mp] += sklearn.metrics.precision_score(Y_test, Y_predict)
            avg_recall[mp] += sklearn.metrics.recall_score(Y_test, Y_predict)
            avg_f1[mp] += sklearn.metrics.f1_score(Y_test, Y_predict)

        avg_precision[mp] /= 10.
        avg_recall[mp] /= 10.
        avg_f1[mp] /= 10.
