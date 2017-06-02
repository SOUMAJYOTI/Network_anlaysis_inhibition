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
from sklearn import ensemble
import random


def getFolds(Y):
    train_folds = []
    test_folds = []
    kf = KFold(n_splits=10)
    for train_index, test_index in kf.split(Y):
        train_folds.append(train_index)
        test_folds.append(test_index)

    return train_folds, test_folds


def plot_line(x, y1, y2, y3, y4):
    fig = plt.figure()
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    plt.xticks(x, size=30)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)
    ax.plot(x, y1,  marker='o', linestyle='-', linewidth=3)
    ax.plot(x, y2, marker='o', linestyle='X', linewidth=3)
    ax.plot(x, y3, marker='o', linestyle='-', linewidth=3)
    ax.plot(x, y4, marker='o', linestyle='-', linewidth=3)

    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    # plt.title('Network_cover vs Motif edge density ', color='#000000', weight="bold", size=50)
    plt.ylabel('% of edges covered by motif', size=40)
    plt.xlabel('Motif edge density', size=50)

    # ax.legend(loc='upper right', fancybox=True, shadow=True, fontsize=15)
    plt.grid(True)

    plt.show()


def motif_X_Y(X_feat, Y_label, intervals, motif_label, mlist):
    # print('Preparing the data.....')
    X = []
    Y = []

    X_mid = {}
    Y_mid = {}
    # stack the features for the intervals to get multiple features
    for idx_f in range(len(X_feat)):
        if motif_label not in X_feat[idx_f]:
            continue
        X_motif = X_feat[idx_f][motif_label] # --> for each motif

        for m in range(len(mlist)):
            mid = mlist[m]
            if mid not in X_motif:
                temp_X = [0] * (intervals)
            else:
                temp_X = X_motif[mid]

            # consider intervals as a parameter
            if len(temp_X) > intervals:
                temp_X = X_motif[mid][:intervals]
            # Zero pad those which do not have k intervals
            temp_X += [0] * (intervals - len(temp_X))
            # print(temp_X)
            if mid not in X_mid:
                X_mid[mid] = []
            X_mid[mid].extend(temp_X)

    Y_motif = Y_label[motif_label]
    # find the mean of the coverage for this motif
    mean_cover = np.mean(np.array(list(Y_motif.values())))
    for m in range(len(mlist)):
        mid = mlist[m]
        Y_mid[mid] = 0.

        if mid not in Y_motif:
            continue
        if Y_motif[mid] > mean_cover:
            Y_mid[mid] = 1.

    for mid in X_mid:
        X.append(X_mid[mid])
        Y.append(Y_mid[mid])

    return X, Y


if __name__ == "__main__":
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    # Load the interim computed features
    X_counts = pickle.load(open('interim_feat/X_count.pickle', 'rb'))
    X_weights = pickle.load(open('interim_feat/X_weights.pickle', 'rb'))
    X_trans = pickle.load(open('interim_feat/X_trans.pickle', 'rb'))
    X_temporal = pickle.load(open('interim_feat/X_temporal.pickle', 'rb'))

    Y_cover = pickle.load(open('interim_feat/Y_cover_min_4.pickle', 'rb'))
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))

    # aggregate all the mids
    mid_1 = set(list(X_counts['M10'].keys()))
    mid_2 = set(list(X_weights['M10'].keys()))
    mid_3 = set(list(X_trans['M10'].keys()))
    mid_4 = set(list(X_temporal['M10'].keys()))

    mlist = set([])
    mlist = mlist.union(mid_1)
    mlist = mlist.union(mid_2)
    mlist = mlist.union(mid_3)
    mlist = mlist.union(mid_4)
    mlist = list(mlist)

    X_string = ['Motif counts', 'Motif transitions', 'Motif weights', 'Temporal motifs', 'All']
    X_total = [[X_counts], [X_trans], [X_weights], [X_temporal], [X_counts, X_trans, X_temporal]] # start with one feature
    # X_total = [[X_weights]]  # start with one feature
    avg_precision = {}
    avg_recall = {}
    avg_f1 = {}
    avg_random = {}

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']

    intervals = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    avg_precision = {}
    avg_recall = {}
    avg_f1 = {}
    avg_random = {}
    for idx_int in range(len(intervals)):
        interval = intervals[idx_int]
        for idx_feat in range(len(X_total)):
            X_all = X_total[idx_feat]
            for mp in motif_pattern:
                print("Interval: ", interval, " ,Feature: ", X_string[idx_feat], " , Motif pattern: ", mp)
                if mp not in avg_precision:
                    avg_precision[mp] = {}
                    avg_recall[mp] = {}
                    avg_f1[mp] = {}
                    avg_random[mp] = {}

                if X_string[idx_feat] not in avg_precision[mp]:
                    avg_precision[mp][X_string[idx_feat]] = {}
                    avg_recall[mp][X_string[idx_feat]] = {}
                    avg_f1[mp][X_string[idx_feat]] = {}
                    avg_random[mp][X_string[idx_feat]] = {}

                avg_p = 0.
                avg_r = 0.
                avg_fs = 0.
                avg_rn = 0

                X, Y = motif_X_Y(X_all, Y_cover, interval, mp, mlist)
                X = np.array(X)
                Y = np.array(Y)

                if len(Y) < 10:
                    avg_precision[mp][X_string[idx_feat]][interval] = avg_p / 10.
                    avg_recall[mp][X_string[idx_feat]][interval] = avg_r / 10.
                    avg_f1[mp][X_string[idx_feat]][interval] = avg_fs / 10.
                    avg_random[mp][X_string[idx_feat]][interval] = avg_rn / 10.
                    continue

                train_fold, test_fold = getFolds(Y)

                """ Logistic Regression """
                clf = linear_model.LogisticRegression()
                clf = ensemble.RandomForestClassifier()

                for idx_fold in range(len(train_fold)):
                    X_train = X[train_fold[idx_fold]]
                    Y_train = Y[train_fold[idx_fold]]

                    X_test = X[test_fold[idx_fold]]
                    Y_test = Y[test_fold[idx_fold]]

                    Y_random = []
                    for idx_r in range(X_test.shape[0]):
                        Y_random.append(random.sample([-1, 1], 1))

                    Y_random = np.array(Y_random)

                    clf.fit(X_train, Y_train)
                    Y_predict = clf.predict(X_test)
                    # print(Y_predict)

                    avg_p += sklearn.metrics.precision_score(Y_test, Y_predict)
                    avg_r += sklearn.metrics.recall_score(Y_test, Y_predict)
                    avg_fs += sklearn.metrics.f1_score(Y_test, Y_predict)
                    avg_rn += sklearn.metrics.f1_score(Y_test, Y_random)

                # print(avg_fs/10)
                avg_precision[mp][X_string[idx_feat]][interval] = avg_p/10.
                avg_recall[mp][X_string[idx_feat]][interval] = avg_r/10.
                avg_f1[mp][X_string[idx_feat]][interval] = avg_fs/10.
                avg_random[mp][X_string[idx_feat]][interval] = avg_rn/10.

    pickle.dump(avg_precision, open('../05/v1/avg_precision_RF.pickle', 'wb'))
    pickle.dump(avg_recall, open('../05/v1/avg_recall_RF.pickle', 'wb'))
    pickle.dump(avg_f1, open('../05/v1/avg_f1_RF.pickle', 'wb'))
    # pickle.dump(avg_random, open('../05/avg_random.pickle', 'wb'))