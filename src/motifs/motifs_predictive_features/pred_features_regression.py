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
# Cross validation for the motif features


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


def motif_X_Y(X_feat, Y_label, interval_start, interval_end, motif_label, mlist):
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
                temp_X = []
                # continue
                # temp_X = [0] * 2 #(intervals)
            else:
                temp_X = X_motif[mid][interval_start:interval_end]

            # consider intervals as a parameter
            # if len(temp_X) > intervals:
            #     temp_X = X_motif[mid][:intervals]

            # Zero pad those which do not have k intervals
            if len(temp_X) == 0:
                temp_X += [0] * (2 - len(temp_X))
            else:
                temp_X += [np.mean(np.array(temp_X))] * (2 - len(temp_X))
            # print(temp_X)
            if mid not in X_mid:
                X_mid[mid] = []
            X_mid[mid].extend(temp_X)

    Y_motif = Y_label[motif_label]

    cnt_mids = 0
    # print(len(Y_motif), motif_label)
    for m in range(len(mlist)):
        mid = mlist[m]

        Y_mid[mid] = 0.
        if mid not in Y_motif:
            continue

        Y_mid[mid] = Y_motif[mid]
        # print(Y_motif[mid])
        cnt_mids += 1

    # print(cnt_mids)

    for mid in Y_mid:
        if mid not in X_mid:
            continue
        if Y_mid[mid] == 0.:
            # continue
            Y_mid[mid] = np.mean(np.array(list(Y_mid.values())))
        X.append(X_mid[mid])
        Y.append(Y_mid[mid])

    return X, Y


if __name__ == "__main__":
    dict_patterns = pickle.load(open('dict_patterns.pickle', 'rb'))
    # Load the interim computed features
    X_counts = pickle.load(open('interim_feat_v1/X_count.pickle', 'rb'))
    X_weights = pickle.load(open('interim_feat_v1/X_weights.pickle', 'rb'))
    X_trans = pickle.load(open('interim_feat_v1//X_trans.pickle', 'rb'))
    X_temporal = pickle.load(open('interim_feat_v1/X_temporal.pickle', 'rb'))

    # Y_cover = pickle.load(open('interim_feat/Y_num_edges.pickle', 'rb'))
    Y_cover = pickle.load(open('interim_feat/Y_cover_edges.pickle', 'rb'))

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

    X_string = ['Motif counts', 'Motif transitions', 'Temporal motifs', 'All']
    X_total = [[X_counts], [X_trans], [X_temporal], [X_counts, X_trans, X_temporal]] # start with one feature
    # X_total = [[X_weights]]  # start with one feature
    avg_random = {}

    motif_pattern = ['M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10'
                     , 'M11', 'M12', 'M13', 'M14', 'M15', 'M16', 'M17', 'M18', 'M19', 'M20']

    intervals = [2, 4, 6, 8, 10]
    # avg_precision = {}
    # avg_recall = {}
    # avg_f1 = {}
    # avg_random = {}
    avg_error = {}
    for idx_int in range(len(intervals)):
        interval_start = intervals[idx_int] - 2
        interval_end = intervals[idx_int]
        for idx_feat in range(len(X_total)):
            X_all = X_total[idx_feat]
            for mp in motif_pattern:
                print("Interval: ", interval_start, " ,Feature: ", X_string[idx_feat], " , Motif pattern: ", mp)
                if mp not in avg_error:
                    avg_error[mp] = {}

                if X_string[idx_feat] not in avg_error[mp]:
                    avg_error[mp][X_string[idx_feat]] = {}

                avg_er = 0.

                X, Y = motif_X_Y(X_all, Y_cover, interval_start, interval_end, mp, mlist)
                X = np.array(X)
                Y = np.array(Y)

                print(Y.shape, X.shape)
                if len(Y) < 600:
                    avg_error[mp][X_string[idx_feat]][interval_start] = avg_er / 10.
                    continue

                train_fold, test_fold = getFolds(Y)

                """ Linear LASSO Regression """
                # clf = linear_model.LogisticRegression()
                # clf = ensemble.RandomForestClassifier()

                for idx_fold in range(len(train_fold)):
                    X_train = X[train_fold[idx_fold]]
                    Y_train = Y[train_fold[idx_fold]]

                    X_test = X[test_fold[idx_fold]]
                    Y_test = Y[test_fold[idx_fold]]

                    # print(Y_train)
                    # print(Y_test)
                    Y_random = []
                    for idx_r in range(X_test.shape[0]):
                        Y_random.append(random.sample([0, 1], 1))

                    Y_random = np.array(Y_random)

                    # print(np.array(X_train).shape)
                    # clf.fit(X_train, Y_train)
                    # Y_predict = clf.predict(X_test)
                    # print(Y_predict)

                    alphas = [0.01, 0.02, 0.03, 0.04]
                    regr = linear_model.Lasso()
                    scores = [regr.set_params(alpha=alpha).fit(X_train, Y_train).score(X_test, Y_test) for alpha in
                              alphas]
                    best_alpha = alphas[scores.index(max(scores))]
                    regr.alpha = best_alpha
                    regr.fit(X_train, Y_train)
                    # print(Y_train)
                    # print(regr.coef_)

                    # print(regr.predict(X_test))
                    # print(Y_test)
                    # exit()
                    # The mean square error
                    avg_er += (np.mean(np.absolute((regr.predict(X_test) - Y_test))))

                avg_error[mp][X_string[idx_feat]][interval_start] = avg_er / 10.

    pickle.dump(avg_error, open('../05/outputs/v5/avg_reg_error_coverage.pickle', 'wb'))