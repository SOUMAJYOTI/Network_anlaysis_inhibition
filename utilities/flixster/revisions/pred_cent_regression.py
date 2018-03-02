from pylab import *
import pickle
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


def motif_X_Y(X_feat, Y_label, interval_start, interval_end, mlist):
    # print('Preparing the data.....')
    X = []
    Y = []

    X_mid = {}
    Y_mid = {}
    # stack the features for the intervals to get multiple features
    for idx_f in range(len(X_feat)):
        X_cent = X_feat[idx_f]

        for m in range(len(mlist)):
            mid = mlist[m]
            if mid not in X_cent:
                temp_X = []
                # continue
                # temp_X = [0] * 2 #(intervals)
            else:
                temp_X = X_cent[mid][interval_start:interval_end]

                temp_X_n = []
                for x in temp_X:
                    if x == {}:
                        continue
                    else:
                        temp_X_n.append(x)

                nanmean = np.nanmean(np.array(temp_X_n))

                temp_X_n = [nanmean if math.isnan(x) else x for x in temp_X_n]
            # consider intervals as a parameter
            # if len(temp_X) > intervals:
            #     temp_X = X_motif[mid][:intervals]

            # Zero pad those which do not have k intervals
            if len(temp_X_n) == 0:
                temp_X_n += [0] * (2 - len(temp_X_n))
            else:
                temp_X_n += [np.mean(np.array(temp_X_n))] * (2 - len(temp_X_n))
            # print(temp_X)
            if mid not in X_mid:
                X_mid[mid] = []
            X_mid[mid].extend(temp_X_n)


    cnt_mids = 0
    # print(len(Y_motif), motif_label)
    for m in range(len(mlist)):
        mid = mlist[m]

        Y_mid[mid] = 0.
        if mid not in Y_label:
            continue

        Y_mid[mid] = Y_label[mid]
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
    # Load the interim computed features
    X_pr = pickle.load(open('centralities/feat_intervals/pr.pickle', 'rb'))
    X_en = pickle.load(open('centralities/feat_intervals/en.pickle', 'rb'))
    X_bw = pickle.load(open('centralities/feat_intervals/bw.pickle', 'rb'))
    X_cl = pickle.load(open('centralities/feat_intervals/cl.pickle', 'rb'))
    X_deg = pickle.load(open('centralities/feat_intervals/deg.pickle', 'rb'))

    # Y_cover = pickle.load(open('interim_feat/Y_num_edges.pickle', 'rb'))
    Y_cover = pickle.load(open('actual_num_edges.pickle', 'rb'))

    Y_edges = {}
    for mid in Y_cover:
        for g in Y_cover[mid]:
            Y_edges[mid] = Y_cover[mid][g]
            break

    # aggregate all the mids
    mid_1 = set(list(X_pr.keys()))
    mid_2 = set(list(X_cl.keys()))
    mid_3 = set(list(X_en.keys()))
    mid_4 = set(list(X_deg.keys()))
    mid_5 = set(list(X_bw.keys()))


    mlist = set([])
    mlist = mlist.union(mid_1)
    mlist = mlist.union(mid_2)
    mlist = mlist.union(mid_3)
    mlist = mlist.union(mid_4)
    mlist = mlist.union(mid_5)
    mlist = list(mlist)

    X_string = ['Pagerank', 'Degree Entropy', 'Clustering', 'Nodal Degree', 'Betweenness', 'All']
    X_total = [[X_pr], [X_en], [X_cl], [X_deg], [X_bw], [X_pr, X_en, X_cl, X_deg, X_bw]] # start with one feature
    # X_total = [[X_weights]]  # start with one feature
    avg_random = {}

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
            print("Interval: ", interval_start, " ,Feature: ", X_string[idx_feat], )

            if X_string[idx_feat] not in avg_error:
                avg_error[X_string[idx_feat]] = {}

            avg_er = 0.

            X, Y = motif_X_Y(X_all, Y_edges, interval_start, interval_end,  mlist)
            X = np.array(X)
            Y = np.array(Y)

            print(Y.shape, X.shape)
            # if len(Y) < 600:
            #     avg_error[mp][X_string[idx_feat]][interval_start] = avg_er / 10.
            #     continue

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
                # Y_random = []
                # for idx_r in range(X_test.shape[0]):
                #     Y_random.append(random.sample([0, 1], 1))
                #
                # Y_random = np.array(Y_random)

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

            avg_error[X_string[idx_feat]][interval_start] = avg_er / 10.

    pickle.dump(avg_error, open('prediction/avg_reg_error_edges.pickle', 'wb'))