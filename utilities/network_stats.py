# Find the steep and the inhibition regions of the cascades

__author__ = 'ssarka18'

import pandas as pd
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
import math
from math import factorial
from scipy.misc import factorial
from scipy.optimize import minimize
import csv
import os
import statsmodels.api as sm
import pickle

path = 'rt_df.csv'


class FormCascadesOrder(object):
    def __init__(self):
        self.flag_cascade = []

    def cascades_form(self):
        flag_scanned = {}
        df = pd.read_csv(path)
        cnt = 0
        size_cascades = []
        window_size = []
        first_time = []
        steep_points_growth = []
        steep_points_size = []
        growth = []
        time_diff_points = []
        time_values = []
        cnt_300 = 0
        cnt_total = 0
        # df_temp = df[df['3'] == 8890888985501908]
        # print(df_temp)
        time_cascades = []
        size_cascades = []
        for index, row in df.iterrows():
            if str(row['3']) not in flag_scanned:
                flag_scanned[str(row['3'])] = 1
                x_time = []
                time_inst = []
                reaction_times = []
                node_rt_time = {}
                node_reaction_time = []
                map_size_time = {}
                cascade_set = df[(df['3'] == row['3'])]
                cnt += 1
                cascade_set_1 = list(cascade_set['0'])
                cascade_set_2 = list(cascade_set['1'])
                cascade_tgt = list(cascade_set['0'])
                cascade_src = list(cascade_set['1'])
                cascade_time_idx = []
                len_cascade = len(list(set(cascade_set_1 + cascade_set_2)))
                if len_cascade <= 300: # DO NOT CONSIDER CASCADES OF SIZE LESS THAN 300
                    continue
                cnt_300 += 1
                print(cnt_300, cnt)
                if cnt_300 < 2000:
                    # with open('8.txt', 'w') as f:
                    #     a = csv.writer(f, lineterminator= '\n')
                    #     all = []
                    #     for i in range(len(cascade_src)):
                    #         f.write(cascade_src[i])
                    #         f.write('\t')
                    #         f.write(cascade_tgt[i])
                    #         f.write(('\n'))
                    #     #a.writerows(all)
                    # f.close()
                    # break

                    #print("Cacsade ", cnt, "of final size: ", len(cascade_set))
                    cascade_rt_time = list(cascade_set['2'])
                    min_time = str(row['4']) # minimum time is the start timestamp of the cascade

                    #min_time = datetime.datetime.strptime(min_time, '%Y-%m-%d %H:%M:%S')
                    rec_time = str(cascade_rt_time[0])
                    rec_date = rec_time[:10]
                    rec_t = rec_time[11:19]
                    record_time = rec_date + ' ' + rec_t
                    #print(record_time)
                    min_time = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')

                    # Store the reaction times for each retweet in the cascade
                    for idx in range(len(cascade_set_1)):
                        rec_time = str(cascade_rt_time[idx])
                        rec_date = rec_time[:10]
                        rec_t = rec_time[11:19]
                        record_time = rec_date + ' ' + rec_t
                        #print(record_time)
                        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                        d1_ts = time.mktime(time_x.timetuple())
                        d2_ts = time.mktime(min_time.timetuple())
                        diff = int(int(d1_ts-d2_ts) / 60)
                        reaction_times.append(diff)
                        x_time.append(diff)
                        time_inst.append(time_x)
                        cascade_time_idx.append(diff)
                        time_values.append(time_x)
                        diff_x = diff

                    time_cascades.append(x_time[len(x_time)-1])
                    size_cascades.append(len(x_time))
                    for idx in range(len(cascade_set_1)):
                        if cascade_set_2[idx] not in node_rt_time:
                            if idx == 0:
                                rt_time = datetime.datetime.strptime(cascade_rt_time[idx], '%Y-%m-%d %H:%M:%S')
                                d1_ts = time.mktime(rt_time.timetuple())
                                d2_ts = time.mktime(min_time.timetuple())
                                diff = int(int(d1_ts-d2_ts) / 60)
                                node_rt_time[cascade_set_2[idx]] = min_time
                                node_rt_time[cascade_set_1[idx]] = rt_time
                                node_reaction_time.append(diff)
                            elif cascade_set_1[idx] not in node_rt_time and idx != 0:
                                rt_time = datetime.datetime.strptime(cascade_rt_time[idx], '%Y-%m-%d %H:%M:%S')
                                rt_time_prev = datetime.datetime.strptime(cascade_rt_time[idx-1], '%Y-%m-%d %H:%M:%S')
                                d1_ts = time.mktime(rt_time.timetuple())
                                d2_ts = time.mktime(rt_time_prev.timetuple())
                                diff = int(int(d1_ts-d2_ts) / 60)
                                node_rt_time[cascade_set_2[idx]] = rt_time_prev
                                node_rt_time[cascade_set_1[idx]] = rt_time
                                node_reaction_time.append(diff)
                            else:
                                diff = 0
                                node_reaction_time.append(diff)
                        else:
                            orig_time = node_rt_time[cascade_set_2[idx]]
                            rt_time = datetime.datetime.strptime(cascade_rt_time[idx], '%Y-%m-%d %H:%M:%S')
                            d1_ts = time.mktime(rt_time.timetuple())
                            d2_ts = time.mktime(orig_time.timetuple())
                            diff = int(int(d1_ts-d2_ts) / 60)
                            node_rt_time[cascade_set_1[idx]] = rt_time
                            node_reaction_time.append(diff)
                            #print(diff)

                    # for i in range(len(node_reaction_time)):
                    #     if node_reaction_time[i] == 0:
                    #         continue
                    #     total_reaction_time.append(node_reaction_time[i])
                    #
                    # newpath = 'E:/paper_ASONAM/figures/'
                    # if not os.path.isdir(newpath):
                    #     os.makedirs(newpath)
                    #
                    # len_rt = len(total_reaction_time)
                    # in_values = sorted(set(total_reaction_time))
                    # list_values = list(total_reaction_time)
                    # in_hist = [list_values.count(x) for x in in_values]
                    # for i in range(len(in_hist)):
                    #     in_hist[i] = (in_hist[i] / len_rt)

                   # Group the cascade retweets based on intervals of size 'freq_cd'
                   #    that is calculated (and hence dynamic)
                    max_offset = max(x_time)
                    cumul = []
                    y_freq = []
                    node_out_degree = {}
                    node_out_degree_idx = {}
                    node_degree_time = []
                    for idx in range(len(cascade_set_1)):
                        node_out_degree[cascade_set_1[idx]] = 0
                        node_out_degree[cascade_set_2[idx]] = 0
                        node_out_degree_idx[cascade_set_2[idx]] = {}
                    for idx in range(len(cascade_set_2)):
                        node_out_degree[cascade_set_2[idx]] += 1

                    for idx in range(len(cascade_set_2)):
                        node_degree_time.append(node_out_degree[cascade_set_2[idx]])

                    list_pd = []
                    for idx in range(len(x_time)):
                        list_pd.append((cascade_set_1[idx], cascade_set_2[idx], time_inst[idx], x_time[idx],
                                        node_reaction_time[idx]))

                    df_cd_size = pd.DataFrame(list_pd, columns=['target', 'source', 'time',
                                                                'Offset', 'reaction_time'])
                    df_cd_size['time'] = pd.to_datetime(df_cd_size['time'])
                    df_cd_size = df_cd_size[df_cd_size['reaction_time'] != 0]
                    index_values = list(range(len(list(df_cd_size['target']))))
                    df_cd_size = df_cd_size.set_index([index_values])
                    freq_cd = str(5*int(math.log10(int(max_offset)))) + 'Min' # Window size
                    window_size.append((5*math.log10(int(max_offset))))
                    df_cd_grouped = df_cd_size.groupby(pd.TimeGrouper(freq=freq_cd, key='time'),
                                                       as_index=False).apply(lambda x: x['Offset'])
                    df_cd_size['period'] = df_cd_grouped.index.get_level_values(0)
                    src_nodes = df_cd_size['source'].tolist()
                    tgt_nodes = df_cd_size['target'].tolist()
                    time_nodes = df_cd_size['time'].tolist()

                    # After grouping, calculate the Hawkes intensity
                    cascade_set_1 = list(df_cd_size['target'])
                    cascade_set_2 = list(df_cd_size['source'])
                    x_time = list(df_cd_size['Offset'])
                    node_reaction_time = list(df_cd_size['reaction_time'])
                    period_val = list(df_cd_size['period'])
                    time_actual = list(df_cd_size['time'])
                    intensity = [0 for i in range(len(cascade_set_1))]
                    for idx in range(len(cascade_set_1)):
                        node_out_degree[cascade_set_1[idx]] = 0
                        node_out_degree[cascade_set_2[idx]] = 0
                        node_out_degree_idx[cascade_set_2[idx]] = {}
                    for idx in range(len(cascade_set_1)):
                        cumul.append(cascade_set_1[idx])
                        cumul.append(cascade_set_2[idx])
                        cumul = list(set(cumul))
                        cumul_len = len(cumul)
                        #print(cascade_set_3[idx], cumul_len)
                        y_freq.append(cumul_len)
                        map_size_time[x_time[idx]] = cumul_len
                    newpath = 'Cascade_figures/07/' + str(cnt)
                    if not os.path.isdir(newpath):
                        os.makedirs(newpath)
                    x_time_new = []
                    node_time_new = []
                    # for idx in range(len(node_degree_time)):
                    #     if node_degree_time[idx] != 0:
                    #         node_time_new.append(node_degree_time[idx])
                    #         x_time_new.append(x_time[idx])
                    print(len(x_time), len(node_degree_time))
                    plt.close()
                    plt.figure(figsize=(20, 12))
                    [y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                    plt.rc('text', usetex=True)
                    plt.rc('axes')
                    plt.rc('font', family='serif', weight='bold')
                    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                    plt.plot(x, y, linewidth=2.5)
                    plt.grid(True)
                    plt.xlabel(r'\textbf{Time offset (minutes)}', fontsize=90)
                    plt.ylabel(r'\textbf{Cascade size}', fontsize=90)
                    plt.xticks(np.arange(0, max(x), 8000))
                    #plt.xticks(np.arange(0,40000,5000))
                    plt.tick_params(axis='x', labelsize=50)
                    plt.tick_params(axis='y', labelsize=50)
                    # plt.tick_params('x', labelsize=20)
                    # plt.tick_params('y', labelsize=20)
                    #plt.title('Cascade lifecyle')
                    #plt.savefig(newpath + '/logistic_curve.png')
                    plt.subplots_adjust(left=0.13, bottom=0.13)
                    plt.show()
                    plt.close()
                   #  #
                    for idx in range(len(cascade_set_2)):
                        node_out_degree[cascade_set_2[idx]] += 1
                        node_out_degree_idx[cascade_set_2[idx]][idx] = node_out_degree[cascade_set_2[idx]]

                    print(len(x_time))
                    time_new = []
                    intensity_new = []
                    max_offset = max(x_time)
                    for idx in range(len(cascade_set_1)):
                        idx_set = math.exp(5*x_time[idx]/max_offset)
                        print(idx_set)
                        # if x_time[math.ceil(idx_set)] <= 5000:
                        #     print(x_time[math.ceil(idx_set)])
                        # else:
                        #     break

                        if idx_set >= 90:
                            idx_set = 90
                        if idx >= idx_set:
                            #t_st = idx - math.ceil(idx_set)
                            t_st = idx - math.ceil(idx_set)
                            while t_st < idx:
                                #out_degree = node_out_degree_idx[cascade_set_2[t_st]][t_st]
                                out_degree = node_out_degree[cascade_set_2[idx]]
                                #intensity[idx] += ((abs(y_freq[idx]-y_freq[idx-1])/y_freq[idx]) * out_degree * (node_reaction_time[t_st])/1000)
                                #intensity[idx] += (out_degree * (node_reaction_time[t_st])/x_time[idx])
                                intensity[idx] += (out_degree * (node_reaction_time[t_st])/x_time[idx])
                                t_st += 1
                            intensity_new.append(intensity[idx])
                            time_new.append(x_time[idx])
                    df_cd_size['intensity'] = intensity
                    df_cd_size['intensity'] = df_cd_size['intensity'].groupby(df_cd_size['period']).transform('sum')
                    df_cd_size = df_cd_size.drop_duplicates('intensity')
                    time_hawkes = df_cd_size['Offset'].tolist()
                    np_hawkes = list(df_cd_size['intensity'])
                    x_time_new = list(df_cd_size['Offset'])
                    time_actual_new = list(df_cd_size['time'])
                    rng_max = int(len(np_hawkes)/6)
                    print(rng_max)

                    plt.close()

                    plt.figure(figsize=(20, 12))
                    #[y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                    plt.rc('text', usetex=True)
                    plt.rc('axes')
                    plt.rc('font', family='serif', weight='bold')
                    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                    plt.plot(time_new, intensity_new, linewidth=2.5)
                    plt.grid(True)
                    plt.xlabel(r'\textbf{Time offset (minutes)}', fontsize=90)
                    plt.ylabel(r'\textbf{Hawkes intensity}', fontsize=90)
                    plt.xticks(np.arange(0, max(time_mean), 8000))
                    #plt.xticks(np.arange(0,40000,5000))
                    plt.tick_params(axis='x', labelsize=50)
                    plt.tick_params(axis='y', labelsize=50)
                    # plt.tick_params('x', labelsize=20)
                    # plt.tick_params('y', labelsize=20)
                    #plt.title('Cascade lifecyle')
                    #plt.savefig(newpath + '/logistic_curve.png')
                    plt.subplots_adjust(left=0.16, bottom=0.16)
                    plt.show()
                    plt.close()


                    plt.figure(figsize=(20, 12))
                    #[y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                    plt.rc('text', usetex=True)
                    plt.rc('axes')
                    plt.rc('font', family='serif', weight='bold')
                    plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                    plt.plot(time_mean, mean_hawkes_fall, linewidth=2.5)
                    plt.grid(True)
                    plt.xlabel(r'\textbf{Time offset (minutes)}', fontsize=90)
                    plt.ylabel(r'\textbf{Hawkes intensity}', fontsize=90)
                    plt.xticks(np.arange(0, max(time_mean), 8000))
                    #plt.xticks(np.arange(0,40000,5000))
                    plt.tick_params(axis='x', labelsize=50)
                    plt.tick_params(axis='y', labelsize=50)
                    # plt.tick_params('x', labelsize=20)
                    # plt.tick_params('y', labelsize=20)
                    #plt.title('Cascade lifecyle')
                    #plt.savefig(newpath + '/logistic_curve.png')
                    plt.subplots_adjust(left=0.16, bottom=0.16)
                    plt.show()
                    plt.close()


                    # break
                else:
                    break
                    # continue

        #
        # f = open('Growth', 'wb')
        # pickle.dump(growth, f)
        # f.close()
        #
        # f = open('Time_difference', 'wb')
        # pickle.dump(time_diff_points, f)
        # f.close()

        # plt.close()
        # plt.figure(1, figsize=(10, 6))
        # n, bins, patches = plt.hist(size_cascades, bins=40, facecolor='b')
        # plt.xlabel('Cascade size', size=35)
        # plt.ylabel('Frequency', size=35)
        # # plt.title('')
        # plt.grid(True)
        # plt.tick_params('y', labelsize=17)
        # plt.tick_params('x', labelsize=17)
        # plt.show()
        # plt.savefig('Cascade_figures/histograms/cascade_size.png')

        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(growth, bins=20, facecolor='g')
        # plt.xlabel('Growth (Decay Size / Steep Size)')
        # plt.ylabel('Frequency')
        # plt.title('')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/histograms/growth.png')
        #
        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(time_diff_points, bins=20, facecolor='g')
        # plt.xlabel('Time difference')
        # plt.ylabel('Frequency')
        # plt.title('')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/histograms/time_diff.png')


        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(first_time, bins=50, facecolor='g')
        # plt.xlabel('Time')
        # plt.ylabel('Frequency')
        # plt.title('')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/trash/first_time_points.png')
        #
        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(steep_points_growth, bins=50, facecolor='g')
        # plt.xlabel('Time')
        # plt.ylabel('Frequency')
        # plt.title('')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/trash/steep_points_growth.png')
        #
        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(steep_points_size, bins=50, facecolor='g')
        # plt.xlabel('Time')
        # plt.ylabel('Frequency')
        # plt.title('')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/trash/peak_indexes.png')
        #
        # # PERFORM MAXIMUM LIKELIHOOD ESTIMATION ON THE STEEP POINTS USING SCIPY OPTIMIZE METHOD
        #
        # def gaussian(x, mu, sig):
        #     return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
        #
        # def poisson(k, lamb):
        #     """poisson pdf, parameter lamb is the fit parameter"""
        #     return (lamb**k/factorial(k)) * np.exp(-lamb)
        #
        # def lognorm(x, mu=0,sigma=1):
        #     a = (math.log(x) - mu)/math.sqrt(2*sigma**2)
        #     p = 0.5 + 0.5*math.erf(a)
        #     return p
        #
        # def negLogLikelihood(params, data):
        #     """ the negative log-Likelohood-Function"""
        #     # print(params)
        #     #lnl = - np.sum(np.log(poisson(data, params[0])))
        #     lnl = - np.sum(np.log(gaussian(data, params[0], params[0])))
        #     #lnl = - np.sum(-np.log(data))
        #     return lnl
        #
        #
        # # minimize the negative log-Likelihood
        #
        # # result = minimize(negLogLikelihood,  # function to minimize
        # #                   x0=np.ones(1),     # start value
        # #                   args=(growth,),      # additional arguments for function
        # #                   method='Powell',   # minimization method, see docs
        # #                   )
        # # # result is a scipy optimize result object, the fit parameters
        # # # are stored in result.x
        # # print(result)
        # #
        # # result = minimize(negLogLikelihood,  # function to minimize
        # #                   x0=np.ones(1000),     # start value
        # #                   args=(time_diff_points,),      # additional arguments for function
        # #                   method='Powell',   # minimization method, see docs
        # #                   )
        # # # result is a scipy optimize result object, the fit parameters
        # # # are stored in result.x
        # # print(result)
        #
        # # plt.close()
        # # plt.figure()
        # # X = np.asarray(time_mean)
        # # Y = np.asarray(mean_hawkes)
        # # X = sm.add_constant(X)
        # # model = sm.OLS(Y, X)
        # # fitted = model.fit()
        # # print(fitted.rsquared)
        # # print(fitted.summary())
        # #
        # # X_pred = np.linspace(X.min(), X.max())
        # # X_pred2 = sm.add_constant(X_pred)
        # # Y_pred = fitted.predict(X_pred2)
        # #
        # # plt.plot(X, Y, '-', color='blue', linewidth=2, label='Original Hawkes curve')
        # # plt.plot(X_pred, Y_pred, '-', color='red', linewidth=2, label='Fitted curve')
        # # plt.xlabel('Time (Offset in minutes')
        # # plt.ylabel('Hawkes intensity')
        # # plt.legend(loc='upper left')
        # # p = newpath + '/fitted_regression.png'
        # # plt.savefig(p)
        #
        # # plt.plot(time_mean, mean_hawkes, '-', color='green', linewidth=2, label='Moving mean - Hawkes amp.')
        # # plt.xlabel('Time (Offset in minutes')
        # # plt.ylabel('Mean Hawkes intensity')
        # # plt.legend(loc='upper left')
        # # p = newpath + '/move_mean.png'
        # # plt.savefig(p)
        #
        # # plt.close()
        # # plt.figure()
        # # n, bins, patches = plt.hist(fitted.resid, bins=30, facecolor='g')
        # # plt.xlabel('Residuals')
        # # plt.ylabel('Frequency')
        # # plt.title('Normalized residuals')
        # # plt.grid(True)
        # # plt.savefig(newpath + '/residuals.png')
        #
        # # n, bins, patches = plt.hist(window_size, 50, facecolor='g')
        # # plt.xlabel('Window sizes')
        # # plt.ylabel('Frequency')
        # # plt.title('Histogram of Window sizes')
        # # plt.grid(True)
        # # plt.savefig('Cascade_figures/window_sizes.png')
        # # print(len(window_size))
        #
        # '''
        # n, bins, patches = plt.hist(node_reaction_time, 50, facecolor='g')
        # plt.xlabel('User activity times (offset in minutes')
        # plt.ylabel('Frequency')
        # plt.title('Histogram of User activity Times')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/trash/user_activity_histogram.png')
        #
        # plt.close()
        # plt.figure()
        # n, bins, patches = plt.hist(node_reaction_time, 50, facecolor='g')
        # plt.xlabel('Reaction times (offset in minutes')
        # plt.ylabel('Frequency')
        # plt.title('Histogram of Reaction Times')
        # plt.grid(True)
        # plt.savefig('Cascade_figures/trash/reaction_histogram.png')
        #
        # plt.close()
        # fig = plt.figure(figsize=(10, 8))
        # [y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
        # #for i in range(len(x)):
        # #    print(x[i], y[i])
        # plt.plot(x, y)
        # #plt.xlim([min(x), max(x)+datetime.timedelta(0, 360)])
        # plt.grid(True)
        #
        # # Create an Axes object.
        # ax = fig.add_subplot(1,1,1) # one row, one column, first plot
        # # Plot the data.
        # ax.scatter(x_time, node_reaction_time, color="red", marker="o")
        # # Add a title.
        # ax.set_title("Time vs node reaction time")
        # # Add some axis labels.
        # ax.set_xlabel("Time offset (minutes)")
        # ax.set_ylabel("Node reaction times")
        # plt.ylim([0,10000])
        # plt.savefig('Cascade_figures/trash/reaction vs time.png')
        #
        # plt.close()
        # plt.figure()
        # plt.plot(x_time, np_hawkes, color='blue')
        # [y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
        # plt.plot(x, y, color='red')
        # plt.grid(True)
        # plt.xlabel('Time offset(minutes)')
        # plt.ylabel('hawkes Intensity')
        # plt.savefig("Cascade_figures/trash/hawkes_curve.png")
        #
        # plt.close()
        # plt.figure()
        # bins = bins[1:]
        # #print(len(bins), len(n))
        # plt.plot(bins, n)
        # plt.grid(True)
        # plt.xlabel('Time offset (minutes)')
        # plt.ylabel('Frequency')
        # plt.savefig('Cascade_figures/trash/reaction_curve.png')
        #
        # plt.close()
        # plt.figure()
        # [y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
        # #for i in range(len(x)):
        # #    print(x[i], y[i])
        # plt.plot(x, y, color='blue')
        # #plt.xlim([min(x), max(x)+datetime.timedelta(0, 360)])
        # plt.grid(True)
        # plt.xlabel('Time offset (minutes')
        # plt.ylabel('Cumulative Cascade Size')
        # plt.title('Cascade lifecyle ')
        # ind = detect_peaks(n,bins, mph=0, mpd=1, show=True)
        # '''

cascades_instance = FormCascadesOrder()
cascades_instance.cascades_form()