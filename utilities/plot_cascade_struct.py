# Find the steep and the inhibition regions of the cascades

__author__ = 'ssarka18'

import pandas as pd
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks
import math
from math import factorial
from scipy.misc import factorial
from scipy.optimize import minimize
import csv
import os
import statsmodels.api as sm
import pickle

path = '../data//rt_df.csv'

# Savitzky Golay filtering for smoothing a curve
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        print("error")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


reaction_global = []


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
                if cnt <=5000:

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
                                reaction_global.append(diff)
                            elif cascade_set_1[idx] not in node_rt_time and idx != 0:
                                rt_time = datetime.datetime.strptime(cascade_rt_time[idx], '%Y-%m-%d %H:%M:%S')
                                rt_time_prev = datetime.datetime.strptime(cascade_rt_time[idx-1], '%Y-%m-%d %H:%M:%S')
                                d1_ts = time.mktime(rt_time.timetuple())
                                d2_ts = time.mktime(rt_time_prev.timetuple())
                                diff = int(int(d1_ts-d2_ts) / 60)
                                node_rt_time[cascade_set_2[idx]] = rt_time_prev
                                node_rt_time[cascade_set_1[idx]] = rt_time
                                node_reaction_time.append(diff)
                                reaction_global.append(diff)
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
                            reaction_global.append(diff)
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
                   #  hfont = {'fontname': 'Arial'}
                   #  csfont = {'fontname': 'Arial'}
                   #  # print(len(x_time), len(node_degree_time))
                   #  plt.close()
                   #  plt.figure(figsize=(12, 8))
                   #  [y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                   #  # plt.rc('text', usetex=True)
                   #  # plt.rc('axes')
                   #  # plt.rc('font', family='Arial')
                   #  # plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                   #  plt.plot(x, y, linewidth=2, color='blue')
                   #  plt.grid(True)
                   #  plt.xlabel('Time offset (min.)', fontsize=40, **hfont)
                   #  plt.ylabel('Cascade size', fontsize=40, **csfont)
                   #  plt.xticks(np.arange(0, max(x), 8000))
                   #  #plt.xticks(np.arange(0,40000,5000))
                   #  plt.tick_params(axis='x', labelsize=25)
                   #  plt.tick_params(axis='y', labelsize=25)
                   #  # plt.tick_params('x', labelsize=20)
                   #  # plt.tick_params('y', labelsize=20)
                   #  #plt.title('Cascade lifecyle')
                   #  #plt.savefig(newpath + '/logistic_curve.png')
                   #  plt.subplots_adjust(left=0.13, bottom=0.13)
                   #  plt.show()
                   #  plt.close()
                   # #  #
                   #  for idx in range(len(cascade_set_2)):
                   #      node_out_degree[cascade_set_2[idx]] += 1
                   #      node_out_degree_idx[cascade_set_2[idx]][idx] = node_out_degree[cascade_set_2[idx]]
                   #
                   #  # print(len(x_time))
                   #  time_new = []
                   #  intensity_new = []
                   #  max_offset = max(x_time)
                   #  for idx in range(len(cascade_set_1)):
                   #      idx_set = math.exp(5*x_time[idx]/max_offset)
                   #      # print(idx_set)
                   #      # if x_time[math.ceil(idx_set)] <= 5000:
                   #      #     print(x_time[math.ceil(idx_set)])
                   #      # else:
                   #      #     break
                   #
                   #      if idx_set >= 90:
                   #          idx_set = 90
                   #      if idx >= idx_set:
                   #          #t_st = idx - math.ceil(idx_set)
                   #          t_st = idx - math.ceil(idx_set)
                   #          while t_st < idx:
                   #              #out_degree = node_out_degree_idx[cascade_set_2[t_st]][t_st]
                   #              out_degree = node_out_degree[cascade_set_2[idx]]
                   #              #intensity[idx] += ((abs(y_freq[idx]-y_freq[idx-1])/y_freq[idx]) * out_degree * (node_reaction_time[t_st])/1000)
                   #              #intensity[idx] += (out_degree * (node_reaction_time[t_st])/x_time[idx])
                   #              intensity[idx] += (out_degree * (node_reaction_time[t_st])/x_time[idx])
                   #              t_st += 1
                   #          if intensity[idx] < 0:
                   #              continue
                   #          intensity_new.append(intensity[idx])
                   #          time_new.append(x_time[idx])
                   #  df_cd_size['intensity'] = intensity
                   #  df_cd_size['intensity'] = df_cd_size['intensity'].groupby(df_cd_size['period']).transform('sum')
                   #  df_cd_size = df_cd_size.drop_duplicates('intensity')
                   #  time_hawkes = df_cd_size['Offset'].tolist()
                   #  np_hawkes = list(df_cd_size['intensity'])
                   #  x_time_new = list(df_cd_size['Offset'])
                   #  time_actual_new = list(df_cd_size['time'])
                   #  rng_max = int(len(np_hawkes)/6)
                   #  # print(rng_max)
                   #  if rng_max >= 30:
                   #      rng_max = 30
                   #  if rng_max % 2 == 0:
                   #      rng_max += 1
                   #  try:
                   #      np_hawkes = savitzky_golay(np_hawkes, rng_max, 2)
                   #  except:
                   #      continue
                   #  np_hawkes = np_hawkes.tolist()
                   #
                   #  np_hawkes_temp = []
                   #  for idx_np in range(len(np_hawkes)):
                   #      if np_hawkes[idx_np] < 0:
                   #          np_hawkes_temp.append(0)
                   #          continue
                   #      np_hawkes_temp.append(np_hawkes[idx_np])
                   #
                   #  np_hawkes = np_hawkes_temp
                   #  # CALCULATE THE MOVING MEANS FOR EACH TIME POINT ALONG THE CURVE OF HAWKES INTESNITY
                   #  mean_hawkes = []
                   #  time_mean = []
                   #  rng_max = 5
                   #  for idx in range(len(np_hawkes)):
                   #      if idx > rng_max:
                   #          #time_mean.append(time_hawkes[idx])
                   #          mean_hawkes.append(np.mean(np_hawkes[idx-rng_max:idx]))
                   #  # CALCULATE THE MOVING MEANS FOR EACH TIME POINT ALONG THE CURVE OF HAWKES INTESNITY
                   #  # FOR DETECTING THE LOCAL MINIMAS
                   #  mean_hawkes_fall = []
                   #  for idx in range(len(np_hawkes)):
                   #      if idx > rng_max and idx + rng_max < len(np_hawkes):
                   #          if time_hawkes[idx] < 28000 and time_hawkes[idx] > 23000:
                   #              time_mean.append(time_hawkes[idx])
                   #              mean_hawkes_fall.append(np.mean(np_hawkes[idx-rng_max:idx+rng_max])-600)
                   #              continue
                   #          val = np.mean(np_hawkes[idx-rng_max:idx+rng_max])-200
                   #          if val < 0:
                   #              continue
                   #          time_mean.append(time_hawkes[idx])
                   #          mean_hawkes_fall.append(np.mean(np_hawkes[idx-rng_max:idx+rng_max])-200)
                   #
                   #  newpath = 'Cascade_figures/07/' + str(cnt)
                   #  if not os.path.exists(newpath):
                   #      os.makedirs(newpath)
                   #
                   #  plt.close()
                   #
                   #  print(np_hawkes)
                   #  plt.figure(figsize=(12, 8))
                   #  #[y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                   #  # plt.rc('text', usetex=True)
                   #  # plt.rc('axes')
                   #  # plt.rc('font', family='helvetica')
                   #  # plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                   #  # plt.plot(x_time_new, np_hawkes, linewidth=2, color='#000000')
                   #  # plt.plot(x_time_new, np_hawkes, linewidth=2, color='#000000')
                   #  # plt.stem(x_time_new, np_hawkes, linefmt='b-', markerfmt='bo', markerfacecolor='black',
                   #  #          markeredgecolor='#000000')
                   #  markerline, stemlines, baseline = plt.stem(time_new, intensity_new, '-.')
                   #  plt.setp(stemlines, 'color', 'black')
                   #  plt.setp(markerline, 'markerfacecolor', 'black')
                   #  plt.grid(True)
                   #  plt.xlabel('Time offset (min.)', fontsize=50, **hfont)
                   #  plt.ylabel('Hawkes Intensity', fontsize=50, **hfont)
                   #  plt.xticks(np.arange(0, max(x_time_new), 8000))
                   #  #plt.xticks(np.arange(0,40000,5000))
                   #  plt.tick_params(axis='x', labelsize=25)
                   #  plt.tick_params(axis='y', labelsize=25)
                   #  # plt.tick_params('x', labelsize=20)
                   #  # plt.tick_params('y', labelsize=20)
                   #  #plt.title('Cascade lifecyle')
                   #  #plt.savefig(newpath + '/logistic_curve.png')
                   #  plt.subplots_adjust(left=0.16, bottom=0.16)
                   #  plt.show()
                   #  plt.close()
                   #
                   #
                   #  plt.figure(figsize=(12, 8))
                   #  #[y, x] = zip(*sorted(zip(y_freq, x_time), key=lambda x: x[0]))
                   #  # plt.rc('text', usetex=True)
                   #  # plt.rc('axes')
                   #  # plt.rc('font', family='arial', weight='bold')
                   #  # plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
                   #  # plt.plot(time_mean, mean_hawkes_fall, linewidth=2, color='#000000')
                   #  # # plt.stem(time_mean, mean_hawkes_fall, linefmt='b-', markerfmt='bo', markerfacecolor='#000000',
                   #  # #          markeredgecolor='#000000')
                   #  markerline, stemlines, baseline = plt.stem(time_mean, mean_hawkes_fall, '-.')
                   #  plt.setp(stemlines, 'color', 'black')
                   #  plt.setp(markerline, 'markerfacecolor', 'black')
                   #  plt.grid(True)
                   #  plt.xlabel('Time offset (min.)', fontsize=50, **hfont)
                   #  plt.ylabel('Hawkes Intensity', fontsize=50, **hfont)
                   #  plt.xticks(np.arange(0, max(x_time_new), 8000))
                   #  #plt.xticks(np.arange(0,40000,5000))
                   #  plt.tick_params(axis='x', labelsize=25)
                   #  plt.tick_params(axis='y', labelsize=25)
                   #  # plt.tick_params('x', labelsize=20)
                   #  # plt.tick_params('y', labelsize=20)
                   #  #plt.title('Cascade lifecyle')
                   #  #plt.savefig(newpath + '/logistic_curve.png')
                   #  plt.subplots_adjust(left=0.16, bottom=0.16)
                   #  plt.show()
                   #  plt.close()

                    #
                    # #
                    # peak_values = detect_peaks(np_hawkes, mean_hawkes, time_hawkes, newpath, count=cnt, show=True)
                    #peak_values_down = detect_peaks(np_hawkes, mean_hawkes_fall, time_hawkes, newpath, count=cnt, edge='falling', show=True)
                    # # #print(peak_values_down)
                    # # peak_indexes_fall = []
                    # # peak_amp_fall = []
                    # # peak_min = 10000000
                    # # peak_min_time = 0
                    #
                    # peak_indexes_rise = []
                    # peak_amp_rise = []
                    #
                    # max_peak = -1
                    # max_x = 0
                    # max_x_2 = 0
                    # for i in range(len(peak_values)):
                    #     ind = np_hawkes.index(peak_values[i])
                    #     if True: #x_time_new[ind] < (peak_down):
                    #         if peak_values[i] > max_peak and x_time_new[ind] < 10000:
                    #             max_peak = peak_values[i]
                    #             max_x = x_time_new[ind]
                    #             time_max = time_actual_new[ind]
                    #         # peak_indexes_rise.append(x_time_new[ind])
                    #         # y_ind = x_time.index(x_time_new[ind])
                    #         # peak_amp_rise.append(y_freq[y_ind])
                    #         # #first_time.append(y_freq[y_ind])
                    #         # steep_points_growth.append(y_freq[y_ind]/base_time)
                    #         # steep_points_size.append(y_freq[y_ind])
                    #
                    # if max_x == 0:
                    #     continue
                    #
                    # peak_indexes_rise.append(max_x)
                    # y_ind = x_time.index(max_x)
                    # peak_amp_rise.append(y_freq[y_ind])
                    # #first_time.append(y_freq[y_ind])
                    # #steep_points_growth.append(y_freq[y_ind]/base_time)
                    # steep_points_size.append(y_freq[y_ind])
                    # steep_size = y_freq[y_ind]
                    #

                    # for i in range(len(peak_values_down)):
                    #     ind = np_hawkes.index(peak_values_down[i])
                    #     x_value = x_time_new[ind]
                    #     #print(peak_min, peak_values_down[i])
                    #     #print(map_size_time[x_value], steep_size)
                    #     #if x_time_new[ind] < 15000:
                    #     # if (map_size_time[x_value] >= 1.5*steep_size) \
                    #     #      and ((x_value-max_x) >= 500) and (peak_min > peak_values_down[i])\
                    #     #         and x_value < 20000:
                    #     if x_value > max_x:
                    #         peak_min = peak_values_down[i]
                    #         peak_min_time = x_time_new[ind]
                    #         time_min = time_actual_new[ind]
                    #         #break
                    #
                    #     #peak_indexes_fall.append(x_time_new[ind])
                    # if peak_min_time == 0:
                    #     continue
                    # try:
                    #     #print(len(peak_indexes_fall))
                    #     peak_down = peak_min_time
                    #
                    #     #print("Peak min. time", time_min)
                    #     #peak_down = max(peak_indexes_fall)
                    # except:
                    #     continue

                    # for i in range(len(peak_values_down)):
                    #     ind = np_hawkes.index(peak_values_down[i])
                    #     if x_time_new[ind] > max_x: # peak_down:
                    #         y_ind = x_time.index(x_time_new[ind])
                    #         peak_amp_fall.append(y_freq[y_ind])
                    #         inhibit_size = y_freq[y_ind]
                    #         growth.append(inhibit_size/steep_size)
                    #         time_diff_points.append(x_time_new[ind] - max_x)




                    # break
                else:
                    break
                    # continue

cascades_instance = FormCascadesOrder()
cascades_instance.cascades_form()
#
rec_new = []
rec_new.extend(reaction_global)
rec_new.extend(reaction_global)
hfont = {'fontname': 'Arial'}
in_values = sorted(set(rec_new))
list_degree_values = list(rec_new)

in_hist = [list_degree_values.count(x)/len(rec_new) for x in in_values]
plt.figure()
plt.semilogx(in_values, in_hist, 'o', basex=2, color='black')
plt.ylim([-0.001, 0.05])
plt.xlim([0, 1000])
plt.xlabel('Reaction time', fontsize=50, **hfont)
plt.ylabel('Density', fontsize=50, **hfont)
plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)
plt.grid(True)
plt.show()
