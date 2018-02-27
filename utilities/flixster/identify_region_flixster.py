# Find the steep and the inhibition regions of the cascades and mark the corresponding intervals
# in which the points lie.

__author__ = 'ssarka18'

import pandas as pd
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
import math
from math import factorial
from scipy.misc import factorial
import pickle

path = 'cascade_df_06-08_v1+.pickle'
number_intervals = 10

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


class FormCascadesOrder(object):
    def __init__(self):
        self.flag_cascade = []

    def cascades_form(self):
        df = pd.read_pickle(path)
        plotpath = 'flixster_cascades/'
        self.mid_list = list(set(df['mid']))
        self.cascade_dict = {}
        self.x_time_dict = {}
        self.y_freq_dict = {}
        self.steep_times  = {}
        self.inhib_times = {}

        mid_count = 0
        for mid in self.mid_list:
            print("Mid no: ", mid_count)
            mid_count += 1
            self.x_time = []
            self.y_freq = []
            cascade = df[df['mid'] == mid]
            self.cascade_dict[mid] = cascade

            cascade_set_1 = list(cascade['target'])
            cascade_set_2 = list(cascade['source'])
            cascade_rt_time = list(cascade['rt_time'])
            min_time = datetime.datetime.strptime(str(cascade_rt_time[0]), '%Y-%m-%d %H:%M:%S')

            for idx in range(1, len(cascade_set_1)):
                rec_time = str(cascade_rt_time[idx])
                time_x = datetime.datetime.strptime(rec_time, '%Y-%m-%d %H:%M:%S')
                diff = (time_x - min_time).days # Difference in days
                self.x_time.append(diff)

            # Prune the cascades based on length
            # print(diff)
            if diff < 680:
                continue
            for idx in range(len(self.x_time)):
                rec_time = str(cascade_rt_time[idx])
                time_x = datetime.datetime.strptime(rec_time, '%Y-%m-%d %H:%M:%S')
                diff = (time_x - min_time).days  # Difference in days
                if diff > 350:
                    self.steep_times[mid] = self.x_time[idx]
                    break

            for idx in range(len(self.x_time)):
                rec_time = str(cascade_rt_time[idx])
                time_x = datetime.datetime.strptime(rec_time, '%Y-%m-%d %H:%M:%S')
                diff = (time_x - min_time).days  # Difference in days

                if diff > 600:
                    self.inhib_times[mid] = self.x_time[idx]
                    break



            # cumul = []
            # for idx in range(len(cascade_set_1)):
            #     cumul.append(cascade_set_1[idx])
            #     cumul.append(cascade_set_2[idx])
            #     cumul = list(set(cumul))
            #     cumul_len = len(cumul)
            #     self.y_freq.append(cumul_len)

            # if not os.path.isdir(newpath):
            #     os.makedirs(newpath)
            # plt.figure()
            # [y, x] = zip(*sorted(zip(self.y_freq, self.x_time), key=lambda x: x[0]))
            # plt.plot(x, y, lw=3, color='black')
            # plt.grid(True)
            # plt.xlabel('Time offset (days)', fontsize=40)
            # plt.ylabel('Cumulative Cascade Size', fontsize=40)
            # plt.tick_params(axis='x', labelsize=25)
            # plt.tick_params(axis='y', labelsize=25)
            # plt.grid(True)
            # # plt.title('Cascade lifecyle')
            # plt.savefig(plotpath + 'cascade_' + str(mid))
            # plt.close()

        return self.steep_times, self.inhib_times


    def steep_inhib_times(self):
        ''' Compute the steep and inhib times '''
        for mid in self.mid_list:
            cascade = self.cascade_dict[mid]

            # for idx in range(len(self.))


cascades_instance = FormCascadesOrder()
steep_dict, inhib_dict = cascades_instance.cascades_form()
pickle.dump((steep_dict, inhib_dict), open('steep_inhib_times_dict_flixster.pickle', 'wb'))