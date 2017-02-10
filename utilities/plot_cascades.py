__author__ = 'ssarka18'

# FIT A LOGISTIC CURVE TO THE TYPE II CASCADES - THE MAIN OBJECTIVE IS TO SEE HOW THE TYPE II CASCADES DIFFER WITHIN THEMSELVES
# AND THEN SEE HOW PREDICTION PERFORMS FOR THESE TYPES

import pandas as pd
import numpy as np
import time
import datetime
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
from scipy.stats.distributions import t
import os
from mpl_toolkits.axes_grid.inset_locator import inset_axes

path = '../data/rt_df.csv'

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=12)#x_labelsize)
    subax.yaxis.set_tick_params(labelsize=12)#y_labelsize)
    return subax


class FitCascades(object):
    def logistic_function(self, x, k, x0):
        return (691/(1+np.exp(-k*(x-x0))))

    def concave_up_increase(self, x, k, l):
        return l*(x**k)

    def straight_line(self, x, a, b):
        return (b/a)*x

    def cascades_fit(self):
        flag_scanned = {}
        df = pd.read_csv(path)
        cnt = 0
        size_cascades = []
        cnt_300 = 0

        # path_files = ''
        # filename_list = list(os.listdir(path_files))
        for index, row in df.iterrows():
            if str(row['3']) not in flag_scanned:
                flag_scanned[str(row['3'])] = 1
                x_time = []
                reaction_times = []

                map_size_time = {}
                cascade_set = df[(df['3'] == row['3'])]
                cnt += 1
                cascade_set_1 = list(cascade_set['0'])
                cascade_set_2 = list(cascade_set['1'])

                len_cascade = len(list(set(cascade_set_1 + cascade_set_2)))
                #print(len_cascade)
                size_cascades.append(len_cascade)
                if len_cascade <= 300: # DO NOT CONSIDER CASCADES OF SIZE LESS THAN 300
                    continue
                cnt_300 += 1
                #print(len_cascade)
                print(cnt_300, cnt)
                # 1883
                # 2424
                # 1306
                # 141
                cd = pd.DataFrame()
                if cnt == 2424:
                    # cd['source'] = cascade_set_2
                    # cd['target'] = cascade_set_1
                    # cd.to_csv('cascade_temp.csv')
                    #continue
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
                    prev_time = datetime.datetime.strptime('1992-01-01', '%Y-%m-%d')
                    time_plot = []
                    shares_plot = []
                    temp_shares = 1
                    for idx in range(len(cascade_set_1)):
                        rec_time = str(cascade_rt_time[idx])
                        rec_date = rec_time[:10]
                        rec_t = rec_time[11:19]
                        record_time = rec_date + ' ' + rec_t
                        #print(record_time)
                        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                        time_date = datetime.datetime.strptime(rec_date, '%Y-%m-%d')
                        d1_ts = time.mktime(time_x.timetuple())
                        d2_ts = time.mktime(min_time.timetuple())
                        diff = int(int(d1_ts-d2_ts) / 60)
                        reaction_times.append(diff)
                        x_time.append(diff)
                        if time_date != prev_time:
                            shares_plot.append(temp_shares)
                            time_plot.append(prev_time)
                            temp_shares = 1
                            prev_time = time_date
                        else:
                            temp_shares += 1


                    cumul = []
                    y_freq = []
                    for idx in range(len(cascade_set_1)):
                        cumul.append(cascade_set_1[idx])
                        cumul.append(cascade_set_2[idx])
                        cumul = list(set(cumul))
                        cumul_len = len(cumul)
                        #print(cascade_set_3[idx], cumul_len)
                        y_freq.append(cumul_len)
                        map_size_time[x_time[idx]] = cumul_len

                    orig_x = np.array(x_time)
                    orig_y = np.array(y_freq)

                    max_x = max(orig_x)
                    max_y = max(orig_y)
                    # pars, pcov = curve_fit(self.logistic_function, orig_x, orig_y)
                    # print(pars)
                    alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

                    # n = len(orig_y)    # number of data points
                    # p = len(pars) # number of parameters

                    # dof = max(0, n - p) # number of degrees of freedom

                    # student-t value for the dof and confidence level
                    # tval = t.ppf(1.0-alpha/2., dof)

                    # for i, p,var in zip(range(n), pars, np.diag(pcov)):
                    #     sigma = var**0.5
                        # print('p{0}: {1} [{2}  {3}]'.format(i, p,
                        #                               p - sigma*tval,
                        #                               p + sigma*tval))
                    fig = plt.figure(figsize=(12,8))
                    ax = fig.add_subplot(111)
                    plt.plot(orig_x, orig_y, 'b', linewidth=1.5)
                    # yfit = self.concave_up_increase(orig_x, 0.45, 8.35)
                    # yfit = self.straight_line(orig_x, max_x, max_y)
                    # yfit = self.logistic_function(orig_x, 0.000221, 10000)
                    # plt.plot(orig_x, yfit,'bo ', color='black', linewidth=2)

                    lower = []
                    upper = []
                    # for p,var in zip(pars, np.diag(pcov)):
                    #     sigma = var**0.5
                    #     lower.append(p - sigma*tval)
                    #     upper.append(p + sigma*tval)

                    # yfit = self.straight_line(orig_x, *lower)
                    # plt.plot(orig_x, yfit, '--', color='black')
                    # yfit = self.straight_line(orig_x, *upper)

                    # plt.rc('text', usetex=True)
                    # plt.rc('axes')
                    # plt.rc('font', family='arial')
                    # plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

                    # plt.plot(orig_x, yfit, '--', color='black', linewidth=3)
                    hfont = {'fontname': 'Arial'}
                    # plt.legend(['data', 'curve fit'], loc='best', fontsize=27)
                    plt.xlabel('Time offset (min.)', fontsize=45, **hfont)
                    plt.ylabel('Cascade   size', fontsize=45, **hfont)
                    plt.xticks(np.arange(0, max(orig_x), 8000))
                    plt.tick_params('x', labelsize=25)
                    plt.tick_params('y', labelsize=25)
                    plt.grid(True)
                    plt.subplots_adjust(left=0.13, bottom=0.13)
                    # plt.show()
                    dir = 'cascades/' + str(cnt)
                    # if not os.path.exists(dir):
                    #     os.makedirs(dir)
                    # plt.savefig(dir + '/' + str(cnt) + '.png')
                    # residuals_logistic = orig_y - self.concave_up_increase(orig_x, pars[0], pars[1])
                    # fres_logistic = sum(residuals_logistic**2)
                    # print(fres_logistic)

                    subpos = [0.5, 0.1, 0.45, 0.35]
                    subax1 = add_subplot_axes(ax, subpos)


                    hfont = {'fontname': 'Arial'}
                    width = 0.35
                    # ind = np.arange(len(data_mean))  # the x locations for the groups
                    ind = np.arange(len(shares_plot))
                    # fig = plt.figure(figsize=(12, 8))
                    # ax = fig.add_subplot(111)
                    ## the bars
                    # rects1 = ax.bar(ind, data_mean, width,
                    #                 color='#C0C0C0')
                    rects1 = subax1.bar(ind, shares_plot, width,
                                    color='#000000')  # axes and labels
                    subax1.set_xlim(-width, len(ind) + width)
                    # ax.set_ylim(87, 95)
                    subax1.set_ylabel('Shares per day', size=15)
                    subax1.set_xlabel('Days since start', size=15)
                    # ax.set_title('Scores by group and gender')
                    # xTickMarks = titles
                    # subax1.set_xticks(ind, fontsize=10)
                    # xtickNames = ax.set_xticklabels(xTickMarks)
                    # plt.setp(xtickNames, rotation=45, fontsize=5)
                    subax1.grid(True)

                    # subax1.s xticks(size=12)
                    # subax1.yticks(size=12)
                    # in_ax.subplots_adjust(left=0.13, bottom=0.25)
                    ## add a legend
                    # ax.legend( (rects1[0], ('Men', 'Women') )

                    plt.show()

                    #
                    # orig_x = np.array(orig_x)
                    # orig_y = np.array(orig_y)
                    # residuals_straight = orig_y - self.straight_line(orig_x, max_x, max_y)
                    # residuals_concave = orig_y - self.concave_up_increase(orig_x, 0.45, 8.35)
                    # residuals_logistic = orig_y - self.logistic_function(max_y, orig_x, 0.000221, 10000)
                    # fres_concave = sum(residuals_concave**2)
                    # fres_logistic = sum(residuals_logistic**2)
                    # fres_straight = sum(residuals_straight**2)
                    # val = [fres_concave, fres_logistic, fres_straight]
                    # ind_min = val.index(min(val))
                    #
                    # # plt.plot(orig_x, orig_y, 'b-')
                    # # plt.plot(orig_x, self.straight_line(orig_x, max_x, max_y), 'b-', color='cyan')
                    # # plt.plot(orig_x, self.concave_up_increase(orig_x, 0.45, 8.35), 'b-', color='green')
                    # # plt.plot(orig_x, self.logistic_function(max_y, orig_x, 0.000221, 10000), 'b-', color='black')
                    # # plt.show()
                    #
                    # with open('types_cascades_08.txt', 'a') as fp:
                    #     a = csv.writer(fp, lineterminator= '\n')
                    #     all = []
                    #     data = []
                    #     data.append(cnt)
                    #     data.append(ind_min)
                    #     all.append(data)
                    #     a.writerows(all)
                    # fp.close()

                    break
                else:
                    # break
                    continue

cascades_instance = FitCascades()
cascades_instance.cascades_fit()