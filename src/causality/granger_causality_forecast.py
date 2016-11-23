import numpy as np
import matplotlib.pyplot as plt
import pylab
import pickle
import os
import pandas as pd
import statsmodels.tsa.api as sta
import math
import statsmodels.tsa.stattools as sts
from datetime import datetime
import time
import math
import statistics as st
import statsmodels.stats.stattools as ssts
import seaborn
import itertools
from pandas.stats.api import ols
import statsmodels.api as sm

# Using statsmodels granger causality api.
# Test for co-integration between the time series is not implemented here.

def test_stationarity(timeseries):

    #Determing rolling statistics
    rolmean = pd.rolling_mean(timeseries, window=12)
    rolstd = pd.rolling_std(timeseries, window=12)

    #Plot rolling statistics:
    # plt.close()
    # orig = plt.plot(timeseries, color='blue',label='Original')
    # mean = plt.plot(rolmean, color='red', label='Rolling Mean')
    # std = plt.plot(rolstd, color='black', label = 'Rolling Std')
    # plt.legend(loc='best')
    # plt.title('Rolling Mean & Standard Deviation')
    # plt.show()

    #Perform Dickey-Fuller test:
    # print('Results of Dickey-Fuller Test:')
    dftest = sts.adfuller(timeseries, autolag='AIC')
    dfoutput = pd.Series(dftest[0:4], index=['Test Statistic','p-value','#Lags Used','Number of Observations Used'])
    # for key,value in dftest[4].items():
    #     dfoutput['Critical Value (%s)'%key] = value
    # print(dfoutput)
    return dfoutput[2]

if __name__ == '__main__':
    measure_file_path = 'F://Inhibition//VAR_causality//data_files//measure_time_series//'
    steep_inhib_times = pickle.load(open('F://Inhibition//VAR_causality//data_files//steep_inhib_times.pickle', 'rb'))

    dependent_variables = []
    # print(cascade_int_time)

    # Load the measure files
    cnt_measures = 0
    model_df = pd.DataFrame()
    control_variables = {}
    dependent_variables = {}

    measures = []
    model_order = 10
    f_list_null = []
    p_list_null = []

    f_list_full = []
    p_list_full = []
    measure_time_df = []

    m_steps = 5
    for files in os.listdir(measure_file_path):
        full_path = measure_file_path + '//' + files
        mname = files[:len(files)-7]
        measure_time_df.append(pickle.load(open(full_path, 'rb')))
        measures.append(mname)

    mae_model = {}
    statistic_values = {}
    critical_values = {}
    p_values = {}
    granger_cause_count = {}

    for L in range(0, len(range(len(measures))) + 1):
        for subset in itertools.combinations(range(len(measures)), L):
            if len(subset) > 1 and len(subset) < len(measures):
                continue
            num_measures = len(subset)
            if num_measures == 0:
                continue
            measure_series = []
            cascade_df_feat = pd.DataFrame()
            cnt_mids = 0
            measures_causality = []
            measures_string = ''
            for idx_sub in range(len(subset)):
                measures_causality.append(measures[subset[idx_sub]])
                measures_string += (measures[subset[idx_sub]] + ' + ')
            measures_string = measures_string[:len(measures_string) - 3]
            for mid in measure_time_df[0]:
                # print(mid)
                steep_time = pd.to_datetime(steep_inhib_times[mid]['steep'])
                inhib_time = pd.to_datetime(steep_inhib_times[mid]['decay'])
                # print(steep_time)
                # print(inhib_time)
                # Combine all the features in a dataframe for sorting them by time.
                cascade_df = measure_time_df[0][mid]
                cascade_df['time'] = pd.to_datetime(cascade_df['time'])
                cascade_df_feat = cascade_df

                # print(cascade_df_feat)
                for idx in range(num_measures):
                    # if measures[subset[idx]] == 'bw' or measures[subset[idx]] == 'cc':
                    cascade_df_feat[measures[subset[idx]]] = measure_time_df[subset[idx]][mid][measures[subset[idx]]]
                    # cascade_df_feat[measures[subset[idx]]] = measure_time_df[subset[idx]][mid]['nbr_deg']
                cascade_df_feat = cascade_df_feat.sort('time')

                ### THIS IS IMPORTANT - CLIPPING THE CASCADE BASED ON STEEP AND INHIBITION REGIONS
                cascade_df_feat = cascade_df_feat[cascade_df_feat['time'] < pd.to_datetime(inhib_time)]

                # print(cascade_df_feat)
                time_series = list(cascade_df_feat['time'])
                measure_series_all = []
                for idx in range(num_measures):
                    measure_series_all.append(list(cascade_df_feat[measures[subset[idx]]]))

                Y_act = [0 for idx in range(len(time_series))]

                for idx in range(1, len(time_series)):
                    rt_time = str(time_series[idx])
                    rt_date = rt_time[:10]
                    rt_t = rt_time[11:19]
                    record_time = rt_date + ' ' + rt_t
                    time_x = datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                    cur_time = time.mktime(time_x.timetuple())

                    rt_time = str(time_series[idx-1])
                    rt_date = rt_time[:10]
                    rt_t = rt_time[11:19]
                    record_time = rt_date + ' ' + rt_t
                    time_x = datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                    prev_time = time.mktime(time_x.timetuple())

                    Y_act[idx] = (cur_time - prev_time)/60

                X = [[] for i in range(len(subset))]
                Y = []
                time_new = []
                sum_all_points = [0 for idx_points in range(len(subset))]
                num_points = 0
                # Take the mean values of reshares at the same time
                for idx in range(len(Y_act)):
                    if Y_act[idx] == 0 and idx > 0:
                        for idx_sub in range(len(subset)):
                            sum_all_points[idx_sub] += measure_series_all[idx_sub][idx]
                        num_points += 1
                        continue
                    if idx > 0 and Y_act[idx-1] == 0:
                        for idx_sub in range(len(subset)):
                            sum_all_points[idx_sub] += X[idx_sub][len(X[idx_sub])-1]
                            num_points += 1
                            # print(sum_all_points[idx_sub], num_points)
                            X[idx_sub][len(X[idx_sub])-1] = (sum_all_points[idx_sub]/num_points)
                        sum_all_points = [0 for idx_points in range(len(subset))]
                        num_points = 0

                    for idx_sub in range(len(subset)):
                        X[idx_sub].append(measure_series_all[idx_sub][idx])
                    Y.append(Y_act[idx])
                    time_new.append(time_series[idx])


                cascade_ts_df = pd.DataFrame()
                cascade_ts_df['time_diff'] = Y
                cascade_ts_df.index = time_new
                for idx_sub in range(len(subset)):
                    cascade_ts_df[measures[subset[idx_sub]]] = X[idx_sub]
                # print(cascade_ts_df)
                # Step 1: Setup the VAR model to set the lag order based on AIC and BIC
                # try:
                VAR_model = sta.VAR(cascade_ts_df)
                VAR_results = VAR_model.fit(maxlags=model_order+1, ic='aic')
                model_order = VAR_results.k_ar

                Y_diff = []
                X_diff = [[] for i in range(len(subset))]
                for idx in range(len(Y)):
                    if np.isnan(Y[idx]):
                       continue
                    Y_diff.append(Y[idx])
                    for idx_sub in range(len(subset)):
                        X_diff[idx_sub].append(X[idx_sub][idx])

                X = np.asarray(X_diff)
                Y = np.asarray(Y_diff)
                #
                # measures_causality = []
                # measures_string = ''
                # for idx_sub in range(len(subset)):
                #     measures_causality.append(measures[subset[idx_sub]])
                #     measures_string += (measures[subset[idx_sub]] + ' + ')
                # measures_string = measures_string[:len(measures_string) - 3]
                #
                # if measures_string not in mae_model:
                #     mae_model[measures_string] = {}
                #
                # cascade_VAR_df = pd.DataFrame()
                # cascade_VAR_df['time_diff'] = Y
                # for idx_sub in range(len(subset)):
                #     cascade_VAR_df[measures[subset[idx_sub]]] = X[idx_sub, :]

                # print(cascade_VAR_df)

                # Granger causality regression
                X = np.asarray(X)
                Y = np.asarray(Y)
                Z = np.transpose(np.vstack((Y, X)))
                try:
                    granger_results = sts.grangercausalitytests(Z, maxlag=model_order)
                except:
                    continue
                for idx_lag in range(1, model_order + 1):
                    for objects in granger_results[idx_lag]:
                        # Get the OLS regression results
                        if type(objects) is list:
                            results_OLS_granger = objects[1]  # Granger results for the full model
                            print(results_OLS_granger.summary())
                #
                # for m_s in range(1, m_steps+1):
                #     if m_s not in mae_model[measures_string]:
                #         mae_model[measures_string][m_s] = []
                #
                #     x = sm.add_constant(cascade_VAR_df[measures_causality])
                #     X_train = x.iloc[:len(Y) - m_s]
                #     Y_train = Y[:len(Y) - m_s]
                #
                #     test_data = x.iloc[len(Y)-m_s:]
                #     test_values = Y[len(Y)-1-m_s:]
                #     try:
                #         # print(X_train)
                #         OLSmodel = sm.OLS(Y_train, X_train)
                #         res_model = OLSmodel.fit()
                #         estimated_values = res_model.predict(test_data)
                #         diff = abs(estimated_values - test_values)
                #         mae_model[measures_string][m_s].append(diff)
                #
                #         # post_label = res.predict(test_data)
                #         # some issue with the plots....
                #         # VAR_results.plot()
                #         # plt.show()
                #     except ValueError:
                #         continue

                # # print(cascade_ts_df)
                # try:
                #     cascade_ts_df = cascade_ts_df.fillna(0)
                #     VAR_model = sta.VAR(cascade_ts_df)
                #     VAR_results = VAR_model.fit(5)
                #     # print(VAR_results.summary())
                #     model_order = VAR_results.k_ar
                #     causality_results = VAR_results.test_causality('time_diff', measures_causality, kind='wald', verbose=False)
                #     if measures_string not in statistic_values:
                #         statistic_values[measures_string] = []
                #         critical_values[measures_string] = []
                #         p_values[measures_string] = []
                #         granger_cause_count[measures_string] = 0
                #
                #     statistic_values[measures_string].append(causality_results['statistic'])
                #     critical_values[measures_string].append(causality_results['crit_value'])
                #     p_values[measures_string].append(causality_results['pvalue'])
                #     if causality_results['conclusion'] != 'reject':
                #         granger_cause_count[measures_string] += 1
                # except:
                #     # print(cnt_mids)
                #     continue

                    # exit()

                cnt_mids += 1
                print(cnt_mids)
                if cnt_mids > 0:
                    break
    # print(cnt_mids)
    # pickle.dump(mae_model, open('..//..//data_files//results_granger//11_23//mae_OLS.pickle', 'wb'))
    # pickle.dump(statistic_values, open('../../data_files//results_granger//11_23//statistics.pickle', 'wb'))
    # pickle.dump(critical_values,
    #             open('../../data_files//results_granger//11_23//critical.pickle', 'wb'))
    # pickle.dump(p_values,
    #             open('../../data_files//results_granger//11_23//p_values.pickle', 'wb'))
    # pickle.dump(granger_cause_count,
    #             open('../../data_files//results_granger//11_23//cause_count.pickle', 'wb'))
    #

