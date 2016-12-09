# Calculate the log-likelihood of some functions which are assumed from the distributions.
# Find the log-likelihood which is the maximum among them.
# First use a fit function to guess the parameters, then use the same function to measure
# log likelihood.

# Make a function of both age of node and degree of the source node.
# the problem is the probabiltiy distribution in the log likelihood function
# is multivariate - need to resolve this !!!! - *********************************

import matplotlib.pyplot as plt
import pandas as pd
from scipy.misc import factorial
from scipy.optimize import minimize
import statsmodels.api as sm
import numpy as np
import pickle
from scipy.optimize import curve_fit
import math
from scipy.optimize import minimize
import scipy.stats as stats
from matplotlib import gridspec
from pylab import *

# def exponential_decay(x, tau):
#     return (tau * np.power(beta, x))
#
# # this function just returns the age parameter
# # of the product model
# def age_tau(a, tau):
#     return np.power(a, tau)
#     #return np.power(a, tau)
#
# def deg_tau(d, tau):
#     return np.power(d, 1-tau)
#
#
# def LogLikelihood(age, deg, params):
#     """ the negative log-Likelihood-Function """
#     # print(params)
#     # the first variable is independent of the parameter tau so can just be summed up
#     lnl = 0
#     for idx in range(len(age)):
#         try:
#             lnl += ( math.log(deg[idx]) + (params*math.log(age[idx])))
#         except:
#             #print(deg[idx], age[idx])
#             continue
#     #lnl =  np.sum(-np.log((params*age) + np.log(deg)))
#     return lnl
#
# def regressLL(params, age, deg):
#     # Resave the initial parameter guesses
#     b0 = params[0]
#     b1 = params[1]
#     sd = params[2]
#
#     # Calculate the predicted values from the initial parameter guesses
#     yPred = math.pow(age, params) * deg
#     # Calculate the negative log-likelihood as the negative sum of the log of a normal
#     # PDF where the observed values are normally distributed around the mean (yPred)
#     # with a standard deviation of sd
#     logLik = -np.sum(stats.norm.logpdf(yObs, loc=yPred, scale=sd) )
#
#     # Tell the function to return the NLL (this is what will be minimized)
#     return(logLik)
#
# # Make a list of initial parameter guesses (b0, b1, sd)
# initParams = [1, 1, 1]
#
# # Run the minimizer
# results = minimize(regressLL, initParams, method='nelder-mead')

if __name__ == "__main__":

    time_end_degree = pickle.load(
        open('src_age_tgt_deg_v2_t07.pickle', 'rb'))

    # lambdas = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    lambdas = list(np.arange(0, 0.4, 0.03))
    tgt_indegree = []
    node_age = []
    sum_lambdas = [0 for i in range(len(lambdas))]
    for l in range(len(lambdas)):
        count = 0
        for nd, na in time_end_degree:
            na = na/(60*24)
            try:
                # print(math.pow(na, lambdas[l]) * nd)
                sum_lambdas[l] -= math.log(math.pow(na, lambdas[l]) * nd)
                tgt_indegree.append(nd)
                node_age.append(na)
                count += 1
                # if count > 1000000:
                #     break
            except:
                continue

    # node_age = node_age[:1000]
    # tgt_degree = tgt_indegree[:1000]

    # orig_x = np.array(node_age)
    # orig_y = np.array(node_degree)
    #
    beta_list = list(np.arange(0, 1, 0.05))
    residuals = []
    param_values = []
    lnl_list = []
    #
    # Find the age probability from the data
    # in_values = sorted(set(node_age))
    # list_values = list(node_age)
    # in_hist = [list_values.count(x) for x in list_values]
    # in_hist = [in_hist[i] / len(node_age) for i in range(len(in_hist))]
    #
    #
    # # Find the node degree probability from the data
    # in_values = sorted(set(tgt_indegree))
    # list_values = list(tgt_indegree)
    # in_hist_deg = [list_values.count(x) for x in list_values]
    # in_hist_deg = [in_hist_deg[i] / len(tgt_indegree) for i in range(len(in_hist_deg))]

    # print(len(in_hist))
    # print(len(in_hist_deg))
    #
    # for i in beta_list:
    #     beta = i
    #     #popt, pcov = curve_fit(exponential_decay, orig_x, orig_y)
    #
    #     # plt.plot(orig_x, orig_y, 'o')
    #     # yfit = exponential_decay(orig_x, popt[0])
    #     # plt.plot(orig_x, yfit, 'o', color='green')
    #
    #
    #     # plt.legend(['data', 'fit', 'CI 95%'], loc='best', fontsize=20)
    #     # plt.show()
    #
    #     #print(beta)
    #     param_values.append(beta)
    #     #residuals_exp = 0
    #     #for idx in range(len(orig_y)):
    #     #    residuals_exp += orig_y[idx] - exponential_decay(orig_x[idx], beta)
    #
    #     #residuals.append(residuals_exp)
    #     #print(residuals)
    #
    #     lnl_value = LogLikelihood(in_hist, in_hist_deg, beta)
    #     lnl_list.append(lnl_value)
    # #     #print(lnl_list)
    # #
    # #     # plt.figure()
    # #     # plt.plot(in_values, in_hist, 'ro')
    # #     # plt.ylim([-0.1, 0.2])
    # #     # plt.xlabel('Degree (log)')
    # #     # plt.ylabel('Number of vertices(log)')
    # #     # plt.title('Users network')
    # #     # plt.show()
    # # # plt.plot(beta_list, residuals, 'o')
    # # # plt.xlabel('Beta coefficients')
    # # # plt.ylabel('Residuals (LS)')
    # # # plt.show()
    # # #
    # plt.close()
    # plt.plot(param_values, lnl_list, 'o')
    # plt.ylabel('Log Likelihood')
    # plt.xlabel('Taus')
    # plt.show()
    # #
    # #
    # #
    # # # plt.figure()
    # # # plt.plot(in_values, in_hist, 'ro')
    # # # plt.ylim([-0.1, 0.2])
    # # # plt.xlabel('Degree (log)')
    # # # plt.ylabel('Number of vertices(log)')
    # # # plt.title('Users network')
    # # # plt.show()

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

    mpl.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',  # i need upright \micro symbols, but you need...
        r'\sisetup{detect-all}',  # ...this to force siunitx to actually use your fonts
        r'\usepackage{helvet}',  # set the normal font here
        r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
        r'\sansmath'  # <- tricky! -- gotta actually tell tex to use!
    ]

    plt.plot(lambdas, sum_lambdas, marker='o', linestyle='--', color='r', lw=5, markersize=20)
    plt.xlabel(r"$ \lambda $", fontsize=55)
    plt.ylabel('Log-likelihood', fontsize=55)
    plt.tick_params('y', labelsize=35)
    plt.tick_params('x', labelsize=35)
    plt.grid(True)
    plt.show()
