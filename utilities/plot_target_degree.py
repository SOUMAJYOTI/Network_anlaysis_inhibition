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
from pylab import *


if __name__ == "__main__":

    time_end_degree = pickle.load(
        open('edges_degree_t07.pickle', 'rb'))

    tgt_indegree = []
    node_age = []
    node_count = []
    sum_nodes = 0
    for age in time_end_degree:
        sum_nodes += time_end_degree[age]
    for age in time_end_degree:
        # print(age, time_end_degree[age])
        node_age.append(age)
        node_count.append(time_end_degree[age]/sum_nodes)
    # for nd, na in time_end_degree:
    #     # if nd > 1000:
    #     #     continue
    #     if na<=0:
    #         continue
    #     tgt_indegree.append(nd)
    #     node_age.append(na)

    # node_age = node_age[:1000]
    # tgt_degree = tgt_indegree[:1000]

    # orig_x = np.array(node_age)
    # orig_y = np.array(node_degree)
    #

    #     #print(lnl_list)
    #
    # Create an Axes object.
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
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    # Plot the data.
    ax.scatter(node_age, node_count, color="red", marker="o")
    # Add a title.
    plt.xlabel(r" Target node in-degree ", fontsize=55)
    plt.ylabel(r"$ p_e(d) $", fontsize=55)
    # Add some axis labels.
    plt.tick_params('y', labelsize=35)
    plt.tick_params('x', labelsize=35)
    plt.ylim([-0.001, 0.05])
    plt.grid(True)
    plt.show()