import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle

length_series = pickle.load(open('..//..//data_files//results_granger//11_23//length_series.pickle', 'rb'))
length_new = []
for i in range(len(length_series)):
    length_new.append(length_series[i])

n, bins, patches = plt.hist(length_new, 50, facecolor='g')
plt.xlabel('Length of the time series', size=30)
plt.ylabel('Frequency', size=30)
# plt.title('Histogram of')
plt.grid(True)
plt.xticks(size=20)
plt.yticks(size=20)
# file_save = dir_save + '/' + 'count_motif_' + str(m) + '.png'
# plt.savefig(file_save)
plt.show()
plt.close()

