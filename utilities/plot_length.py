import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle

length_series = pickle.load(open('..//data//network_stats//length_series.pickle', 'rb'))
length_new = []
for i in range(len(length_series)):
    length_new.append(length_series[i])

plt.figure(figsize=(12, 8))
hfont = {'fontname': 'Arial'}
n, bins, patches = plt.hist(length_new, 50, facecolor='b')
plt.yscale('log', nonposy='clip', basey=2)
plt.xlabel('Length of the time series', size=40, **hfont)
plt.ylabel('Frequency', size=40, **hfont)
# plt.title('Histogram of')
plt.xlim([0, 700])
plt.ylim([0, 2**(12)])
plt.grid(True)
plt.xticks(size=25)
plt.yticks(size=25)
# file_save = dir_save + '/' + 'count_motif_' + str(m) + '.png'
# plt.savefig(file_save)
plt.subplots_adjust(left=0.16, bottom=0.16)
plt.show()

plt.close()

