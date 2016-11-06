# N.B. : the intervals stored in the pickle files are reverse in order.
# reverse the order in this program to get the original order of the intervals.

# this is the general utility plot program
__author__ = 'ssarka18'

import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import seaborn

import datetime
import pickle

number_intervals = 20
path = '11_05/entropy/entropy_velocity/v1/inhib'
cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
data_to_plot = []
titles = []

for filename in os.listdir(path):
    print("Reading file...", filename)
    full_path = path + '/' + filename
    global_list = pickle.load(open(full_path, 'rb'))
    print(global_list)
    for interval in range(1, number_intervals):
        int_reverse = number_intervals - interval
        measure_list = global_list[int_reverse]
        cent_values = []
        cent = []
        val_list = []
        sum_cent = 0
        cnt = 0
        for mid in measure_list:
            measure, num_edges = measure_list[mid]
            print(measure)
            if measure == []:
                continue
            if measure == NAN:
                continue
            # print(idx)
            # for idx_2 in range(len(idx)):
            val_list.append(measure)
            # val_list.append(idx)
        data_to_plot.append(val_list)
        titles.append('I' + str(interval))

#data_to_plot = data_to_plot.reverse()


fig = plt.figure(1, figsize=(10, 6))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_to_plot, patch_artist=True)

for box in bp['boxes']:
    # change outline color
    box.set(color='#7570b3', linewidth=2)
    # change fill color
    box.set(facecolor='#1b9e77')

## change color and linewidth of the whiskers
for whisker in bp['whiskers']:
    whisker.set(color='#7570b3', linewidth=2)

## change color and linewidth of the caps
for cap in bp['caps']:
    cap.set(color='#7570b3', linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#b2df8a', linewidth=2)

# change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#e7298a', alpha=0.5)

third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
third_quartile = max(third_quartile)

first_quartile = [item.get_ydata()[1] for item in bp['whiskers']]
first_quartile = max(first_quartile)

ax.set_title('Entropy', fontsize=55)
#ax.set_title(r'\textbf{Shortest path - Newly appeared nodes by interval}', fontsize=55)
ax.set_xlabel('Intervals before steep region', fontsize=50)
plt.ylim([-third_quartile - 0.5*math.pow(10, int(math.log10(third_quartile))),
          third_quartile + math.pow(10, int(math.log10(third_quartile)))])
# plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
# plt.ylim([-0.01, 0.1])
plt.tick_params('y', labelsize=50)
# ax.set_xticklabels(titles)
plt.show()
