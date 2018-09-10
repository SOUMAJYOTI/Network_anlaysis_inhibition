import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *

import datetime
import pickle
import numpy as np
import matplotlib.image as image
from operator import itemgetter
import math

path = 'edges_list.pickle'
cnt_rec = 0

data_list = pickle.load(open(path, 'rb'))
data_plot = [[] for _ in range(20)]
for d in range(len(data_list)):
    cid = data_list[d]
    num_int = len(cid)-1
    print(num_int)
    count_int = 0
    temp = []
    while(True):
        # print(count_int)
        temp.append(cid[num_int])
        data_plot[count_int].append(cid[num_int])
        count_int  += 1
        num_int -= 1

        if count_int >= 20 or num_int == -1:
            break
    # data_plot.append(temp)

data_plot = list(reversed(data_plot))

print(data_plot)
# for idx in range()
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
#
# mpl.rcParams['text.latex.preamble'] = [
#     r'\usepackage{siunitx}',  # i need upright \micro symbols, but you need...
#     r'\sisetup{detect-all}',  # ...this to force siunitx to actually use your fonts
#     r'\usepackage{helvet}',  # set the normal font here
#     r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#     r'\sansmath'  # <- tricky! -- gotta actually tell tex to use!
# ]

# Create the box_plots
# Create the box_plots
fig = plt.figure(1, figsize=(10, 8))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_plot, patch_artist=True)

for box in bp['boxes']:
    # change outline color
    box.set(color='#0000FF', linewidth=2)
    # change fill color
    box.set(facecolor='#FFFFFF')

## change color and linewidth of the whiskers
# for whisker in bp['whiskers']:
#     whisker.set(color='#7570b3', linewidth=2)

## change color and linewidth of the caps
# for cap in bp['caps']:
#     cap.set(color='#7570b3', linewidth=2)

## change color and linewidth of the medians
for median in bp['medians']:
    median.set(color='#FF0000', linewidth=4)

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#e7298a', alpha=0.5)

#ax.set_ylim([0, 0.05])

third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
third_quartile = max(third_quartile)

dir_save = '../plots/motif_count_plots/12_08/v2/inhib'
if not os.path.exists(dir_save):
    os.makedirs(dir_save)
# file_save = dir_save + '/' + 'mc_inhib_' + str(m) + '.png'
# plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
# plt.ylim([0, limits_y_steep[m]])
# major_ticks = np.arange(0, 21, 5)
# plt.xticks(major_ticks)
plt.tick_params('y', labelsize=20)
plt.tick_params('x', labelsize=20)
plt.xlabel('Network subsequences leading to inhibition', fontsize=25)
plt.ylabel('Edge cardinality', fontsize=25)
# limits_y_steep[m] = third_quartile + math.pow(10, int(math.log10(third_quartile)))
# plt.ylim([0, third_quartile + 4*math.pow(10, int(math.log10(third_quartile)))])
plt.grid(True)
plt.show()
plt.close()
    # plt.show()