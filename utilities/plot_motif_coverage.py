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
import numpy as np
import matplotlib.image as image
from operator import itemgetter
import math

path = 'motifs_coveragt_min_4.pickle'
cnt_rec = 0

mdict = pickle.load(open(path, 'rb'))
mdict = sorted(mdict.items())
motifs_values = []
for m, val in mdict:
    val = list(filter((0.0).__ne__, val))
    val = [x / 2 for x in val]
    if val == []:
        continue
    motifs_values.append(val)

for x in motifs_values:
    print(max(x))

# motifs_values = sorted(motifs_values, key=lambda k: median(k))
# Create the box_plots
fig = plt.figure(1, figsize=(10, 8))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(motifs_values, patch_artist=True)

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

## change the style of fliers and their fill
for flier in bp['fliers']:
    flier.set(marker='o', color='#e7298a', alpha=0.5)

#ax.set_ylim([0, 0.05])

third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
third_quartile = max(third_quartile)

dir_save = 'motif_weights_plots/generative/v4/5/steep'
if not os.path.exists(dir_save):
    os.makedirs(dir_save)
file_save = dir_save + '/' + 'mw_steep_' + str(m) + '.png'
# plt.ylim([0, third_quartile + 0.2*math.pow(10, int(math.log10(third_quartile)))])
# major_ticks = np.arange(0, 21)
# plt.xticks(major_ticks)
plt.tick_params('y', labelsize=15)
plt.tick_params('x', labelsize=15)
plt.xlabel(r'\textbf{Edge Densities (-> Increasing)}', fontsize=15)
plt.ylabel(r'\textbf{Percentage cover}', fontsize=15)
# limits_y_steep[m] = third_quartile + 0.2*math.pow(10, int(math.log10(third_quartile)))
# plt.ylim([0, y])
# plt.savefig(file_save)
plt.show()
plt.close()
