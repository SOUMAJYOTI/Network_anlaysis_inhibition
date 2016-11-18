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
import datetime
import pickle
from pylab import *
from mpl_toolkits.mplot3d import Axes3D


number_intervals = 21
path = '../data/centralities/11_14/deg/deg_values/v1/inhib/deg.pickle'
path_2 = '../data/centralities/11_14/entropy/entropy_values/inhib/entropy.pickle'

cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
data_to_plot = []
titles = []
#
# filename = os.listdir(path)
# filename_2 = os.listdir(path_2)
# print("Reading file...", filename)
# full_path = fil
# full_path_2 = filename_2
global_list = pickle.load(open(path, 'rb'))
global_list_2 = pickle.load(open(path_2, 'rb'))
# print(global_list)
x=  []
y = []
z = []
colors = ['silver', 'red', 'green', 'blue']
for interval in range(1, number_intervals):
    int_reverse = number_intervals - interval
    measure_list = global_list[int_reverse]
    measure_list_2 = global_list_2[int_reverse]
    cent_values = []
    cent = []
    val_list = []
    sum_cent = 0
    cnt = 0
    print(len(measure_list))
    for mid in measure_list:
        try:
            measure, num_edges = measure_list[mid]
            measure_2, num_edges_2 = measure_list_2[mid]
        except:
            continue
        # print(measure)
        if measure == [] or measure_2 == []:
            continue
        if measure == NAN or measure_2 == NAN:
            continue
        # print(idx)
        # for idx_2 in range(len(idx)):
        x.append(int_reverse)
        y.append(measure)
        z.append(measure_2)
        val_list.append(measure)
    print(len(x))
        # val_list.append(idx)
    data_to_plot.append(val_list)
    titles.append('I' + str(interval))

print(len(x))


fig = plt.figure()
ax = plt.axes(projection='3d')
rgb = np.ones((z.shape[0], z.shape[1], 3))


# Get lighting object for shading surface plots.
from matplotlib.colors import LightSource

# Get colormaps to use with lighting object.
from matplotlib import cm

# Create an instance of a LightSource and use it to illuminate the surface.
light = LightSource(90, 45)

green = np.array([0,1.0,0])
green_surface = light.shade_rgb(rgb * green, z)
ax = Axes3D(plt.figure())
ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, antialiased=False,
                facecolors=green_surface)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
#data_to_plot = data_to_plot.reverse()
#
#
# fig = plt.figure(1, figsize=(10, 8))
#
# # Create an axes instance
# ax = fig.add_subplot(111)
#
# # Create the boxplot
# bp = ax.boxplot(data_to_plot, patch_artist=True)
#
# for box in bp['boxes']:
#     # change outline color
#     box.set(color='#0000FF', linewidth=2)
#     # change fill color
#     box.set(facecolor='#FFFFFF')
#
#     ## change color and linewidth of the whiskers
#     # for whisker in bp['whiskers']:
#     #     whisker.set(color='#7570b3', linewidth=2)
#
#     ## change color and linewidth of the caps
#     # for cap in bp['caps']:
#     #     cap.set(color='#7570b3', linewidth=2)
#
#     ## change color and linewidth of the medians
# for median in bp['medians']:
#     median.set(color='#FF0000', linewidth=4)
#
#     ## change the style of fliers and their fill
# for flier in bp['fliers']:
#     flier.set(marker='o', color='#e7298a', alpha=0.5)
#
# third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
# third_quartile = max(third_quartile)
#
# first_quartile = [item.get_ydata()[1] for item in bp['whiskers']]
# first_quartile = max(first_quartile)
#
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
#
# # ax.set_title(r'\textbf{Entropy}', fontsize=25)
# #ax.set_title(r'\textbf{Shortest path - Newly appeared nodes by interval}', fontsize=55)
# ax.set_xlabel(r'\textbf{Intervals before steep region}', fontsize=25)
# ax.set_ylabel(r'\textbf{Degree entropy change}', fontsize=25)
#
# # plt.ylim([-third_quartile - 0.5*math.pow(10, int(math.log10(third_quartile))),
# #           third_quartile + math.pow(10, int(math.log10(third_quartile)))])
# # plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
# # plt.ylim([-0.2, 0.2])
# plt.tick_params('y', labelsize=20)
# plt.tick_params('x', labelsize=20)
# plt.grid(True)
# # ax.set_xticklabels(titles)
# plt.show()
