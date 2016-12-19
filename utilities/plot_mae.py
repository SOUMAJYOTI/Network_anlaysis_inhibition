import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle


mae_model = pickle.load(open('..//..//data_files//results_granger//11_23//mae_OLS.pickle', 'rb'))

# map_titles ={'pr': 'clustering coefficient', 'entropy': 'nodal degree', 'cc': 'Pagerank', 'nbr_deg': 'entropy', 'bw': 'betweeneness'}
map_titles ={'pr': 'Pagerank', 'entropy': 'Nodal Degree', 'cc': 'Clustering coefficient', 'nbr_deg': 'Degree entropy', 'bw': 'betweeneness'}

titles = []
mae_new = []
for ms in mae_model:
    temp = []
    for idx in range(len(mae_model[ms])):
        temp.append(math.sqrt(mae_model[ms][idx]))
    mae_new.append(np.mean(temp))
    print(ms, np.mean(temp))
    titles.append(map_titles[ms])

width = 0.35  # the width of the bars
# print(data_mean)
ind = np.arange(len(mae_new))  # the x locations for the groups
print(mae_new)
fig = plt.figure()
ax = fig.add_subplot(111)
## the bars
rects1 = ax.bar(ind, mae_new, width,
                color='#C0C0C0',
                error_kw=dict(elinewidth=2,ecolor='black'))

# axes and labels
ax.set_xlim(-width,len(ind)+width)
# ax.set_ylim(0,45)
ax.set_ylabel('Mean Absolute error', size=30)
# ax.set_title('Scores by group and gender')
xTickMarks = titles
ax.set_xticks(ind)
xtickNames = ax.set_xticklabels(xTickMarks)
plt.setp(xtickNames, rotation=25, fontsize=5)
plt.grid(True)
plt.xticks(size=20)
plt.yticks(size=20)

## add a legend
# ax.legend( (rects1[0], ('Men', 'Women') )

plt.show()

# # Create the box_plots
# fig = plt.figure(1, figsize=(10, 8))
#
# # Create an axes instance
# ax = fig.add_subplot(111)
#
# # Create the boxplot
# bp = ax.boxplot(mae_new, patch_artist=True)
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
#     for median in bp['medians']:
#         median.set(color='#FF0000', linewidth=4)
#
#     ## change the style of fliers and their fill
#     for flier in bp['fliers']:
#         flier.set(marker='o', color='#e7298a', alpha=0.5)
#
# ax.set_title('Absolute errors ', fontsize=30)
# # ax.set_xlabel('Features', fontsize=20)
# #ax.set_ylim([0, 0.05])
# ax.set_xticklabels(titles, size=15)
#
# # for tick in ax.get_xticklabels():
# #     tick.set_rotation(45)
# for tick in ax.yaxis.get_major_ticks():
#     tick.label.set_fontsize(30)
#
# third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
# third_quartile = max(third_quartile)
#
# # # dir_save = 'motif_count_plots/5/inhib'
# # if not os.path.exists(dir_save):
# #     os.makedirs(dir_save)
# # file_save = dir_save + '/' + 'count_motif_' + str(m) + '.png'
# plt.ylim([0, third_quartile + 4*math.pow(10, int(math.log10(third_quartile)))])
# plt.subplots_adjust(bottom=0.15)
# # plt.savefig(file_save)
# plt.grid(True)
# plt.show()
# plt.close()
#
