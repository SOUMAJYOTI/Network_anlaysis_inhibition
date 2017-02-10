import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator as o

import numpy as np
#
dpoints_inhib = np.array([['Model 1-full series', 'Degree entropy', 140.834],
                    ['Model 1-full series', 'Clustering coeff.', 257.166],
                    ['Model 1-full series', 'Pagerank', 168.279],
                    ['Model 1-full series', 'Betweenness', 161.062],
                    ['Model 1-full series', 'Nodal degree', 163.884],
                    ['Model 2-full series', 'Degree entropy', 100.062],
                    ['Model 2-full series', 'Clustering coeff.', 375.577],
                    ['Model 2-full series', 'Pagerank', 140.55],
                    ['Model 2-full series', 'Betweenness', 166.292],
                    ['Model 2-full series', 'Nodal degree', 131.027],
                    ['Model 1-clipped series', 'Degree entropy', 81.182],
                    ['Model 1-clipped series', 'Clustering coeff.', 224.043],
                    ['Model 1-clipped series', 'Pagerank', 145.251],
                    ['Model 1-clipped series', 'Betweenness', 175.769],
                    ['Model 1-clipped series', 'Nodal degree', 146.163],
                    ['Model 2-clipped series', 'Degree entropy', 180.568],
                    ['Model 2-clipped series', 'Clustering coeff.', 314.015],
                    ['Model 2-clipped series', 'Pagerank', 213.198],
                    ['Model 2-clipped series', 'Betweenness', 201.728],
                    ['Model 2-clipped series', 'Nodal degree', 228.667]])

# dpoints = np.array([['Model 1-full series', 'Degree entropy', 33.65],
#                     ['Model 1-full series', 'Clustering coeff.', 400.00],
#                     ['Model 1-full series', 'Pagerank', 52.738],
#                     ['Model 1-full series', 'Betweenness', 32.131],
#                     ['Model 1-full series', 'Nodal degree', 42.99],
#                     ['Model 2-full series', 'Degree entropy', 83.002],
#                     ['Model 2-full series', 'Clustering coeff.', 267.575],
#                     ['Model 2-full series', 'Pagerank', 116.66],
#                     ['Model 2-full series', 'Betweenness', 78.20],
#                     ['Model 2-full series', 'Nodal degree', 92.50]])

hfont = {'fontname': 'Arial'}
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)


def barplot(ax, dpoints):
    '''
    Create a barchart for data across different categories with
    multiple conditions for each category.

    @param ax: The plotting axes from matplotlib.
    @param dpoints: The data set as an (n, 3) numpy array
    '''

    # Aggregate the conditions and the categories according to their
    # mean values
    conditions = [(c, np.mean(dpoints[dpoints[:, 0] == c][:, 2].astype(float)))
                  for c in np.unique(dpoints[:, 0])]
    categories = [(c, np.mean(dpoints[dpoints[:, 1] == c][:, 2].astype(float)))
                  for c in np.unique(dpoints[:, 1])]

    colors = ['#C0C0C0', '#000000', '#00FFFF', '#0000FF']
    hatch_pat=['-', '+', 'x', '\\']
    # sort the conditions, categories and data so that the bars in
    # the plot will be ordered by category and condition
    conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
    categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]

    dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))

    # the space between each set of bars
    space = 0.3
    n = len(conditions)
    width = (1 - space) / (len(conditions))

    # Create a set of bars at each position
    for i, cond in enumerate(conditions):
        indeces = range(1, len(categories) + 1)
        vals = dpoints[dpoints[:, 0] == cond][:, 2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in indeces]
        ax.bar(pos, vals, width=width, label=cond,
               color='#C0C0C0',
               edgecolor='black',
               # fill=False,
               hatch=hatch_pat[i])#cm.Accent(float(i) / n))

    # Set the x-axis tick labels to be equal to the categories
    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=50)

    # Add the axis labels
    ax.set_ylabel("Mean Absolute error", size=40, **hfont)
    # ax.set_xlabel("Structure")

    # Add a legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left', fontsize=20)
plt.tick_params('x', labelsize=25)
plt.tick_params('y', labelsize=25)
plt.grid(True)
plt.rc('legend', **{'fontsize':8})
plt.subplots_adjust(left=0.13, bottom=0.30, top=0.9)

barplot(ax, dpoints_inhib)
# savefig('barchart_3.png')
plt.show()