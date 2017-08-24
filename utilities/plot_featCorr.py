from pylab import *
from math import *
from numpy.random import *
import pickle
from sklearn.metrics import jaccard_similarity_score


def plot_line(x, y1, y2, y3, y4, y5, l1, l2, l3, l4, l5, ye1, ye2, ye3, ye4, ye5, title, ):
    fig = plt.figure(1, figsize=(24, 16))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)
    plt.xticks(arange(len(x)), x, size=35)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)

    x = np.array(range(len(x))) + 1
    plt.plot(x, y1,  marker='o', linestyle='--', label=l1, markersize=5, linewidth=4)
    plt.plot(x, y2, marker='^', linestyle='--', label=l2, markersize=5, linewidth=4)
    plt.plot(x, y3, marker='p', linestyle='--', label=l3, markersize=5, linewidth=4)
    plt.plot(x, y4, marker='*', linestyle='--', label=l4, markersize=5, linewidth=4)
    plt.plot(x, y5, marker='3', linestyle='--', label=l5, markersize=5, linewidth=4)

    plt.xlim(0, len(x) + 1)
    # plt.tight_layout()  # showing xticks (as xticks have long names)
    # ax.grid()

    plt.title(title, color='#000000', weight="bold", size=30)
    plt.xlabel('Intervals before steep region', size=50)
    plt.ylabel('Mean Jaccard Similarity ', size=50)

    ax.legend(loc='top left', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    plt.ylim([0, 0.18])
    # plt.xlim([0, 6])
    plt.savefig('plots/feat_corr/png/steep/'+ title[9:]+ '.png')
    plt.savefig('plots/feat_corr/pdf/steep/'+ title[9:]+ '.pdf')
    plt.close()


if __name__ == "__main__":
    # Load the feature values
    measures = {}
    measures['Pagerank'] = pickle.load(open('11_13/pr/pr_top/v1/inhib/pr_top_nodes.pickle', 'rb'))
    measures['Nodal Degree'] = pickle.load(open('11_13/deg/deg_top/v1/steep/deg_top_nodes.pickle', 'rb'))
    measures['Betweenness'] = pickle.load(open('11_13/bw/bw_top/v1/steep/bw_top_nodes.pickle', 'rb'))
    measures['Entropy'] = pickle.load(open('11_13/entropy/entropy_top/v1/inhib/ent_top_nodes.pickle', 'rb'))
    measures['Clustering coefficient'] = pickle.load(open('11_13/cc/cc_top/v1/inhib/cc_top_nodes.pickle', 'rb'))
    measures['Alpha Cent.'] = pickle.load(open('11_13/katz/katz_top/v1/inhib/katz_top_nodes.pickle', 'rb'))

    number_intervals = 21
    # Find the top rank jaccard similarity between the nodes
    for m in measures:
        print('Computing for ', m)
        val_list = [{} for _ in range(2000)]
        for interval in range(1, number_intervals):
            # Top nodes in the measure considered
            int_reverse = number_intervals - interval
            measure_list = measures[m][int_reverse]

            for mid in measure_list:
                top_nodes = measure_list[mid][:20]
                # val_list[interval-1][mid] = []
                if top_nodes == []:
                    continue
                val_list[interval-1][mid] = top_nodes

        jacSim_vals = {}
        for m_other in measures:
            print("Feature correlation wth: ", m_other)
            if m == m_other:
                continue
            jacSim_vals[m_other] = {}
            val_list_others = [{} for _ in range(2000)]
            for interval in range(1, number_intervals):
                # Top nodes in the other measure considered
                int_reverse = number_intervals - interval
                measure_list = measures[m_other][int_reverse]

                for mid in measure_list:
                    top_nodes = measure_list[mid][:20]
                    # print(top_nodes)
                    # val_list_others [interval-1][mid] = []
                    if top_nodes == []:
                        continue
                    val_list_others[interval - 1][mid] = top_nodes

            for i in range(20):
                jacSim_vals[m_other][i] = []
                for mid in val_list[i]:
                    if mid not in val_list_others[i]:
                        continue

                    jaccard_sim = jaccard_similarity_score(val_list[i][mid], val_list_others[i][mid])
                    jacSim_vals[m_other][i].append(jaccard_sim)
                    # if jaccard_sim == 0.0:
                    #     print(i, m, val_list[i][mid])
                    #     print(i, m_other, val_list_others[i][mid])
            # print(jacSim_vals[m_other])


        labels = []
        y_vals = [[] for _ in range(len(jacSim_vals))]
        y_err = [[] for _ in range(len(jacSim_vals))]
        cnt_m = 0
        for mval in jacSim_vals:
            labels.append(mval)
            for i in range(20):#jacSim_vals[mval]:
                y_vals[cnt_m].append(np.mean(np.array(jacSim_vals[mval][i])))
                y_err[cnt_m].append(np.std(np.array(jacSim_vals[mval][i])))

            cnt_m += 1

        x = range(20)
        # print(labels)
        plot_line(x, y_vals[0], y_vals[1], y_vals[2], y_vals[3], y_vals[4], labels[0], labels[1],
                  labels[2], labels[3], labels[4], y_err[0], y_err[1], y_err[2], y_err[3], y_err[4],
                  "Feature: " + m)







