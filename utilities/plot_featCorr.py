from pylab import *
from math import *
from numpy.random import *
import pickle
from sklearn.metrics import jaccard_similarity_score


def plot_line(x, y1, y2,  l1, l2, title, ):
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.subplot(111)  # row x col x position (here 1 x 1 x 1)

    # max_val = np.max(y1)

    # ax.imshow(subim, aspect='auto', extent=(4.8, 7.3, max_val+5.8, max_val + 9), zorder=-1)
    # ax.add_patch(
    #     patches.Rectangle(
    #         (5.2, 0.77),
    #         2.5,
    #         0.28,
    #         fill=False,  # remove background
    #         linewidth=3
    #     )
    # )
    plt.xticks(arange(len(x)), x, size=30)  # rotate x-axis labels to 75 degree
    plt.yticks(size=30)

    x = np.array(range(len(x)-2)) + 1
    # print(x.shape, len(y1))
    ax.plot(x, y1,  marker='o', linestyle='-', label=l1, linewidth=3)
    ax.plot(x, y2, marker='o', linestyle='-', label=l2, linewidth=3)
    # ax.plot(x, y3, marker='o', linestyle='-', label=l3, linewidth=3)
    # ax.plot(x, y4, marker='o', linestyle='-', label=l4, linewidth=3)
    # ax.plot(x, y5, marker='o', linestyle='-', label=l5, linewidth=3)

    # plt.xlim(0, len(var) + 1)
    plt.tight_layout()  # showing xticks (as xticks have long names)
    ax.grid()

    plt.title(title, color='#000000', weight="bold", size=30)
    plt.ylabel('Intervals leading to inhibition', size=40)
    plt.xlabel('Jaccard Similarity (Normalized)', size=40)

    ax.legend(loc='lower right', fancybox=True, shadow=True, fontsize=25)
    plt.grid(True)
    plt.subplots_adjust(left=0.12, bottom=0.15, top=0.90)

    # plt.ylim([0, 1])
    # plt.xlim([0, 6])
    # plt.savefig('outputs/v4/edges_reg_motifs/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/edges_reg_motifs/Motif_' + m + '.pdf')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.png')
    # plt.savefig('outputs/v4/edges/f1/Motif_' + m + '.pdf')
    # plt.close()
    plt.show()
    plt.close()

if __name__ == "__main__":
    # Load the feature values
    measures = {}
    # measures['Pagerank'] = pickle.load(open('11_13/pr/pr_top/v1/inhib/pr_top_nodes.pickle', 'rb'))
    measures['Nodal Degree'] = pickle.load(open('11_13/deg/deg_top/v1/inhib/deg_top_nodes.pickle', 'rb'))
    measures['Betweenness'] = pickle.load(open('11_13/bw/bw_top/v1/inhib/bw_top_nodes.pickle', 'rb'))
    # measures['Entropy'] = pickle.load(open('11_13/entropy/entropy_top/v1/inhib/ent_top_nodes.pickle', 'rb'))
    # measures['Clustering coefficient'] = pickle.load(open('11_13/cc/cc_top/v1/inhib/cc_top_nodes.pickle', 'rb'))
    # measures['Alpha Cent.'] = pickle.load(open('11_13/katz/katz_top/v1/inhib/katz_top_nodes.pickle', 'rb'))

    number_intervals = 20
    # Find the top rank jaccard similarity between the nodes
    for m in measures:
        val_list = [[] for _ in range(20)]
        for interval in range(1, number_intervals):
            # Top nodes in the measure considered
            int_reverse = number_intervals - interval
            measure_list = measures[m][int_reverse]

            for mid in measure_list:
                top_nodes = measure_list[mid]
                if top_nodes == []:
                    continue
                if top_nodes == NAN:
                    continue
                val_list[interval-1].append(top_nodes)

        jacSim_vals = {}
        for m_other in measures:
            if m == m_other:
                continue
            jacSim_vals[m_other] = {}
            val_list_others = [[] for _ in range(20)]
            for interval in range(1, number_intervals):
                # Top nodes in the other measure considered
                int_reverse = number_intervals - interval
                measure_list = measures[m_other][int_reverse]

                for mid in measure_list:
                    top_nodes = measure_list[mid]
                    if top_nodes == []:
                        continue
                    if top_nodes == NAN:
                        continue
                    val_list_others[interval - 1].append(top_nodes)

                jaccard_sim = jaccard_similarity_score(val_list, val_list_others)
                jacSim_vals[m_other][interval-1] = jaccard_sim

        labels = []
        y_vals = [[] for _ in range(len(jacSim_vals))]
        cnt_m = 0
        for mval in jacSim_vals:
            labels.append(mval)
            for i in jacSim_vals[mval]:
                y_vals[cnt_m].append(jacSim_vals[mval][i])

            cnt_m += 1

        print(m)
        x = range(20)
        exit()
        # plot_line(x, y[0], y[1], labels[0], labels[1], m)








