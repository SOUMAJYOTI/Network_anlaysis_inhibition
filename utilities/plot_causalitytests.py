import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import math

if __name__ == "__main__":
    statistic_values = pickle.load(open('..//..//data_files//results_granger//11_15//statistics.pickle', 'rb'))
    crit_values = pickle.load(
        open('..//..//data_files//results_granger//11_15//critical.pickle', 'rb'))
    p_values = pickle.load(
        open('..//..//data_files//results_granger//11_15//p_values.pickle',  'rb'))
    cause_count = pickle.load(
        open('..//..//data_files//results_granger//11_15//cause_count.pickle', 'rb'))

    map_cent = ['cc', 'entropy', 'nbr_deg', 'bw', 'pr']
    data_to_plot = []
    titles = []
    for idx in range(len(map_cent)):
        sub = map_cent[idx]
        temp = []
        # num_feat = len(sub.split('+'))
        # if num_feat != 1:
        #     continue
        for idx in range(len(p_values[sub])):
            if statistic_values[sub][idx] <= 0:
                continue
            temp.append(p_values[sub][idx])
        data_to_plot.append(temp)
        titles.append(sub)
    fig = plt.figure(1, figsize=(12, 8))

    # Create an axes instance
    ax = fig.add_subplot(111)

    # Create the boxplot
    bp = ax.boxplot(data_to_plot, patch_artist=True)

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

    third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
    third_quartile = max(third_quartile)

    first_quartile = [item.get_ydata()[1] for item in bp['whiskers']]
    first_quartile = max(first_quartile)

    # ax.set_title('Entropy', fontsize=55)
    # ax.set_title(r'\textbf{Shortest path - Newly appeared nodes by interval}', fontsize=55)
    # ax.set_xlabel(, fontsize=25)
    # plt.ylim([-third_quartile - 0.5 * math.pow(10, int(math.log10(third_quartile))),
    #           third_quartile + math.pow(10, int(math.log10(third_quartile)))])
    # plt.ylim([0, third_quartile + math.pow(10, int(math.log10(third_quartile)))])
    # plt.ylim([0, 3])
    plt.tick_params('y', labelsize=25)
    ax.set_xticklabels(titles, size=25)
    plt.grid(True)
    plt.show()

    # for sub in statistic_values:
    #     # print(len(statistic_values[sub]))
    #     plt.figure()
    #     n, bins, patches = plt.hist(statistic_values[sub], bins=20, facecolor='g')
    #     plt.xlabel('Statistic values - Wald test')
    #     plt.ylabel('Frequency')
    #     plt.title(sub)
    #     plt.grid(True)
    #     plt.show()
    #     # plt.savefig('F://Inhibition//VAR_causality//plots//granger_plots//statistics_' + sub + '.png')
    #     plt.close()
    #
    # cause_percentage = []
    # idx = 0
    # cause_idx = []
    # titles = []
    # for sub in statistic_values:
    #     num_feat = len(sub.split('+'))
    #     if num_feat != 1:
    #             continue
    #
    #     cause_percentage.append(cause_count[sub]/len(statistic_values[sub])*100 )
    #     cause_idx.append(idx)
    #     print(sub, cause_count[sub]/len(statistic_values[sub])*100 )
    #     titles.append(sub)
    #     idx += 1
    # print(len(titles))
    # fig = plt.figure()
    # ax = fig.add_subplot(1, 1, 1)  # one row, one column, first plot
    # # Plot the data.
    # ax.scatter(cause_idx, cause_percentage, color="red", marker="o")
    # plt.xticks(range(len(titles)), titles, rotation=75, fontsize=20)
    # # ax.set_xticklabels(titles, size=20)
    # # ax.scatter(p_list_null, f_list_null, color="blue", marker="o")
    # # Add a title.
    # # ax.set_title("P_values vs Chi-square wald")
    # # Add some axis labels.
    # # ax.set_xlabel("P-values")
    # ax.set_ylabel("Percentage of cascades where feat GCaused T", fontsize=20)
    # # for tick in ax.get_xticklabels():
    # #     tick.set_rotation(45)
    # plt.subplots_adjust(bottom=0.15)
    # # plt.ylim([0, 10])
    # plt.show()

    # for sub in statistic_values:
    #     print(sub)
    #     # num_feat = len(sub.split('+'))
    #     # if num_feat != 1:
    #     #     continue
    #     fig = plt.figure()
    #     ax = fig.add_subplot(1, 1, 1)  # one row, one column, first plot
    #     # Plot the data.
    #     ax.scatter(p_values[sub], statistic_values[sub], color="red", marker="o")
    #     # ax.scatter(p_list_null, f_list_null, color="blue", marker="o")
    #     # Add a title.
    #     ax.set_title("P_values vs Chi-square wald")
    #     # Add some axis labels.
    #     ax.set_xlabel("P-values")
    #     ax.set_ylabel("Chi-square wald")
    #     plt.show()
    #     # plt.ylim([0, 10])
    #     plt.show()
    #     # plt.savefig('F://Inhibition//VAR_causality//plots//P_F_entropy_wo_difference_Y.png')
