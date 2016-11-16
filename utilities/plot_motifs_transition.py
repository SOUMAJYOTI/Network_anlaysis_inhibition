# Plot the transitions of motif patterns form 4 sized motifs to 5 sized motifs in
# next interval. Observe which transitions occur in which of the intervals.

import pandas as pd
from pylab import *
import os
import math
import statistics as st
import statistics
from math import  *
import datetime
import pickle
import math

path = '../data/motifs/motifs_transition/10_05/5/08/inhib'
number_intervals = 21
cnt_rec = 0

cnt_1 = 0
cnt_2 = 0
motif_trans_interval = {}
titles = []

for filename in os.listdir(path):
    print("Reading file...", filename)
    full_path = path + '/' + filename
    motif_transitions = pickle.load(open(full_path, 'rb'))

    if filename[-9:-8] == '_':
        interval = int(filename[-8:-7])
    else:
        interval = int(filename[-9:-7])

    int_reverse = number_intervals - interval
    if interval>= 21 or interval ==0:
        continue

    # print(motif_transitions)
    try:
        for m4 in motif_transitions:
            if m4 == {}:
                continue
            # print(m4)
            if m4 not in motif_trans_interval:
                motif_trans_interval[m4] = {}
            for m5 in motif_transitions[m4]:
                # print(m5)
                if m5 not in motif_trans_interval[m4]:
                    motif_trans_interval[m4][m5] = [[] for i in range(21)]
                motif_trans_interval[m4][m5][int_reverse-20].extend(motif_transitions[m4][m5])
        #titles.append('I' + str(interval))
    except:
        continue

print('Saving plots...')
limits_y_steep  = {'M4_5': {'M5_6': 0.85846774193548381, 'M5_14': 0.53847802786709531, 'M5_18': 1.150498338870432, 'M5_10': 0.55142857142857138, 'M5_12': 1.0567967698519516, 'M5_11': 1.3, 'M5_17': 1.3, 'M5_9': 1.0850318471337579, 'M5_7': 0.80316251452175025, 'M5_19': 1.3, 'M5_13': 1.3, 'M5_15': 1.3, 'M5_16': 1.0599833610648919, 'M5_8': 0.64162895927601804, 'M5_5': 1.3}, 'M4_0': {'M5_6': 0.076783487903463399, 'M5_14': 0.45414203290832678, 'M5_2': 1.0509604644078123, 'M5_18': 0.885147247119078, 'M5_12': 0.44620952677459524, 'M5_16': 1.3, 'M5_11': 0.12494485294117647, 'M5_3': 0.086301547270124734, 'M5_17': 0.88283132530120478, 'M5_19': 1.3, 'M5_13': 0.11266473760357754, 'M5_7': 0.083954812615320701, 'M5_4': 0.095647847808174896, 'M5_20': 1.3, 'M5_9': 0.41346470417352571, 'M5_8': 0.1179998349358313, 'M5_10': 0.40744929797191887, 'M5_1': 0.076593065517352288, 'M5_0': 0.066086280606296149, 'M5_15': 1.3, 'M5_5': 0.074603853405067366}, 'M4_4': {'M5_6': 1.3, 'M5_14': 0.43705766121789114, 'M5_18': 0.76375545851528381, 'M5_10': 0.10982388973966309, 'M5_12': 0.62183292533659729, 'M5_11': 0.41707688236072515, 'M5_17': 0.72710092905405399, 'M5_13': 0.5200362318840579, 'M5_7': 0.46529634344720278, 'M5_9': 0.60523644135720012, 'M5_15': 0.4790279179280188, 'M5_16': 0.96993870112657388, 'M5_8': 0.51213549733988817, 'M5_5': 1.0712601468882876}, 'M4_2': {'M5_6': 0.093919262034932532, 'M5_14': 0.47475855849084958, 'M5_2': 0.45963053746699345, 'M5_18': 1.3, 'M5_12': 0.43211223694466094, 'M5_16': 0.58901771198767849, 'M5_11': 0.095485961123110155, 'M5_3': 0.1069257416024849, 'M5_17': 0.44696508504923904, 'M5_13': 0.42532929703748212, 'M5_7': 0.11206538344026887, 'M5_4': 0.10837541517526038, 'M5_19': 1.3, 'M5_9': 0.090201258189981467, 'M5_8': 0.60755536682563527, 'M5_10': 0.076416089680184637, 'M5_1': 0.46349129699949515, 'M5_15': 0.45297899878919456, 'M5_5': 0.057705843508435731}, 'M4_3': {'M5_6': 0.12967829695970767, 'M5_14': 0.60593241916124096, 'M5_2': 1.2639451728247915, 'M5_18': 1.3, 'M5_12': 0.47428174574164778, 'M5_16': 0.53712807244501937, 'M5_11': 0.10570002568901611, 'M5_3': 0.53831907228101783, 'M5_17': 0.47010961214165259, 'M5_19': 1.1740249609984399, 'M5_13': 0.51622983870967742, 'M5_7': 0.063304187818915003, 'M5_4': 0.60770707831325299, 'M5_20': 0, 'M5_9': 0.099195787986489758, 'M5_8': 0.42374256843818525, 'M5_10': 0.088602444637265632, 'M5_1': 1.0071438766771339, 'M5_15': 0.44550875168301596, 'M5_5': 0.61501680827222316}, 'M4_1': {'M5_6': 0.08195057315419696, 'M5_14': 0.072439099126704459, 'M5_2': 0.064356343026934826, 'M5_18': 0.61747973696283831, 'M5_12': 0.44329191609235874, 'M5_16': 1.3, 'M5_11': 0.049109203734848922, 'M5_3': 0.10206471675037082, 'M5_17': 1.0931180968564147, 'M5_19': 1.3, 'M5_13': 0.46761246813060287, 'M5_7': 0.10353606902244827, 'M5_4': 0.11361798536075288, 'M5_20': 1.3, 'M5_9': 0.41202873861195466, 'M5_8': 0.097869415807560131, 'M5_10': 0.087766860503977001, 'M5_1': 0.092299276395570423, 'M5_15': 0.98764222815870184, 'M5_5': 0.060181500757625662}}

limits_y_inhib = {'M4_1': {'M5_18': 0.82737421624101537, 'M5_19': 0, 'M5_14': 0.080263741304181638, 'M5_9': 0.10219277342800175, 'M5_1': 0.099778521075145579, 'M5_8': 0.10321743320525989, 'M5_5': 0.069870050727979444, 'M5_4': 0.11423743856720847, 'M5_13': 0.12667044809982983, 'M5_11': 0.08741097836803341, 'M5_12': 0.42320790784244677, 'M5_2': 0.066268467824029495, 'M5_16': 0.62011055112989755, 'M5_10': 0.062149101282129929, 'M5_7': 0.10888517561845924, 'M5_3': 0.11209945613989435, 'M5_15': 0.8442017069924046, 'M5_6': 0.090298286700357444, 'M5_17': 1.1672472387425659, 'M5_20': 0}, 'M4_4': {'M5_18': 0.74792576419213974, 'M5_14': 0.095097652491822199, 'M5_9': 0.4847050754458162, 'M5_6': 0.60663987753539994, 'M5_5': 0.82636199696279422, 'M5_12': 0.1238854723258393, 'M5_13': 0.40247063588497367, 'M5_11': 0.12963444488553658, 'M5_19': 0, 'M5_2': 0, 'M5_16': 0.10456140350877194, 'M5_10': 0.084911975106120308, 'M5_7': 0.48661150010608956, 'M5_8': 0.087593984962406016, 'M5_15': 0.43176724464299232, 'M5_17': 0.51958790714414249, 'M5_20': 0}, 'M4_0': {'M5_18': 1.0619392992635572, 'M5_19': 0, 'M5_14': 0.5326904454486765, 'M5_9': 0.07912161517925248, 'M5_1': 0.086501479527837627, 'M5_7': 0.086944872737623635, 'M5_5': 0.052824474412263941, 'M5_4': 0.10681307669825894, 'M5_13': 0.092299969310359942, 'M5_11': 0.41703724291272926, 'M5_12': 0.47453193126442678, 'M5_8': 0.080484430186967368, 'M5_2': 0.43934650507906464, 'M5_16': 1.3, 'M5_10': 0.06367238584921478, 'M5_0': 0.074279350653233281, 'M5_3': 0.091924885829186828, 'M5_15': 0.54574876577070763, 'M5_6': 0.083873629522099549, 'M5_17': 1.05, 'M5_20': 0}, 'M4_5': {'M5_18': 0.81795841209829856, 'M5_14': 0.81390546360193006, 'M5_9': 1.3, 'M5_6': 0.92570093457943936, 'M5_5': 1.3, 'M5_13': 1.3, 'M5_11': 0.84193607980786989, 'M5_12': 1.1330714846818539, 'M5_16': 0.80653915317966329, 'M5_10': 0.50371428571428567, 'M5_7': 1.0515812572608751, 'M5_8': 0.71137770897832819, 'M5_15': 0.86221812218122174, 'M5_17': 1.2598121085594989, 'M5_19': 1.0673144876325089}, 'M4_3': {'M5_18': 0.51340479483633006, 'M5_14': 0.11704943076419452, 'M5_9': 0.080969441038913861, 'M5_1': 0.67185983442763342, 'M5_7': 0.073480104923851902, 'M5_5': 0.4390015744457525, 'M5_12': 0.40286217393899598, 'M5_13': 0.45162211163831345, 'M5_11': 0.059464384014092257, 'M5_19': 0, 'M5_2': 1.3, 'M5_16': 0.68711103253182459, 'M5_10': 0.068452106990606248, 'M5_8': 0.068384378175701727, 'M5_3': 0.40817292766221669, 'M5_15': 0.59794191190613577, 'M5_20': 1.05, 'M5_6': 0.11008817046289493, 'M5_17': 0.40710310965630114, 'M5_4': 0.47844126506024098}, 'M4_2': {'M5_18': 0.98484876429361856, 'M5_14': 0.4071382403888415, 'M5_9': 0.092788982663966568, 'M5_1': 0.80531179349854032, 'M5_7': 0.12697525206232813, 'M5_5': 0.062368670497232395, 'M5_4': 0.11927746068599315, 'M5_13': 0.11228496196606202, 'M5_11': 0.1041036717062635, 'M5_12': 0.43765197665230804, 'M5_8': 0.1038246355103932, 'M5_2': 0.49357332968639867, 'M5_16': 0.93571583986074836, 'M5_10': 0.075962028721864222, 'M5_6': 0.10670413825583062, 'M5_3': 0.12368255725924, 'M5_15': 0.41219023822562978, 'M5_17': 0.49986905281536442, 'M5_19': 0}}

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

limits_y_steep = {}
for m4 in motif_trans_interval:
    for m5 in motif_trans_interval[m4]:
        max_value_interval = []
        for idx in range(len(motif_trans_interval[m4][m5])):
            if len(motif_trans_interval[m4][m5][idx]) == 0:
                continue
            max_value_interval.append(max(motif_trans_interval[m4][m5][idx]))

        try:
            max_val = max(max_value_interval)
        except:
            continue
        if max_val == 0:
            continue
        data_to_plot = []
        for idx in range(len(motif_trans_interval[m4][m5])):
            for idx_v in range(len(motif_trans_interval[m4][m5][idx])):
                motif_trans_interval[m4][m5][idx][idx_v] = motif_trans_interval[m4][m5][idx][idx_v] /max_val
            data_to_plot.append(motif_trans_interval[m4][m5][idx])
            titles.append(str(idx+1))

        # Create the box_plots
        fig = plt.figure(1, figsize=(10, 8))

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

        #ax.set_title('Motif transition:' + str(m4) + '-->' + str(m5))
        ax.set_ylabel('Interval')
        #ax.set_ylim([0, 0.05])
        ax.set_xticklabels(titles)

        third_quartile = [item.get_ydata()[0] for item in bp['whiskers']]
        third_quartile = max(third_quartile)

        dir_save = '../plots/motif_transition_plots/10_05/inhib'
        if not os.path.exists(dir_save):
            os.makedirs(dir_save)
        file_save = dir_save + '/' + 'mt_inhib_' + str(m4) + '_' + str(m5) + '.png'

        # try:
        #     plt.ylim([0, third_quartile + 0.3 * math.pow(10, int(math.log10(third_quartile)))])
        # except:
        #     print(third_quartile)
        # plt.ylim([0, limits_y_steep[m4][m5]])
        # major_ticks = np.arange(0, 21, 5)
        # plt.xticks(major_ticks)
        plt.tick_params('y', labelsize=25)
        plt.tick_params('x', labelsize=25)
        plt.xlabel(r'\textbf{Intervals before inhibition region}', fontsize=25)
        plt.ylabel(r'\textbf{Normalized motif counts}', fontsize=25)

        if m4 not in limits_y_steep:
            limits_y_steep[m4] = {}
        try:
            limits_y_steep[m4][m5] = third_quartile + 0.3 * math.pow(10, int(math.log10(third_quartile)))
        except:
            limits_y_steep[m4][m5] = 0
        try:
            plt.ylim([0, limits_y_inhib[m4][m5]])
        except:
            pass
        plt.grid(True)
        plt.savefig(file_save)
        plt.close()

print(limits_y_steep)

