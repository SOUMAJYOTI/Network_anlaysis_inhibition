import matplotlib.pyplot as plt
import pickle

growth = pickle.load(open('../data/network_stats/size_cascades.pickle', 'rb'))
time_diff = pickle.load(open('../data/network_stats/time_cascades.pickle', 'rb'))

g_new = []
t_new = []
for idx in range(len(time_diff)):
    # if time_diff[idx] <= 20000:
        g_new.append(growth[idx])
        t_new.append(time_diff[idx])

g_new.extend(g_new)
t_new.extend(t_new)

g_new.extend(g_new)
t_new.extend(t_new)

g_new.extend(g_new)
t_new.extend(t_new)

# in_values = sorted(set(g_new))
# list_degree_values = list(g_new)
#
#
# in_hist = [list_degree_values.count(x)/len(g_new) for x in in_values]
# idx_max = in_hist.index(max(in_hist))
# print(list_degree_values[idx_max], max(in_hist))
# print(list_degree_values[0], in_hist[0])
# plt.figure()
# plt.plot(in_values, in_hist, '.')
# #plt.xlim([1, 2**18])
# plt.xlabel('Degree (log)')
# plt.ylabel('Number of vertices(log)')
# plt.title('Users network')
# tick_locs = [0, 500, 1000, 1500, 2000, 4000, 6000]
# plt.xticks(tick_locs)
# plt.xlim([0, 6000])
# plt.grid(True)
# plt.show()
#
# exit()
# plt.rc('text', usetex=True)
# plt.rc('axes')
# plt.rc('font', family='arial')
# plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

hfont = {'fontname': 'Arial'}

plt.close()
plt.figure(figsize=(12, 8))
n, bins, patches = plt.hist(g_new, bins=20, facecolor='b')
plt.xlabel('Cascade size', size=40, **hfont)
plt.ylabel('Frequency', size=40, **hfont)
plt.title('')
plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)
plt.subplots_adjust(left=0.16, bottom=0.16)
plt.grid(True)
plt.show()
# plt.savefig('Cascade_figures/growth_1.png')

plt.close()
plt.figure(figsize=(12, 8))
n, bins, patches = plt.hist(t_new, bins=20, facecolor='b')
plt.xlabel('Cascade lifetime', size=40, **hfont)
plt.ylabel('Frequency', size=40, **hfont)
plt.title('')
plt.tick_params(axis='x', labelsize=25)
plt.tick_params(axis='y', labelsize=25)
plt.grid(True)
plt.subplots_adjust(left=0.16, bottom=0.16)
plt.show()
# plt.savefig('Cascade_figures/time_diff_1.png')


