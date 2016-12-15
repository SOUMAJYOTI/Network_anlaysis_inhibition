__author__ = 'ssarka18'


import numpy as np
import networkx as nx
import pickle
import pandas as pd
import os
import csv
import matplotlib.pyplot as plt

g = nx.Graph()
print("creating diffusion network...")

node_diff_map = {}
cnt_nodes = 0
line_count = 1
print("Reading diffusion file...")
path = 'diffusion_files'
for filename in os.listdir(path):
    print(filename)
    full_path = path + '/' + filename
    edges = pd.read_csv(full_path)
    sources = edges.ix[:,0]
    targets = edges.ix[:,1]

    nodes = {}
    for row in range(len(sources)):
            nodes[0] = sources[row]
            nodes[1] = targets[row]

            line_count += 1
            if nodes[0] not in node_diff_map:
                node_diff_map[nodes[0]] = cnt_nodes
                cnt_nodes += 1
                #string_labels[int(v1)] = nodes[0]

            v1 = node_diff_map[nodes[0]]

            if nodes[1] not in node_diff_map:
                node_diff_map[nodes[1]] = cnt_nodes
                cnt_nodes += 1
                #string_labels[int(v2)] = nodes[1]

            v2 = node_diff_map[nodes[1]]

            if g.has_edge(v1, v2):
                    continue
            else:
                g.add_edge(v1, v2, weight = 1)

            line_count += 1
            print("Line_count: ", line_count)
            if line_count > 10000000:
               break

print("Number of graph vertices: ", g.number_of_nodes())
print("Number of graph edges: ", g.number_of_edges())

for v in g.nodes():
    try:
        out= g.degree(v)
        if out == 0:
            g.remove_node(v)
    except ValueError:
        pass

print(g.number_of_nodes(), g.number_of_edges())

degrees = g.degree() # dictionary node:degree
in_values = sorted(set(degrees.values()))
list_degree_values = list(degrees.values())

in_hist = [list_degree_values.count(x) for x in in_values]
plt.figure()
plt.loglog(in_values, in_hist, 'ro', basex=2, basey=2)
#plt.xlim([1, 2**18])
plt.xlabel('Degree (log)')
plt.ylabel('Number of vertices(log)')
plt.title('Users network')
plt.grid(True)
plt.savefig('Degree_distribution_v1.png')

output_dir = './deg_dist/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

pickle.dump(in_values, open(output_dir + 'degree_sort.pickle', 'wb'))
pickle.dump(list_degree_values, open(output_dir + 'degree_list.pickle', 'wb'))
## DIFFUSION NETWORK CREATED .....................................

# sum_degrees = 0
# for i in degrees:
#     sum_degrees += degrees[i]
# avg_degree = sum_degrees / g.number_of_nodes()
# print("Avg degree:", avg_degree)
# g = g.to_undirected()
# avg_cc = nx.average_clustering(g)
# print("Avg CC:", avg_cc)
# num_connected_components = nx.number_connected_components(g)

# print("Num CC:", num_connected_components)


