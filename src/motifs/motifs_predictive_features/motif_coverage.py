# see the network coverge in the inhibition interval by motifs
# the algorithm is similar to the clique percolation method for
# finding "community" of motifs


from graph_tool.all import *
import graph_tool as gt
import graph_tool.stats as gts
import graph_tool.util as gtu
import graph_tool.draw as gtd
from pylab import *
from math import *
from numpy.random import *
import pickle
import pandas as pd
import os
import glob
import csv
import statistics as st
import random
from sqlalchemy import create_engine
from sqlalchemy import create_engine
import graph_tool.topology as gtt
import operator
import graph_tool.centrality as gtc
import time
import itertools
import csv
import resource
import multiprocessing
import threading

print('Loading files .... ')
diff_dict_07 = pickle.load(open('diffusion_dict_v1_t07.pickle', 'rb'))
diff_file = 'rt_df.csv'
df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
df['mid'] = df['mid'].astype('str')
steep_inhib_times = pickle.load(open('steep_inhib_times.pickle', 'rb'))


def checkIsomorphism(graph_list, g):
    for gr in graph_list:
        if gtt.isomorphism(gr, g):
            return gr
    return False


def mcoverage(mid, cnt_mids):
    cascade_set = df[(df['mid'] == str(mid))]

    # print('Cascade: ', cnt_mids, 'of size: ', len_cascade)
    steep_time = pd.to_datetime(steep_inhib_times[mid]['steep'])
    inhib_time = pd.to_datetime(steep_inhib_times[mid]['decay'])

    nodes_src = [[] for interval in range(500)]
    new_nodes = []
    new_nodes_count = 0

    weighted_edges = [[] for interval in range(500)]
    stop_steep_flag = 0
    stop_inhib_flag = 0
    cnt_intervals = 0

    nodes_interval = [[] for interval in range(500)]

    time_interval = [0 for interval in range(500)]
    last_time = 0

    motifs_edges_coverage = {}
    print("Cascade of mid: ", cnt_mids)

    # start_time = time.time()
    for i, r in cascade_set.iterrows():
        if new_nodes_count > 40 or (i == len(cascade_set['source']) - 1):
            # Only consider the inhibition interval
            if stop_inhib_flag != 1:
                # store the attributes of the previous interval
                nodes_interval[cnt_intervals] = new_nodes
                time_interval[cnt_intervals] = last_time
                cnt_intervals += 1
                new_nodes = []
                new_nodes_count = 0
                continue

            else: # if cnt_intervals >= 0:
                nodes_interval[cnt_intervals] = new_nodes
                nodes_cur_interval = list(set(nodes_interval[cnt_intervals] + nodes_interval[cnt_intervals - 1]))
                time_interval[cnt_intervals] = last_time

                # these are the inter-interval diffusion edges
                diff_edges_cur_interval = []
                for v in nodes_cur_interval:
                    try:
                        for uid, tid in diff_dict_07[v]:
                            if str(uid) in nodes_cur_interval:
                                if (v, uid) not in diff_edges_cur_interval:
                                    diff_edges_cur_interval.append(((v, uid), 1))
                                    # if len(list(set(diff_edges_cur_interval))) > 2*len_tree_edges:
                                    #     break
                    except:
                        continue

                edge_list = list(
                    set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals - 1] + diff_edges_cur_interval))

                # create the graph from the edge list
                cascade_graph = gt.Graph(directed=False)
                node_cascades_diff_map = {}
                cnt_nodes = 0
                cascade_vertices = cascade_graph.new_vertex_property("string")

                for (src, tgt), cur_time in edge_list:
                    if src not in node_cascades_diff_map:
                        node_cascades_diff_map[src] = cnt_nodes
                        v1 = cascade_graph.add_vertex()
                        cascade_vertices[v1] = src
                        cnt_nodes += 1
                    else:
                        v1 = cascade_graph.vertex(node_cascades_diff_map[src])

                    if tgt not in node_cascades_diff_map:
                        node_cascades_diff_map[tgt] = cnt_nodes
                        v2 = cascade_graph.add_vertex()
                        cascade_vertices[v2] = tgt
                        cnt_nodes += 1
                    else:
                        v2 = cascade_graph.vertex(node_cascades_diff_map[tgt])

                    if cascade_graph.edge(v1, v2):
                        continue
                    else:
                        cascade_graph.add_edge(v1, v2)

                gts.remove_self_loops(cascade_graph)
                gts.remove_parallel_edges(cascade_graph)

                # FINDING THE MOTIFS IN THE CASCADE GRAPH + DIFFUSION NETWORK
                motifs_graph_filtered, motifs_count_filtered, vertex_maps_filtered = \
                    gt.clustering.motifs(cascade_graph, 5, return_maps=True)

                # Check the coverage of the motifs - for now just the inhibition interval
                # Check for the individual motifs - not cumulative !!!!

                for idx_map in range(len(vertex_maps_filtered)):
                    map = vertex_maps_filtered[idx_map]

                    # find the motif instance with highest sum of degrees - seed motif
                    max_deg_sum = 0
                    for idx in map:
                        sum_deg = 0
                        for vert in idx.a:
                            sum_deg += cascade_graph.vertex(vert).out_degree()

                        if sum_deg > max_deg_sum:
                            max_deg_sum = sum_deg
                            max_motif_inst = idx

                    # print('max degree:', max_deg_sum)
                    # Starting with a seed motif, use CP algo to find the connected structure
                    nodes_covered = list(max_motif_inst.a)
                    edges_curr = []
                    edges_prev = []
                    motifs_covered = [max_motif_inst]
                    ind = 0
                    cnt_maps = 0
                    iter = 0

                    while True:
                        # print('Iteration: ', iter)
                        # In each of the iterations, find the ones common with those covered
                        for idx in map:
                            # Ignore the motifs covered
                            if idx in motifs_covered:
                                continue

                            cnt_maps += 1

                            # 1. The number of free nodes >= k --> at least k=[1, 4] free nodes
                            # 2. There should be at least one common node with nodes covered
                            number_nodes_free = len(idx.a) - len(list(set(idx.a).intersection(set(nodes_covered))))
                            if number_nodes_free < 2 or number_nodes_free == 5:
                                continue

                            # If an edge is already covered, do not count again
                            pairs_vert = list(itertools.combinations(idx.a, 2))
                            for (v1, v2) in pairs_vert:
                                e = cascade_graph.edge(v1, v2)
                                if e == None:
                                    continue
                                if ((v1, v2) in edges_prev) or ((v2, v1) in edges_prev):
                                    continue
                                edges_curr.append((v1, v2))
                                nodes_covered.append(v1)
                                nodes_covered.append(v2)

                            motifs_covered.append(idx)
                            if cnt_maps > 40000:
                                motifs_count_filtered[ind] = 40001
                                break

                        # Break out of this loop when there is no increase in edges covered
                        # in this iteration
                        if len(list(set(edges_prev))) == len(list(set(edges_curr))):
                            break
                        iter += 1
                        edges_prev = edges_curr

                    ind += 1
                    # print('Edges covered: ', len(edges_covered))
                    motifs_edges_coverage[motifs_graph_filtered[idx_map]] = (len(edges_curr) / cascade_graph.num_edges() * 100)
                    # motif_weights_patterns_cascade_list[cnt_intervals][motifs_graph_filtered[idx_map]] =\
                    #     total_weight / cnt_maps

                # At this point, the computation stops for the cascade for steep operation.
                # the values are stored in the reverse order for the plots.
                # if stop_steep_flag == 1:
                #     # print("Inside steep: ", len(motif_patterns_cascade_list))
                #     motif_weights_patterns_steep = motif_weights_patterns_cascade_list[:cnt_intervals + 1]
                #
                #     stop_steep_flag = 2
                #     motif_weights_patterns_cascade_list = [{} for i in range(500)]
                #     # cnt_intervals = 0

                if stop_inhib_flag == 1:
                    return motifs_edges_coverage

                    # motif_weights_patterns_cascade_list = motif_weights_patterns_cascade_list[:cnt_intervals + 1]
                    # return (motif_weights_patterns_steep, motif_weights_patterns_cascade_list)

            cnt_intervals += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        # if stop_inhib_flag == 2:
        #     return (motif_weights_patterns_steep, motif_weights_patterns_cascade_list)

        # for the steep and the inhibition point
        if rt_time >= steep_time and stop_steep_flag == 0:
            stop_steep_flag = 1
        if rt_time >= inhib_time and stop_inhib_flag == 0:
            stop_inhib_flag = 1

        src = r['source']
        tgt = r['target']

        rt_time = str(r['rt_time'])
        rt_date = rt_time[:10]
        rt_t = rt_time[11:19]
        record_time = rt_date + ' ' + rt_t
        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
        cur_time = time.mktime(time_x.timetuple())
        last_time = cur_time

        weighted_edges[cnt_intervals].append(((src, tgt), rt_time))

        nodes_src[cnt_intervals].append(r['source'])
        new_nodes.append(r['source'])
        new_nodes.append(r['target'])

        nodes_src[cnt_intervals] = list(set(nodes_src[cnt_intervals]))
        new_nodes = list(set(new_nodes))
        new_nodes_count = len(new_nodes)


if __name__ == '__main__':
    number_intervals = 500
    count_motifs = 0

    motifs_coverage = {}

    motif_patterns_list = []
    dict_patterns = {}

    numProcessors = 8
    pool = multiprocessing.Pool(numProcessors)

    num_cascades = len(steep_inhib_times.keys())

    print("Loading cascade data...")

    cnt_mids = 0

    # count_motifs = 0
    tasks = []
    for mid in steep_inhib_times:
        tasks.append( (mid, cnt_mids) )
        cnt_mids += 1
        if cnt_mids > 1500:
            break

    results = pool.starmap_async(mcoverage, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            motifs_coverage_values = motif_data[idx]
            # print(len(motifs_coverage_values))
            for m in motifs_coverage_values:
                pat = checkIsomorphism(motif_patterns_list, m)
                if pat == False:
                    motif_patterns_list.append(m)
                    dict_patterns[m] = 'M' + str(count_motifs)
                    count_motifs += 1

                    motifs_coverage[dict_patterns[m]] = []
                    motifs_coverage[dict_patterns[m]].append(motifs_coverage_values[m])
                else:
                    motifs_coverage[dict_patterns[pat]].append(motifs_coverage_values[m])

        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)
    pickle.dump(motifs_coverage, open('motifs_coveragt_min_2'
                                                  + '.pickle', 'wb'))
