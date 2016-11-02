# generative : edge weights are a probabilistic function based on the diffusion network and
#                 therefore are not independent - this file computes it.

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
import gc

diff_dict_07 = pickle.load(open('diffusion_dict_v1_t07.pickle', 'rb'))

print('Loading diffusion file...')
diff_file = 'rt_df.csv'
df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
df['mid'] = df['mid'].astype('str')
steep_inhib_times = pickle.load(open('steep_inhib_times.pickle', 'rb'))



def checkIsomorphism(graph_list, g):
    for gr in graph_list:
        if gtt.isomorphism(gr, g):
            return gr
    return False


def motif_weights_operation(mid):
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

    motif_patterns_list_4 = []
    motif_patterns_list_5 = []
    dict_patterns_4 = {}
    dict_patterns_5 = {}
    count_motifs_4 = 0
    count_motifs_5 = 0
    motifs_transition = [{} for i in range(number_intervals)]
    print("Cascade of mid: ", mid)

    # start_time = time.time()
    for i, r in cascade_set.iterrows():
        if new_nodes_count > 40 or (i == len(cascade_set['source']) - 1):
            if cnt_intervals >= 0:
                nodes_interval[cnt_intervals] = new_nodes

                if cnt_intervals == 0:
                    nodes_cur_interval = list(set(nodes_interval[cnt_intervals]))
                else:
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

                motifs_graph_filtered_prev, motifs_count_filtered_prev, vertex_maps_filtered_prev = \
                    gt.clustering.motifs(cascade_graph, 4, return_maps=True)

                # print("Number of 4 sized motifs: ", sorted(motifs_count_filtered_prev, reverse=True))

                for idx in range(len(motifs_graph_filtered_prev)):
                    pat_global = checkIsomorphism(motif_patterns_list_4, motifs_graph_filtered_prev[idx])
                    if pat_global == False:
                        motif_patterns_list_4.append(motifs_graph_filtered_prev[idx])
                        dict_patterns_4[motifs_graph_filtered_prev[idx]] = 'M4_' + str(count_motifs_4)
                        count_motifs_4 += 1

                motifs_4_previous_interval = {}
                for idx in range(len(vertex_maps_filtered_prev)):
                    map = vertex_maps_filtered_prev[idx]
                    if motifs_graph_filtered_prev[idx] not in motifs_4_previous_interval:
                        motifs_4_previous_interval[motifs_graph_filtered_prev[idx]] = []
                    for i in map:
                        motifs_4_previous_interval[motifs_graph_filtered_prev[idx]].append(set(i.a))

                # since there is no previous interval for first interval, we continue
                if cnt_intervals == 0:
                    cnt_intervals += 1
                    new_nodes = []
                    new_nodes_count = 0
                    motifs_prev_interval = motifs_4_previous_interval
                    continue

                motifs_graph_filtered_cur, motifs_count_filtered_cur, vertex_maps_filtered_cur = \
                    gt.clustering.motifs(cascade_graph, 5, return_maps=True)
                # print("Number of 5 sized motifs: ", sorted(motifs_count_filtered_cur, reverse=True))

                # Add new 5 sized motifs to motif list
                for idx in range(len(motifs_graph_filtered_cur)):
                    pat_global = checkIsomorphism(motif_patterns_list_5, motifs_graph_filtered_cur[idx])
                    if pat_global == False:
                        motif_patterns_list_5.append(motifs_graph_filtered_cur[idx])
                        dict_patterns_5[motifs_graph_filtered_cur[idx]] = 'M5_' + str(count_motifs_5)
                        count_motifs_5 += 1

                motifs_cur_interval = {}  # stores the motif list for current interval
                for idx in range(len(vertex_maps_filtered_cur)):
                    map = vertex_maps_filtered_cur[idx]
                    if motifs_graph_filtered_cur[idx] not in motifs_cur_interval:
                        motifs_cur_interval[motifs_graph_filtered_cur[idx]] = []
                    for i in map:
                        motifs_cur_interval[motifs_graph_filtered_cur[idx]].append(set(i.a))

                # Check for this interval, how many 4 sized motifs in previous interval have transitioned
                # into 5 sized motifs in this interval pattern by pattern.
                # Heavy computation - quadratic complexity - need to prune some motifs
                # Time complexity: # 4sized motifs patterns X # 5sized motifs patterns
                #                  X # instances of each 4sized patterns X # instances of each 4sized
                #                                                           patterns
                #
                # print('Number of motifs in this interval: ', len(motifs_cur_interval))
                # print('Number of motifs in previous interval: ', len(motifs_prev_interval))
                # print('Computing transition....')
                for m_prev in motifs_prev_interval:
                    for m_cur in motifs_cur_interval:
                        prop_map = gtt.subgraph_isomorphism(m_prev, m_cur)
                        # Filter 1: The 4 sized motif pattern must be a subgraph of the 5 sized motif
                        if len(prop_map) == 0:
                            continue
                        # Filter 2: The difference in the number of edges must be lesser than a threshold
                        num_edges_4 = m_prev.num_edges()
                        num_edges_5 = m_cur.num_edges()
                        if num_edges_5 - num_edges_4 > 3:
                            continue
                        list_motifs_cur = motifs_cur_interval[m_cur]
                        list_motifs_prev = motifs_prev_interval[m_prev]
                        # Filter 3: Do not consider some insignificant (by count) motifs
                        if len(list_motifs_cur) < 50 or len(list_motifs_prev) < 50:
                            continue
                        pat_4 = checkIsomorphism(motif_patterns_list_4, m_prev)
                        pat_4 = dict_patterns_4[pat_4]
                        pat_5 = checkIsomorphism(motif_patterns_list_5, m_cur)
                        pat_5 = dict_patterns_5[pat_5]

                        if pat_4 not in motifs_transition[cnt_intervals]:
                            motifs_transition[cnt_intervals][pat_4] = {}
                            motifs_transition[cnt_intervals][pat_4][pat_5] = 0
                        elif pat_5 not in motifs_transition[cnt_intervals][pat_4]:
                            motifs_transition[cnt_intervals][pat_4][pat_5] = 0
                        count_transition = 0
                        for idx_prev in range(len(list_motifs_prev)):
                            # Filter 4: consier only the first k=1000 motifs for each pattern combination both ways
                            if idx_prev > 1000:
                                break
                            for idx_cur in range(len(list_motifs_cur)):
                                if idx_cur > 1000:
                                    break
                                diff = list_motifs_cur[idx_cur] - list_motifs_prev[idx_prev]
                                #Filter 5: Consider transition only if the set interscetion count of nodes >=2 .
                                if len(diff) <= 2:
                                    # print(list_motifs_cur[idx_cur], list_motifs_prev[idx_prev], len(diff))
                                    count_transition += 1

                        motifs_transition[cnt_intervals][pat_4][pat_5] = count_transition

                motifs_prev_interval = motifs_4_previous_interval
                # At this point, the computation stops for the cascade for steep operation.
                # the values are stored in the reverse order for the plots.
                if stop_steep_flag == 1:
                    # print("Inside steep: ", len(motif_patterns_cascade_list))
                    motifs_transitions_steep = motifs_transition[:cnt_intervals + 1]

                    stop_steep_flag = 2
                    #cnt_intervals = 0
                    #motifs_transition = [{} for i in range(number_intervals)]

                if stop_inhib_flag == 1:
                    motifs_transition = motifs_transition[:cnt_intervals + 1]
                    return (motifs_transitions_steep, motifs_transition)

            cnt_intervals += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        if stop_inhib_flag == 2:
            return (motifs_transitions_steep, motifs_transition)

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
    gc.collect()
    number_intervals = 500
    count_motifs = 0

    motif_transition_global_list_steep = [{} for i in range(number_intervals)]
    motif_transition_global_list_inhib = [{} for i in range(number_intervals)]
    motif_patterns_list = []
    dict_patterns = {}

    numProcessors = 6
    pool = multiprocessing.Pool(numProcessors)

    num_cascades = len(steep_inhib_times.keys())

    print("Loading cascade data...")

    cnt_mids = 0

    # count_motifs = 0
    tasks = []
    for mid in steep_inhib_times:
        tasks.append( (mid) )
        cnt_mids += 1
        # if cnt_mids > 1000:
        #     break

    results = pool.map_async(motif_weights_operation, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            (motif_transition_steep, motif_transition_inhib) = motif_data[idx]
            cnt_interval_steep = len(motif_transition_steep)
            cnt_intervals_inhib = len(motif_transition_inhib)

            print(cnt_interval_steep, cnt_intervals_inhib)
            # Steep intervals operation
            for int_prev in range(cnt_interval_steep):
                interval = cnt_interval_steep - int_prev
                for m1 in motif_transition_steep[int_prev]:
                    if m1 not in motif_transition_global_list_steep[interval]:
                        motif_transition_global_list_steep[interval][m1] = {}
                    for m2 in motif_transition_steep[int_prev][m1]:
                        if m2 not in motif_transition_global_list_steep[interval][m1]:
                            motif_transition_global_list_steep[interval][m1][m2] = []
                        # print(motif_transition_steep[int_prev][m1][m2])
                        motif_transition_global_list_steep[interval][m1][m2].append(motif_transition_steep[int_prev][m1][m2])

            # inhib intervals operation
            for int_prev in range(cnt_intervals_inhib):
                interval = cnt_intervals_inhib - int_prev
                for m1 in motif_transition_inhib[int_prev]:
                    if m1 not in motif_transition_global_list_inhib[interval]:
                        motif_transition_global_list_inhib[interval][m1] = {}
                    for m2 in motif_transition_inhib[int_prev][m1]:
                        if m2 not in motif_transition_global_list_inhib[interval][m1]:
                            motif_transition_global_list_inhib[interval][m1][m2] = []
                        motif_transition_global_list_inhib[interval][m1][m2].append(
                            motif_transition_inhib[int_prev][m1][m2])

        except:
            continue

    for idx in range(1, 100):
        pickle.dump(motif_transition_global_list_steep[idx], open('motifs_transition/10_05/5/08/steep/transition_motifs_interval_'
                                                 + str(idx) + '.pickle', 'wb'))
        pickle.dump(motif_transition_global_list_inhib[idx],
                    open('motifs_transition/10_05/5/08/inhib/transition_motifs_interval_'
                         + str(idx) + '.pickle', 'wb'))
