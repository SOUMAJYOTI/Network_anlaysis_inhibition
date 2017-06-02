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

node_first_active_time = pickle.load(open('node_start_time_v2_t07.pickle', 'rb'))
node_tgt_degree = pickle.load(open('node_latest_degree_v2_t07.pickle', 'rb'))
tau = 0.7

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


def motif_weights_operation(mid, cnt_mids):
    tgt_degree_temp = node_tgt_degree
    node_first_active_time_temp = node_first_active_time
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

    motif_weights_patterns_cascade_list = [{} for i in range(500)]
    print("Cascade of mid: ", cnt_mids)

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
                len_tree_edges = len(list(set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals - 1])))
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
                edge_weights = cascade_graph.new_edge_property('float')  # keeps the weights

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

                    if cur_time != 1:
                        rt_time = str(cur_time)
                        rt_date = rt_time[:10]
                        rt_t = rt_time[11:19]
                        record_time = rt_date + ' ' + rt_t
                        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                        cur_time = time.mktime(time_x.timetuple())

                    if tgt not in tgt_degree_temp:
                        tgt_degree_temp[tgt] = 1
                        tgt_deg = 1
                    else:
                        tgt_deg = tgt_degree_temp[tgt]
                        tgt_degree_temp[tgt] += 1

                    if cur_time == 1: # diffusion edge
                        time_gap = 1
                    elif src not in node_first_active_time_temp:
                        node_first_active_time_temp[src] = cur_time
                        time_gap = 1 # this is some \eta that can be adjusted
                    else:
                        src_first_time = node_first_active_time_temp[src]

                        time_gap = ((int(cur_time - src_first_time) / 60) )

                    if cascade_graph.edge(v1, v2):
                        continue
                    else:
                        e = cascade_graph.add_edge(v1, v2)
                        try:
                            if cur_time == 1:
                                edge_weights[e] = 1
                            else:
                                edge_weights[e] = math.pow(time_gap, tau) * tgt_deg  # edge weight
                        except:
                            edge_weights[e] = tgt_deg

                gts.remove_self_loops(cascade_graph)
                gts.remove_parallel_edges(cascade_graph)

                # FINDING THE MOTIFS IN THE CASCADE GRAPH + DIFFUSION NETWORK
                motifs_graph_filtered, motifs_count_filtered, vertex_maps_filtered = \
                    gt.clustering.motifs(cascade_graph, 5, return_maps=True)

                for idx_map in range(len(vertex_maps_filtered)):
                    map = vertex_maps_filtered[idx_map]
                    ind = 0
                    total_weight = 0
                    cnt_maps = 0
                    for idx in map:
                        num_edges = 0
                        w = 1
                        pairs_vert = list(itertools.combinations(idx.a, 2))
                        for (v1, v2) in pairs_vert:
                            if cascade_graph.edge(v1, v2):
                                e = cascade_graph.edge(v1, v2)
                                w *= edge_weights[e]
                                num_edges += 1

                        try:
                            w = math.pow(w, (1./num_edges))
                            total_weight += w
                        except:
                            continue
                        cnt_maps += 1
                        if cnt_maps > 1000:
                            motifs_count_filtered[ind] = 1001
                            break

                    ind += 1
                    motif_weights_patterns_cascade_list[cnt_intervals][motifs_graph_filtered[idx_map]] =\
                        total_weight / cnt_maps
                # At this point, the computation stops for the cascade for steep operation.
                # the values are stored in the reverse order for the plots.
                if stop_steep_flag == 1:
                    # print("Inside steep: ", len(motif_patterns_cascade_list))
                    motif_weights_patterns_steep = motif_weights_patterns_cascade_list[:cnt_intervals + 1]

                    stop_steep_flag = 2
                    motif_weights_patterns_cascade_list = [{} for i in range(500)]
                    # cnt_intervals = 0

                if stop_inhib_flag == 1:
                    motif_weights_patterns_cascade_list = motif_weights_patterns_cascade_list[:cnt_intervals + 1]
                    return (motif_weights_patterns_steep, motif_weights_patterns_cascade_list)

            cnt_intervals += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        if stop_inhib_flag == 2:
            return (motif_weights_patterns_steep, motif_weights_patterns_cascade_list)
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

    motif_weights_patterns_global_list = [{} for i in range(number_intervals)]
    motif_weights_patterns_global_list_inhib = [{} for i in range(number_intervals)]
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
        if cnt_mids > 400:
            break

    results = pool.starmap_async(motif_weights_operation, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            (motif_weights_patterns_steep, motif_weights_patterns_inhib) = motif_data[idx]
            cnt_interval_steep = len(motif_weights_patterns_steep)
            cnt_intervals_inhib = len(motif_weights_patterns_inhib)

            print(cnt_interval_steep, cnt_intervals_inhib)
            # Steep intervals operation
            for int_prev in range(1, cnt_interval_steep+1):
                interval = cnt_interval_steep - int_prev

                # print(motif_patterns_steep[interval])
                for m in motif_weights_patterns_steep[interval]:
                    pat_global = checkIsomorphism(motif_weights_patterns_global_list[int_prev], m)
                    pat = checkIsomorphism(motif_patterns_list, m)
                    if pat == False:
                        motif_patterns_list.append(m)
                        dict_patterns[m] = 'M' + str(count_motifs)
                        count_motifs += 1

                    if pat_global == False:
                        motif_weights_patterns_global_list[int_prev][m] = []
                        motif_weights_patterns_global_list[int_prev][m].append(motif_weights_patterns_steep[interval][
                            m])
                    else:
                        # motif_weights_patterns_global_list[int_prev][pat_global] += \
                        #     motif_weights_patterns_steep[interval][m]
                        motif_weights_patterns_global_list[int_prev][pat_global].append(
                            motif_weights_patterns_steep[interval][m])

            # inhib intervals operation
            for int_prev in range(1, cnt_intervals_inhib+1):
                interval = cnt_intervals_inhib - int_prev
                for m in motif_weights_patterns_inhib[interval]:
                    # flag = 0

                    pat_global = checkIsomorphism(motif_weights_patterns_global_list_inhib[int_prev], m)
                    pat = checkIsomorphism(motif_patterns_list, m)
                    if pat == False:
                        motif_patterns_list.append(m)
                        dict_patterns[m] = 'M' + str(count_motifs)
                        count_motifs+= 1
                    if pat_global == False:
                        motif_weights_patterns_global_list_inhib[int_prev][m] = []
                        motif_weights_patterns_global_list_inhib[int_prev][m].append(motif_weights_patterns_inhib[interval][m])
                    else:
                        # motif_weights_patterns_global_list_inhib[int_prev][pat_global] += \
                        #     motif_weights_patterns_inhib[interval][m]
                        motif_weights_patterns_global_list_inhib[int_prev][pat_global].append(
                            motif_weights_patterns_inhib[interval][m])
        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)

    # print(count_motifs.value)
    weights_motifs_interval_dict = {}
    for int_prev in range(1, 100):
        try:
            for m in motif_weights_patterns_global_list[int_prev]:
                output_dir = 'motifs_patterns/5/interval_' + str(int_prev)
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                graph = checkIsomorphism(motif_patterns_list, m)
                str_val = dict_patterns[graph]
                output_string = output_dir + '/' + str_val + '.png'
                pos = gtd.arf_layout(m, max_iter=1000)
                gtd.graph_draw(m, pos=pos, output=output_string)

            for m in motif_weights_patterns_global_list[int_prev]:
                graph = checkIsomorphism(motif_patterns_list, m)
                str_val = dict_patterns[graph]
                weights_motifs_interval_dict[str_val] = motif_weights_patterns_global_list[int_prev][m]
            pickle.dump(weights_motifs_interval_dict, open('motifs_weights/generative/5/08/steep/weights_motifs_interval_'
                                                         + str(int_prev) + '.pickle', 'wb'))
            weights_motifs_interval_dict = {}
        except:
            continue

    weights_motifs_interval_dict = {}
    for int_prev in range(1, 100):
        try:
            for m in motif_weights_patterns_global_list_inhib[int_prev]:
                graph = checkIsomorphism(motif_patterns_list, m)
                str_val = dict_patterns[graph]
                weights_motifs_interval_dict[str_val] = motif_weights_patterns_global_list_inhib[int_prev][m]
            pickle.dump(weights_motifs_interval_dict, open('motifs_weights/generative/5/08/inhib/weights_motifs_interval_'
                                                     + str(int_prev) + '.pickle', 'wb'))
            weights_motifs_interval_dict = {}
        except:
            continue
