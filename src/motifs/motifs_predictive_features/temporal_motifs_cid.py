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
import graph_tool.topology as gtt
import operator
import graph_tool.centrality as gtc
import time
import itertools
import csv
import resource
import time
import multiprocessing
import threading
import random

# The idea is to implement parallel processing for the cascade motif operations !!!

diff_dict_07 = pickle.load(open('diffusion_dict_v1_t06_07.pickle', 'rb'))

print('Loading diffusion file...')
diff_file = 'rt_df.csv'
df = pd.read_csv(diff_file, names=['target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
df['mid'] = df['mid'].astype('str')
steep_inhib_times = pickle.load(open('steep_inhib_times.pickle', 'rb'))

class Stack:
  def __init__(self):
    self.__storage = []

  def isEmpty(self):
    return len(self.__storage) == 0

  def push(self,p):
    self.__storage.append(p)

  def pop(self):
    return self.__storage.pop()

  def size(self):
    return len(self.__storage)

  def return_storage(self):
      return self.__storage

def init(args):
    global count_motifs
    count_motifs = args

def checkIsomorphism(graph_list, g):
    for gr in graph_list:
        if gtt.isomorphism(gr, g):
            return gr
    return False

def motif_operation(mid):
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
    window = 1

    nodes_interval = [[] for interval in range(500)]

    time_interval = [0 for interval in range(500)]
    last_time = 0
    time_measure_individual = {}
    parent_node_time = {}
    parent_node = {}

    motif_frontiers_count = [{} for i in range(500)]
    motif_non_frontiers_count = [{} for i in range(500)]
    motif_frontiers_gt_count = [{} for i in range(500)]
    return_list = []
    print("Cascade of mid: ", mid)

    steep_intervals = 0
    # start_time = time.time()
    for i, r in cascade_set.iterrows():
        if new_nodes_count > 40:
            if stop_steep_flag != 1:
                # store the attributes of the previous interval
                nodes_interval[cnt_intervals] = new_nodes
                time_interval[cnt_intervals] = last_time
                cnt_intervals += 1
                steep_intervals += 1
                new_nodes = []
                new_nodes_count = 0
                continue
            else:
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
                                rt_time = str(tid)
                                rt_date = rt_time[:10]
                                rt_t = rt_time[11:19]
                                record_time = rt_date + ' ' + rt_t
                                time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                                rt_time = time.mktime(time_x.timetuple())

                                time_measure_individual[v].push(rt_time)
                                if (v, uid) not in diff_edges_cur_interval:
                                    diff_edges_cur_interval.append((v, uid))
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
                frontier_vertices = []
                frontier_gt_vertices = []
                non_frontier_vertices = []

                for (src, tgt) in edge_list:
                    if src == tgt:
                        continue
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

                    # print(parent_node_time[src], parent_node_time[tgt])
                    time_diff = (parent_node_time[tgt] - parent_node_time[src])/60

                    if time_diff >= 0:
                        if time_diff < 15*window: # lambda
                            if str(src) in diff_dict_07:
                                flag_frontier = 0
                                for nbr, tid in diff_dict_07[str(src)]:
                                    if str(tgt) == nbr:
                                        frontier_vertices.append(tgt)
                                        flag_frontier = 1
                                        break
                                if flag_frontier == 0:
                                    non_frontier_vertices.append(tgt)
                        else:
                            if str(src) in diff_dict_07:
                                for nbr, tid in diff_dict_07[str(src)]:
                                    if str(tgt) == nbr:
                                        frontier_gt_vertices.append(tgt)
                                        break

                    if cascade_graph.edge(v1, v2):
                        continue
                    else:
                        cascade_graph.add_edge(v1, v2)

                gts.remove_self_loops(cascade_graph)
                gts.remove_parallel_edges(cascade_graph)

                # pos = gt.draw.sfdp_layout(cascade_graph)
                # gt.draw.graph_draw(cascade_graph, pos=pos, output='c_'+str(window)+'.pdf')

                # FINDING THE MOTIFS IN THE CASCADE GRAPH + DIFFUSION NETWORK
                motifs_graph_filtered, motifs_count_filtered, vertex_maps_filtered = \
                    gt.clustering.motifs(cascade_graph, 5, return_maps=True)

                frontier_motifs = {}
                non_frontier_motifs = {}
                frontiers_gt_motifs = {}
                for idx_g in range(len(motifs_graph_filtered)):
                    map = vertex_maps_filtered[idx_g]
                    cnt_maps = 0
                    ind = 0
                    cnt_frontier_motifs = 0
                    cnt_frontiers_gt_motifs = 0
                    cnt_adopter_motifs = 0

                    for idx in map:
                        flag_frontier = 0
                        flag_frontier_gt = 0
                        cnt_maps += 1
                        for n in list(idx.a):
                            user = cascade_vertices[n]
                            if user in frontier_vertices:
                                flag_frontier = 1
                                break
                            if user in frontier_gt_vertices:
                                flag_frontier_gt = 1
                        if flag_frontier == 1:
                            cnt_frontier_motifs += 1
                        elif flag_frontier_gt == 1:
                            cnt_frontiers_gt_motifs += 1
                        else:
                            cnt_adopter_motifs += 1

                        if cnt_maps > 3000:
                            motifs_count_filtered[ind] = 3001
                            break
                    ind += 1
                    # print(cnt_intervals, 'Motifs: ', motifs_count_filtered[idx_g], cnt_frontier_motifs, cnt_frontiers_gt_motifs, cnt_adopter_motifs)
                    frontier_motifs[idx_g] = cnt_frontier_motifs
                    frontiers_gt_motifs[idx_g] = cnt_frontiers_gt_motifs
                    non_frontier_motifs[idx_g] = cnt_adopter_motifs

                counter_lock = threading.Lock()
                with counter_lock:
                    for idx in range(len(motifs_graph_filtered)):
                        motif_frontiers_count[cnt_intervals][motifs_graph_filtered[idx]] = frontier_motifs[idx]
                        motif_non_frontiers_count[cnt_intervals][motifs_graph_filtered[idx]] = non_frontier_motifs[idx]
                        motif_frontiers_gt_count[cnt_intervals][motifs_graph_filtered[idx]] = frontiers_gt_motifs[idx]

                # print(motif_patterns_cascade_list)

                if stop_inhib_flag == 1:
                    motif_frontiers_inhib = motif_frontiers_count[:cnt_intervals + 1], mid
                    motif_non_frontiers_inhib = motif_non_frontiers_count[:cnt_intervals + 1], mid
                    motif_frontiers_gt_inhib = motif_frontiers_gt_count[:cnt_intervals + 1], mid
                    return_list = [motif_frontiers_inhib, motif_non_frontiers_inhib, motif_frontiers_gt_inhib]
                    return return_list

            cnt_intervals += 1
            window += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        if stop_inhib_flag == 2:
            return return_list
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

        if src not in parent_node_time:
            parent_node_time[src] = cur_time
        parent_node_time[tgt] = cur_time
        parent_node[tgt] = src

        weighted_edges[cnt_intervals].append((src, tgt))

        if src not in time_measure_individual:
            time_measure_individual[src] = Stack()

        time_measure_individual[src].push(cur_time)

        nodes_src[cnt_intervals].append(r['source'])
        new_nodes.append(r['source'])
        new_nodes.append(r['target'])

        nodes_src[cnt_intervals] = list(set(nodes_src[cnt_intervals]))
        new_nodes = list(set(new_nodes))
        new_nodes_count = len(new_nodes)


if __name__ == '__main__':
    global count_motifs
    count_motifs = multiprocessing.Value('i', 0)

    number_intervals = 500


    motif_patterns_global_list_1_inhib = {}
    motif_patterns_global_list_2_inhib = {}
    motif_patterns_global_list_3_inhib = {}
    frontiers_count_global_list_inhib = {}
    frontiers_gt_count_global_list_inhib = {}
    adopters_count_global_list_inhib = {}

    motif_patterns_list = []
    dict_patterns = {}

    numProcessors = 6
    pool = multiprocessing.Pool(numProcessors, initializer=init, initargs=(count_motifs,))

    num_cascades = len(steep_inhib_times.keys())

    print("Loading cascade data...")

    cnt_mids = 0

    # count_motifs = 0
    tasks = []
    for mid in steep_inhib_times:
        tasks.append( (mid) )
        cnt_mids += 1
        # if cnt_mids > 3:
        #     break

    results = pool.map_async(motif_operation, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            mf_inhib, mfa_inhib, mfg_inhib = motif_data[idx]
            motif_frontiers_inhib = mf_inhib[0]
            motif_adopters_inhib = mfa_inhib[0]
            motif_frontiers_gt_inhib = mfg_inhib[0]
            mid = mf_inhib[1]
            cnt_interval_inhib = len(motif_frontiers_inhib)

            motif_patterns_global_list_1_inhib[mid] = [{} for i in range(number_intervals)]
            motif_patterns_global_list_2_inhib[mid] = [{} for i in range(number_intervals)]
            motif_patterns_global_list_3_inhib[mid] = [{} for i in range(number_intervals)]
            frontiers_count_global_list_inhib[mid] = [{} for i in range(number_intervals)]
            frontiers_gt_count_global_list_inhib[mid] = [{} for i in range(number_intervals)]
            adopters_count_global_list_inhib[mid] = [{} for i in range(number_intervals)]

            # Inhib intervals operation
            for int_prev in range(1, cnt_interval_inhib + 1):
                interval = cnt_interval_inhib - int_prev

                # print(motif_patterns_steep[interval])
                # We store the intrevals data in reverse order from inhib interval to the first interval
                for m in motif_frontiers_inhib[interval]:
                    frontiers_count_global_list_inhib[mid][int_prev][m] = motif_frontiers_inhib[interval][m]

            for int_prev in range(1, cnt_interval_inhib + 1):
                interval = cnt_interval_inhib - int_prev

                # print(motif_patterns_steep[interval])
                # We store the intrevals data in reverse order from inhib interval to the first interval
                for m in motif_frontiers_gt_inhib[interval]:
                    frontiers_gt_count_global_list_inhib[mid][int_prev][m] = motif_frontiers_gt_inhib[interval][m]

            for int_prev in range(1, cnt_interval_inhib + 1):
                interval = cnt_interval_inhib - int_prev

                # print(motif_patterns_steep[interval])
                # We store the intrevals data in reverse order from inhib interval to the first interval
                for m in motif_adopters_inhib[interval]:
                    adopters_count_global_list_inhib[mid][int_prev][m] = motif_adopters_inhib[interval][m]

        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)

    pickle.dump(frontiers_count_global_list_inhib, open('frontiers_inhib_count.pickle', 'wb'))
    pickle.dump(frontiers_gt_count_global_list_inhib, open('frontiers_gt_inhib_count.pickle', 'wb'))
    pickle.dump(adopters_count_global_list_inhib, open('adopters_inhib_count.pickle', 'wb'))
