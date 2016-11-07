# Increment the number of nodes at each iteration.
# then check the velocity of the centralities...
# check the magnitude and the direction of the rate of change of velocities

# Stop at the steep and the inhibition points....
# Add the diffusion edges at each stage....

# add the diffusion edges before the node in consideration
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
import time
import itertools
import csv
import networkx as nx
from operator import itemgetter
import operator

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

if __name__ == '__main__':
    #path = 'diffusion_dict_v1_t08.pickle'
    #diff_dict = pickle.load(open(path, 'rb'))
    diff_dict_07 = pickle.load(open('diffusion_dict_v1_t06_07.pickle', 'rb'))
    steep_inhib_times = pickle.load(open('steep_inhib_times.pickle', 'rb'))

    diff_file = 'rt_df.csv'
    df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn' ])
    df['mid'] = df['mid'].astype('str')

    diff_time_end = '2011-07-31 00:00:00'
    rt_date = diff_time_end[:10]
    rt_t = diff_time_end[11:19]
    record_time = rt_date + ' ' + rt_t
    time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
    diff_time_end = time.mktime(time_x.timetuple())

    # DSs for plots
    measure_time = []
    measure_interval = []
    print("Loading cascade data...")

    # DStructs for the cascade frames
    node_diff_map = {}
    cnt_nodes = 0
    mids_scanned = {}
    cnt_mids = 0

    number_intervals = 10

    entropy_global = [{} for interval in range(500)]
    entropy_global_velocity = [{} for interval in range(500)]

    entropy_global_inhib = [{} for interval in range(500)]
    entropy_global_velocity_inhib = [{} for interval in range(500)]

    time_end_interval_steep = {}
    time_end_interval_inhib = {}

    prev_nodes = []
    for mid in steep_inhib_times:
        cascade_set = df[(df['mid'] == str(mid))]
        len_cascade = len(list(set(cascade_set['source'] + cascade_set['target'])))
        # Bin the cascades into three types depending on the sizes.

        steep_time = pd.to_datetime(steep_inhib_times[mid]['steep'])
        inhib_time = pd.to_datetime(steep_inhib_times[mid]['decay'])

        rt_time = str(steep_time)
        rt_date = rt_time[:10]
        rt_t = rt_time[11:19]
        record_time = rt_date + ' ' + rt_t
        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
        steep_time_tuple = time.mktime(time_x.timetuple())

        rt_time = str(inhib_time)
        rt_date = rt_time[:10]
        rt_t = rt_time[11:19]
        record_time = rt_date + ' ' + rt_t
        time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
        inhib_time_tuple = time.mktime(time_x.timetuple())

        cnt_mids += 1

        # if cnt_mids > 2:
        #     break
        print('Cascade: ', cnt_mids, 'of size: ', len_cascade)

        nodes_src = [[] for interval in range(500)]
        new_nodes = []
        new_nodes_count = 0

        weighted_edges = [[] for interval in range(500)]
        stop_steep_flag = 0
        stop_inhib_flag = 0
        cnt_intervals = 0

        nodes_interval = [[] for interval in range(500)]
        entropy_velocity_cascade = [[] for interval in range(500)]
        entropy_cascade = [[] for interval in range(500)]
        entropy_accel_cascade = [[] for interval in range(500)]
        entropy_cascade_sum = [[] for interval in range(500)]
        number_edges = [0 for interval in range(500)]

        time_end_interval_steep[mid] = [[] for interval in range(500)]
        time_end_interval_inhib[mid] = [[] for interval in range(500)]
        time_interval = [0 for interval in range(500)] # keeps track of the end of the interval
        last_time = 0

        entropy_individuals = {}  # keeps track of the entropy of individuals by time

        time_measure_individual = {}
        # check the time difference between entropy change in the nodes
        # and then check the rates of the measure changes.
        nodes_time_intervals = [{} for interval in range(500)]
        max_time_nodes = {} # keeps track of the last retweet time of the node (diff. + cascade n/w)

        start_time_interval = [0 for interval in range(500)] # keeps track of the start time of each interval

        total_intervals = 0

        for i, r in cascade_set.iterrows():
            if new_nodes_count > 40:
                if cnt_intervals >= 0:
                    # print(cnt_intervals)
                    cascade_graph = nx.Graph()
                    nodes_interval[cnt_intervals] = new_nodes
                    nodes_src_interval = list(set(nodes_src[cnt_intervals]))

                    if cnt_intervals == 0:
                        nodes_cur_interval = list(set(nodes_interval[cnt_intervals]))
                        # nodes_src_interval = list(set(nodes_src[cnt_intervals]))
                    else:
                        nodes_cur_interval = list(set(nodes_interval[cnt_intervals] + nodes_interval[cnt_intervals-1]))
                        # nodes_src_interval = list(set(nodes_src[cnt_intervals] + nodes_src[cnt_intervals-1]))
                    time_interval[cnt_intervals] = last_time

                    # print(cnt_intervals, time_interval[cnt_intervals])
                    # these are the inter-interval diffusion edges
                    diff_edges_cur_interval = []
                    for v in nodes_cur_interval:
                        try:
                            # max_v_time = time_measure_individual[v].pop()
                            # time_measure_individual[v].push(max_v_time)
                            for uid, tid in diff_dict_07[v]:
                                if str(uid) in nodes_cur_interval:
                                    rt_time = str(tid)
                                    rt_date = rt_time[:10]
                                    rt_t = rt_time[11:19]
                                    record_time = rt_date + ' ' + rt_t
                                    time_x = datetime.datetime.strptime(record_time, '%Y-%m-%d %H:%M:%S')
                                    rt_time = time.mktime(time_x.timetuple())

                                    # diffusion edges before the retweet time of the node
                                    # if rt_time >= diff_end_time:
                                    #     continue

                                    time_measure_individual[v].push(rt_time)
                                    if (v, uid) not in diff_edges_cur_interval:
                                        diff_edges_cur_interval.append((v, uid))
                        except:
                            continue

                    edge_list = list(set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals-1] + diff_edges_cur_interval))
                    edge_list_final = []
                    for src, tgt in edge_list:
                        if src == tgt:
                            continue
                        edge_list_final.append((src, tgt))
                    cascade_graph.add_edges_from(edge_list_final)

                    # print('Cascade data', cascade_graph.number_of_nodes(), cascade_graph.number_of_edges())
                    number_edges[cnt_intervals] = cascade_graph.number_of_edges()
                    avg_neighbor_degree = nx.average_neighbor_degree(cascade_graph)

                    total_degree = {}
                    entropy_node = {}

                    # print('Interval: ', cnt_intervals)
                    #Compute the neighbor degree sum
                    for n in nodes_cur_interval:
                        entropy_node[n] = 0
                        if n not in cascade_graph.nodes():
                            continue
                        num_neighbors = len(cascade_graph.neighbors(n))
                        if num_neighbors == 0:
                            continue
                        neighbor_degree_sum = int(num_neighbors * avg_neighbor_degree[n])
                        if neighbor_degree_sum == 0:
                            continue
                        for nbr in cascade_graph.neighbors(n):
                            if len(cascade_graph.neighbors(nbr)) == 0:
                                continue
                            norm_nbr_deg = (len(cascade_graph.neighbors(nbr)))/neighbor_degree_sum
                            entropy_node[n] += (norm_nbr_deg*math.log(norm_nbr_deg, 10))
                        entropy_node[n] = -entropy_node[n]
                        # total_degree[n] = neighbor_degree_sum + num_neighbors
                        # entropy_node[n] = - num_neighbors*log((num_neighbors/total_degree[n]))

                    sorted_entropy = sorted(entropy_node.items(), key=operator.itemgetter(1), reverse=True)
                    top_k_entropy = sorted_entropy[:]

                    if cnt_intervals == 0:
                        prev_nodes = nodes_cur_interval
                        new_nodes_count = 0
                        new_nodes = []
                        cnt_intervals += 1
                        continue

                    count_individuals = 0
                    measure_sum = 0
                    measure_vel = 0
                    time_diff_sum = 0
                    for k, v in top_k_entropy:
                        if k not in prev_nodes:
                            continue
                        if k not in entropy_individuals:
                            entropy_individuals[k] = v
                            continue

                        prev_entropy = entropy_individuals[k]
                        entropy_individuals[k] = v # update the previous entropy for next interval as the current one

                        vel = (v - prev_entropy)
                        if vel == 0:
                            continue

                        measure_vel += vel
                        measure_sum += v
                        count_individuals += 1

                    if count_individuals == 0:
                        prev_nodes = nodes_cur_interval
                        new_nodes_count = 0
                        new_nodes = []
                        cnt_intervals += 1
                        continue
                    entropy_cascade_sum[cnt_intervals] = measure_sum/count_individuals
                    entropy_velocity_cascade[cnt_intervals] = measure_vel/count_individuals

                    if stop_steep_flag == 1:
                        for counter_prev in range(1, 100):
                            interval = cnt_intervals + 1 - counter_prev
                            try:
                                entropy_global[counter_prev][mid] = (entropy_cascade_sum[interval], number_edges[interval])
                                entropy_global_velocity[counter_prev][mid] = (entropy_velocity_cascade[interval], number_edges[interval])
                                # entropy_global_acceleration[counter_prev_5] = ((entropy_accel_cascade[interval], number_edges), mid)
                                # time_end_interval_steep[mid][counter_prev] = time_interval[interval]
                            except:
                                break
                        stop_steep_flag = 2

                        cnt_intervals = -1
                        entropy_velocity_cascade = [[] for interval in range(500)]
                        entropy_cascade_sum = [[] for interval in range(500)]
                        nodes_interval = [[] for interval in range(500)]
                        weighted_edges = [[] for interval in range(500)]
                        #print(cnt_intervals)

                    if stop_inhib_flag == 1:
                        for counter_prev in range(1, 100):
                            interval = cnt_intervals + 1 - counter_prev
                            try:
                                entropy_global_inhib[counter_prev][mid] = (entropy_cascade_sum[interval], number_edges[interval])
                                entropy_global_velocity_inhib[counter_prev][mid] = (entropy_velocity_cascade[interval], number_edges[interval])
                                # entropy_global_acceleration_inhib[counter_prev_5]= ((entropy_accel_cascade[interval], number_edges), mid)
                            except:
                                break

                        stop_inhib_flag = 2
                        cnt_intervals += 1
                        new_nodes = []
                        new_nodes_count = 0
                        break

                cnt_intervals += 1
                prev_nodes = nodes_cur_interval
                new_nodes = []
                new_nodes_count = 0
                total_intervals += 1

            rt_time = pd.to_datetime(str(r['rt_time']))

            if stop_inhib_flag == 2:
                break
            # for the inhibition point
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
            weighted_edges[cnt_intervals].append((src, tgt))

            if src not in time_measure_individual:
                time_measure_individual[src] = Stack()

            time_measure_individual[src].push(cur_time)

            nodes_src[cnt_intervals].append(r['source'])
            new_nodes.append(r['source'])
            new_nodes.append(r['target'])

            nodes_src[cnt_intervals] = list(set(nodes_src[cnt_intervals]))
            new_nodes = list(set(new_nodes))
            new_nodes_count = len(list(set(new_nodes)))

    output_dir = '11_05/entropy/entropy_change/v2/steep/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global_velocity, open(output_dir + 'entropy_diff' + '.pickle', 'wb'))

    output_dir = '11_05/entropy/entropy_change/v2/inhib/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global_velocity_inhib, open(output_dir + 'entropy_diff' + '.pickle', 'wb'))

    output_dir = '11_05/entropy/entropy_velocity/v2/steep/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global_velocity, open(output_dir + 'entropy_velocity' + '.pickle', 'wb'))

    output_dir = '11_05/entropy/entropy_values/v2/steep/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global, open(output_dir + 'entropy' + '.pickle', 'wb'))

    output_dir = '11_05/entropy/entropy_velocity/v2/inhib/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global_velocity_inhib, open(output_dir + 'entropy_velocity' + '.pickle', 'wb'))

    output_dir = '11_05/entropy/entropy_values/v2/inhib/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(entropy_global_inhib, open(output_dir + 'entropy' + '.pickle', 'wb'))
