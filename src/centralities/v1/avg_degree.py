# Increment the number of nodes at each iteration.
# then check the velocity of the centralities...
# check the magnitude and the direction of the rate of change of velocities
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

class Queue:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def enqueue(self, item):
        self.items.insert(0,item)

    def dequeue(self):
        return self.items.pop()

    def size(self):
        return len(self.items)

if __name__ == '__main__':
    diff_dict_07 = pickle.load(open('diffusion_dict_v1_t06_07.pickle', 'rb'))
    steep_inhib_times = pickle.load(open('steep_inhib_times.pickle', 'rb'))

    diff_file = 'rt_df.csv'
    df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
    df['mid'] = df['mid'].astype('str')

    print("Loading cascade data...")

    # DStructs for the cascade frames
    node_diff_map = {}
    cnt_nodes = 0
    mids_scanned = {}
    cnt_mids = 0

    number_intervals = 10

    measure_global = [{} for interval in range(500)]
    measure_global_velocity = [{} for interval in range(500)]

    measure_global_inhib = [{} for interval in range(500)]
    measure_global_velocity_inhib = [{} for interval in range(500)]

    for mid in steep_inhib_times:
        cascade_set = df[(df['mid'] == str(mid))]
        len_cascade = len(list(set(cascade_set['source'] + cascade_set['target'])))

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

        if cnt_mids > 1500:
            break
        print('Cascade: ', cnt_mids, 'of size: ', len_cascade)

        nodes_src = [[] for interval in range(500)]
        new_nodes = []
        new_nodes_count = 0

        weighted_edges = [[] for interval in range(500)]
        stop_steep_flag = 0
        stop_inhib_flag = 0
        stop_flag = 0
        cnt_intervals = 0

        nodes_interval = [[] for interval in range(500)]
        measure_velocity_cascade = [[] for interval in range(500)]
        measure_cascade = [[] for interval in range(500)]
        measure_cascade_sum = [[] for interval in range(500)]

        time_interval = [0 for interval in range(500)]
        last_time = 0

        measure_individuals = {}  # keeps track of the entropy of individuals by time

        time_measure_individual = {}
        # check the time difference between entropy change in the nodes
        # and then check the rates of the measure changes.
        nodes_time_intervals = [{} for interval in range(500)]
        max_time_nodes = {} # keeps track of the last retweet time of the node (diff. + cascade n/w)

        start_time_interval = [0 for interval in range(500)] # keeps track of the start time of each interval
        for i, r in cascade_set.iterrows():
            if new_nodes_count > 40:
                if cnt_intervals >= 0:
                    #print(cnt_intervals)
                    cascade_graph = nx.Graph()
                    nodes_interval[cnt_intervals] = new_nodes
                    #nodes_src_interval = list(set(nodes_src[cnt_intervals]))

                    if cnt_intervals == 0:
                        nodes_cur_interval = list(set(nodes_interval[cnt_intervals]))
                        nodes_src_interval = list(set(nodes_src[cnt_intervals]))
                    else:
                        nodes_cur_interval = list(set(nodes_interval[cnt_intervals] + nodes_interval[cnt_intervals-1]))
                        nodes_src_interval = list(set(nodes_src[cnt_intervals] + nodes_src[cnt_intervals-1]))
                    time_interval[cnt_intervals] = last_time

                    # these are the inter-interval diffusion edges
                    diff_edges_cur_interval = []
                    for v in nodes_cur_interval:
                        try:
                            max_v_time = time_measure_individual[v].pop()
                            time_measure_individual[v].push(max_v_time)
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
                        except:
                            continue

                    edge_list = list(set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals-1] + diff_edges_cur_interval))
                    cascade_graph.add_edges_from(edge_list)
                    #print(cascade_graph.number_of_nodes(), cascade_graph.number_of_edges(), len(diff_edges_cur_interval))

                    number_edges = cascade_graph.number_of_edges()
                    number_nodes = cascade_graph.number_of_nodes()

                    measure_node = {}
                    for n in nodes_cur_interval:
                        measure_node[n] = len(cascade_graph.neighbors(n))

                    sorted_measure = sorted(measure_node.items(), key=operator.itemgetter(1), reverse=True)
                    top_k_measure = sorted_measure[:]

                    if cnt_intervals == 0:
                        new_nodes_count = 0
                        new_nodes = []
                        cnt_intervals += 1
                        continue

                    count_individuals = 0
                    measure_sum = 0
                    measure_vel = 0
                    for k, v in top_k_measure:

                        # nodes which are leaves are not considered for the measures
                        if k not in nodes_src_interval:
                            measure_individuals[k] = v
                            continue

                        # for velocity, we need at least two comparisons
                        if time_measure_individual[k].size() < 2:
                            measure_individuals[k] = v
                            continue

                        if k not in measure_individuals:
                            measure_individuals[k] = v
                            continue

                        prev_measure = measure_individuals[k]
                        measure_individuals[k] = v # update the previous entropy for next interval as the current one

                        vel = (v - prev_measure)
                        # accel = vel/time_diff
                        #print(vel)
                        if vel == 0:
                            continue

                        measure_sum += v
                        measure_vel += vel
                        count_individuals += 1
                        # measure_velocity_cascade[cnt_intervals].append(vel)

                    if count_individuals == 0:
                        new_nodes_count = 0
                        new_nodes = []
                        cnt_intervals += 1
                        continue
                    measure_cascade_sum[cnt_intervals] = measure_sum/count_individuals
                    measure_velocity_cascade[cnt_intervals] = measure_vel/count_individuals

                    if stop_steep_flag == 1:
                        for counter_prev_5 in range(1, 100):
                            interval = cnt_intervals + 1 - counter_prev_5
                            try:
                                measure_global[counter_prev_5][mid].append((measure_cascade_sum[interval], number_edges))
                                measure_global_velocity[counter_prev_5][mid].append((measure_velocity_cascade[interval], number_edges))
                            except:
                                break
                        stop_steep_flag = 2

                        cnt_intervals = 0
                        measure_velocity_cascade = [[] for interval in range(500)]
                        measure_cascade_sum = [[] for interval in range(500)]

                    if stop_inhib_flag == 1:
                        for counter_prev_5 in range(1, 100):
                            interval = cnt_intervals + 1 - counter_prev_5
                            try:
                                measure_global_inhib[counter_prev_5][mid].append((measure_cascade_sum[interval], number_edges))
                                measure_global_velocity_inhib[counter_prev_5][mid].append((measure_velocity_cascade[interval], number_edges))
                            except:
                                break
                        stop_inhib_flag = 2
                        new_nodes_count = 0
                        new_nodes = []
                        cnt_intervals += 1
                        break


                cnt_intervals += 1
                new_nodes = []
                new_nodes_count = 0

            rt_time = pd.to_datetime(str(r['rt_time']))

            # for the inhibition point
            # if rt_time <= steep_time:
            #     continue
            # if rt_time == inhib_time:
            #      stop_flag = 1

            if stop_inhib_flag == 2:
                break
            # for the steep point
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

            if start_time_interval[cnt_intervals] == 0:
                start_time_interval[cnt_intervals] = cur_time

            last_time = cur_time
            new_nodes_count = len(list(set(new_nodes)))
            weighted_edges[cnt_intervals].append((src, tgt))

            if src not in time_measure_individual:
                time_measure_individual[src] = Stack()

            time_measure_individual[src].push(cur_time)

            nodes_src[cnt_intervals].append(r['source'])
            new_nodes.append(r['source'])
            new_nodes.append(r['target'])

            nodes_src[cnt_intervals] = list(set(nodes_src[cnt_intervals]))
            new_nodes = list(set(new_nodes))

    output_dir = '11_05/deg/deg_change/steep/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(measure_global_velocity, open(output_dir + 'deg_diff' + '.pickle', 'wb'))

    output_dir = '11_05/deg/deg_values/steep/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(measure_global, open(output_dir + 'deg' + '.pickle', 'wb'))

    output_dir = '11_05/deg/deg_change/inhib/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(measure_global_velocity_inhib, open(output_dir + 'deg_diff' + '.pickle', 'wb'))

    output_dir = '11_05/deg/deg_values/inhib/'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    pickle.dump(measure_global_inhib, open(output_dir + 'deg' + '.pickle', 'wb'))