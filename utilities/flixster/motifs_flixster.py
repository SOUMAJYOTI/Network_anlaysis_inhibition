import graph_tool as gt
import graph_tool.all as gta
import graph_tool.stats as gts
from pylab import *
from numpy.random import *
import pickle
import pandas as pd
import graph_tool.topology as gtt
import time
import multiprocessing

# Store the motif counts with the cascade ID for feature analysis and prediction problem !!!

diff_dict = pickle.load(open('data/flixster_sn_dict.pickle', 'rb'))

print('Loading diffusion file...')
diff_file = 'data/cascades_flixster_train.pickle'
df = pd.read_pickle(diff_file)
# df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
df['mid'] = df['mid'].astype('str')
steep_times, inhib_times = pickle.load(open('data/steep_inhib_times_dict_flixster.pickle', 'rb'))

print('Motifs for flixster....')

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

    # print(cascade_set)
    # print('Cascade: ', cnt_mids, 'of size: ', len_cascade)
    steep_time = pd.to_datetime(steep_times[mid])
    inhib_time = pd.to_datetime(inhib_times[mid])

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
    time_measure_individual = {}

    steep_intervals = 0 # keeps track of the index of  steep interval
    zscores_list = [{} for i in range(500)]
    motifs_cascade_list = [{} for i in range(500)]

    print("Cascade of mid: ", mid)

    # start_time = time.time()
    for i, r in cascade_set.iterrows():
        if new_nodes_count > 80 or (i == len(cascade_set['source']) - 1):
            # continue only after the steep interval
            if stop_steep_flag != 1:
                # store the attributes of the previous interval
                nodes_interval[cnt_intervals] = new_nodes
                time_interval[cnt_intervals] = last_time
                cnt_intervals += 1
                steep_intervals += 1
                new_nodes = []
                new_nodes_count = 0
                continue

            else : #cnt_intervals >= 0:
                nodes_interval[cnt_intervals] = new_nodes

                if cnt_intervals == 0:
                    nodes_cur_interval = list(set(nodes_interval[cnt_intervals]))
                else:
                    nodes_cur_interval = list(set(nodes_interval[cnt_intervals] + nodes_interval[cnt_intervals - 1]))
                time_interval[cnt_intervals] = last_time

                # print("Number of nodes: ", len(nodes_cur_interval))

                # these are the inter-interval diffusion edges
                diff_edges_cur_interval = []
                # len_tree_edges = len(list(set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals - 1])))
                for v in nodes_cur_interval:
                    try:
                        for uid in diff_dict[v]:
                            if uid in nodes_cur_interval:
                                # print('Hello')

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

                for (src, tgt) in edge_list:
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
                        e = cascade_graph.add_edge(v1, v2)

                gts.remove_self_loops(cascade_graph)
                gts.remove_parallel_edges(cascade_graph)

                print(len(list(cascade_graph.vertices())), len(list(cascade_graph.edges())))

                # FINDING THE MOTIFS IN THE CASCADE GRAPH + DIFFUSION NETWORK
                start_time = datetime.datetime.now()
                motifs_graph, motifs_count = \
                    gta.motifs(cascade_graph, 5)
                end_time = datetime.datetime.now()

                print((end_time - start_time).seconds/60)
                ''' Pre-emption time constraint '''
                if ((end_time - start_time).seconds / 60) >= 3: # minutes
                    return [], [], mid

                for idx_pat in range(len(motifs_graph)):
                    motifs_cascade_list[cnt_intervals][motifs_graph[idx_pat]] = motifs_count[idx_pat]

                ''' The following steps substitute the significance method in
                    the graph_tool and that method has some issues with rewiring
                 '''

                mpat_count = {}
                for idx_rw in range(3):
                    cascade_graph_copy = cascade_graph.copy()
                    ret = gta.random_rewire(cascade_graph_copy, "uncorrelated")
                    mt, mcount = gta.motifs(cascade_graph_copy, 5)

                    for idx_mt in range(len(mt)):
                        pat = checkIsomorphism(mpat_count.keys(), mt[idx_mt])

                        if pat == False:
                            mpat_count[mt[idx_mt]] = []
                            pat = mt[idx_mt]

                        mpat_count[pat].append(mcount[idx_mt])

                gr_list = mpat_count.keys()
                zscore = {}
                for idx_pat in range(len(motifs_graph)):
                    pat = checkIsomorphism(gr_list, motifs_graph[idx_pat])

                    if pat == False:
                        continue

                    mean_ran = np.mean(mpat_count[pat])
                    std_ran = np.std(mpat_count[pat])

                    if std_ran > 0:
                        zscore[pat] = (motifs_count[idx_pat] - mean_ran) / std_ran
                    else:
                        zscore[pat] = 0.

                zscores_list[cnt_intervals] = zscore
                if stop_inhib_flag == 1:
                    # Add the cascade mid for indexing
                    # return motifs_cascade_list[steep_intervals:cnt_intervals+1], mid
                    return motifs_cascade_list[steep_intervals:cnt_intervals+1], \
                           zscores_list[steep_intervals:cnt_intervals+1], mid

            cnt_intervals += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        if stop_inhib_flag == 2:
            return
        # for the steep and the inhibition point
        # print(rt_time, steep_time, inhib_time)
        if rt_time >= steep_time and stop_steep_flag == 0:
            stop_steep_flag = 1
        if rt_time >= inhib_time and stop_inhib_flag == 0:
            stop_inhib_flag = 1

        src = r['source']
        tgt = r['target']

        rt_time = str(r['rt_time'])
        cur_time = datetime.datetime.strptime(rt_time, '%Y-%m-%d %H:%M:%S')
        last_time = cur_time

        weighted_edges[cnt_intervals].append((src, tgt))

        if src not in time_measure_individual:
            time_measure_individual[src] = Stack()

        # time_measure_individual[src].push(cur_time)

        nodes_src[cnt_intervals].append(r['source'])
        new_nodes.append(r['source'])
        new_nodes.append(r['target'])

        nodes_src[cnt_intervals] = list(set(nodes_src[cnt_intervals]))

        ''' Check the cumulative of the number of new nodes in \mathcal{N} >= 80 ---> break'''
        if cnt_intervals > 0:
            new_nodes = list(set(new_nodes).union(set(nodes_src[cnt_intervals-1])))
            new_nodes_count = len(new_nodes)
        else:
            new_nodes = list(set(new_nodes))
            new_nodes_count = len(new_nodes)


if __name__ == '__main__':
    global count_motifs
    count_motifs = multiprocessing.Value('i', 0)

    number_intervals = 500
    zscores_global_list_inhib = {}
    mcounts_global_list_inhib = {}

    # motif_patterns_list = []
    # dict_patterns = {}

    numProcessors = 5
    pool = multiprocessing.Pool(numProcessors, initializer=init, initargs=(count_motifs,))

    num_cascades = len(steep_times.keys())

    print("Loading cascade data...")

    cnt_mids = 0

    # count_motifs = 0
    tasks = []

    for mid in steep_times:
        if mid == 9218:
            tasks.append( (mid) )
            break
        cnt_mids += 1
        # if cnt_mids > 400:
        #     break

    results = pool.map_async(motif_operation, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            mcounts, zscores,  mid = motif_data[idx]\

            ''' Condition for pre-emption because of time limits'''
            if len(mcounts) == 0:
                continue
            cnt_intervals_inhib = len(mcounts)

            zscores_global_list_inhib[mid] = [{} for i in range(number_intervals)]
            mcounts_global_list_inhib[mid] = [{} for i in range(number_intervals)]
            # inhib intervals operation

            ''' Make the patterns global - add them to a dict '''
            for int_prev in range(1, cnt_intervals_inhib+1):
                interval = cnt_intervals_inhib - int_prev
                for m in zscores[interval]:
                    # if not checkIsomorphism(motif_patterns_list, m):
                    #     motif_patterns_list.append(m)
                    zscores_global_list_inhib[mid][int_prev][m] = zscores[interval][m]

                for m in mcounts[interval]:
                    mcounts_global_list_inhib[mid][int_prev][m] = mcounts[interval][m]
        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)

    pickle.dump(zscores_global_list_inhib, open('data/motifs_z_scores_degreeSeq_k5_flixster.pickle', 'wb'))
    pickle.dump(mcounts_global_list_inhib, open('data/motif_counts_k5_flixster.pickle', 'wb'))