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

diff_dict_07 = pickle.load(open('diffusion_dict_v1_t07.pickle', 'rb'))

print('Loading diffusion file...')
diff_file = 'rt_df.csv'
df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
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

    nodes_interval = [[] for interval in range(500)]

    time_interval = [0 for interval in range(500)]
    last_time = 0
    time_measure_individual = {}

    steep_intervals = 0 # keeps track of the index of  steep interval
    zscores_list = [{} for i in range(500)]
    motifs_cascade_list = [{} for i in range(500)]

    pr_cascade_list = [{} for i in range(500)]
    bw_cascade_list = [{} for i in range(500)]
    degree_cascade_list = [{} for i in range(500)]
    entropy_cascade_list = [{} for i in range(500)]
    cl_cascade_list = [{} for i in range(500)]


    print("Cascade of mid: ", mid)

    # start_time = time.time()
    for i, r in cascade_set.iterrows():
        if new_nodes_count > 40 or (i == len(cascade_set['source']) - 1):
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

                # these are the inter-interval diffusion edges
                diff_edges_cur_interval = []
                len_tree_edges = len(list(set(weighted_edges[cnt_intervals] + weighted_edges[cnt_intervals - 1])))
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


                ''' The idea for computing the below centralities is this
                    Take the top 10 most central users and use their mean in that interval for prediction problems
                '''

                ''' Pagerank '''
                pr_val = gta.pagerank(cascade_graph)
                pr_list = pr_val.get_array()
                pr_sorted = sorted(pr_list, reverse=True)[:10]
                pr_cascade_list[cnt_intervals] = np.mean(np.array(pr_sorted))

                ''' Betweenness '''
                bw_val, bw_edge_val = gta.betweenness(cascade_graph)
                bw_list = bw_val.get_array()
                bw_sorted = sorted(bw_list, reverse=True)[:10]
                bw_cascade_list[cnt_intervals] = np.mean(np.array(bw_sorted))

                ''' Clustering '''
                cl_val = gta.local_clustering(cascade_graph)
                cl_list = cl_val.get_array()
                cl_sorted = sorted(cl_list, reverse=True)[:10]
                cl_cascade_list[cnt_intervals] = np.mean(np.array(cl_sorted))

                ''' Nodal degree '''
                deg_list = []
                for v in cascade_graph.vertices():
                    deg_list.append(v.out_degree())
                deg_sorted = sorted(deg_list, reverse=True)[:10]
                degree_cascade_list[cnt_intervals] = np.mean(np.array(deg_sorted))

                ''' Degree entropy '''
                deg_list = {}
                for v in cascade_graph.vertices():
                    deg_list[v]= int(v.out_degree())
                deg_sum = np.sum(np.array(list(deg_list.values())))

                deg_ent = []
                for v in cascade_graph.vertices():
                    deg_ent.append(-deg_list[v] * np.log(deg_list[v]/deg_sum))
                deg_ent_sorted = sorted(deg_ent, reverse=True)[:10]
                entropy_cascade_list[cnt_intervals] = np.mean(np.array(deg_ent_sorted))

                if stop_inhib_flag == 1:
                    # Add the cascade mid for indexing
                    return pr_cascade_list[steep_intervals:cnt_intervals+1], bw_cascade_list[steep_intervals:cnt_intervals+1], \
                           cl_cascade_list[steep_intervals:cnt_intervals+1], degree_cascade_list[steep_intervals:cnt_intervals+1],\
                           entropy_cascade_list[steep_intervals:cnt_intervals+1], mid

            cnt_intervals += 1
            new_nodes = []
            new_nodes_count = 0

        rt_time = pd.to_datetime(str(r['rt_time']))

        if stop_inhib_flag == 2:
            return
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
    pr_inhib = {}
    bw_inhib = {}
    cl_inhib = {}
    deg_inhib = {}
    en_inhib = {}

    # motif_patterns_list = []
    # dict_patterns = {}

    numProcessors = 5
    pool = multiprocessing.Pool(numProcessors, initializer=init, initargs=(count_motifs,))

    num_cascades = len(steep_inhib_times.keys())

    print("Loading cascade data...")

    cnt_mids = 0

    # count_motifs = 0
    tasks = []
    for mid in steep_inhib_times:
        tasks.append( (mid) )
        cnt_mids += 1
        if cnt_mids > 500:
            break

    results = pool.map_async(motif_operation, tasks)
    pool.close()
    pool.join()

    motif_data = results.get()

    count_invalid = 0
    for idx in range(len(motif_data)):
        try:
            pr_cascade_list, bw_cascade_list, cl_cascade_list, degree_cascade_list, \
            entropy_cascade_list, mid = motif_data[idx]
            cnt_intervals_inhib = len(pr_cascade_list)

            pr_inhib[mid] = [{} for i in range(number_intervals)]
            bw_inhib[mid] = [{} for i in range(number_intervals)]
            cl_inhib[mid] = [{} for i in range(number_intervals)]
            deg_inhib[mid] = [{} for i in range(number_intervals)]
            en_inhib[mid] = [{} for i in range(number_intervals)]
            # inhib intervals operation

            ''' Make the patterns global - add them to a dict '''
            for int_prev in range(1, cnt_intervals_inhib+1):
                interval = cnt_intervals_inhib - int_prev
                pr_inhib[mid][int_prev] = pr_cascade_list[interval]
                bw_inhib[mid][int_prev] = bw_cascade_list[interval]
                cl_inhib[mid][int_prev] = cl_cascade_list[interval]
                deg_inhib[mid][int_prev] = degree_cascade_list[interval]
                en_inhib[mid][int_prev] = entropy_cascade_list[interval]


        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)

    pickle.dump(pr_inhib, open('centralities/pr.pickle', 'wb'))
    pickle.dump(bw_inhib, open('centralities/bw.pickle', 'wb'))
    pickle.dump(cl_inhib, open('centralities/cl.pickle', 'wb'))
    pickle.dump(deg_inhib, open('centralities/deg.pickle', 'wb'))
    pickle.dump(en_inhib, open('centralities/en.pickle', 'wb'))
