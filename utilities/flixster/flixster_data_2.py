import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import datetime
from dateutil.relativedelta import relativedelta
import multiprocessing

movie_file = 'data/Ratings_timed_v1.txt'
movie_df = pd.read_csv(movie_file, sep="\s+", encoding="ISO-8859-1", index_col=False, )
movie_df['date'] = pd.to_datetime(movie_df['date'])


def cascade_op(mid, cnt_mids):
    ''' Storage structures for cascade dataframe'''
    srcList = []
    tgtList = []
    rtList = []
    midList = []
    postTimeList = []
    isOnList = []

    movie_cascade = movie_df[movie_df['movieid'] == mid]

    movie_cascade = movie_cascade.sort_values(by='date')

    ''' Create temporary lists for traversal - THINK OF A BETTER WAY OF DOING THIS '''
    uid_tempList = list(movie_cascade['userid'])
    time_tempList = list(movie_cascade['date'])

    ''' !. Create a link between current node and a future node if the time difference is less than 7 days '''
    ''' 2. If the link between the node and the future node is repeated, ignore it. '''
    users = set([])
    for idx_cur in range(len(uid_tempList) - 1, 0, -1):
        curr_time = time_tempList[idx_cur]
        cur_time_str = datetime.datetime.strptime(curr_time.strftime("%Y-%m-%d"), '%Y-%m-%d')

        users.add(uid_tempList[idx_cur])
        cnt_idx_link = 0
        for idx_link in range(idx_cur - 1, -1, -1):
            cnt_idx_link += 1
            prev_time = time_tempList[idx_link]
            prev_time_str = datetime.datetime.strptime(prev_time.strftime("%Y-%m-%d"), '%Y-%m-%d')

            diff = cur_time_str - prev_time_str

            # print(cur_time_str, prev_time_str)

            if diff.days > 1 :  # This is the delta parameter fixed manually
                break

            users.add(uid_tempList[idx_link])

            tgtList.append(uid_tempList[idx_cur])
            srcList.append(uid_tempList[idx_link])

            midList.append(mid)
            rtList.append(time_tempList[idx_cur])
            postTimeList.append(time_tempList[idx_link])
            isOnList.append(1.)

        if len(users) > 2000:
            break
    print("Movie no.: ", cnt_mids, " with mid: ", mid, " with size: ", len(movie_cascade), ' with source size: ', len(srcList))

    return srcList, tgtList, rtList, midList, postTimeList, isOnList

def create_cascades():
    '''

    Create the cascade files for all movie ids
    :return:
    '''

    ''' Create a mapping of cascade source times '''

    print('Create start mappings.....')
    ''' These are the start and end tiems for analysis'''
    START_TIME = pd.to_datetime('2005-03-01')
    END_TIME = pd.to_datetime('2010-01-01')

    movie_cascades = movie_df[movie_df['date'] >= START_TIME]
    # movie_cascades = movie_cascades[movie_cascades['date'] < END_TIME]

    movie_id_list = list(set(movie_cascades['movieid']))  # Set of movie ids

    cnt_movies = 0
    movie_id_filter = []

    for mid in movie_id_list:
        movie_cascade = movie_df[movie_df['movieid'] == mid]
        movie_cascade = movie_cascade.sort_values(by='date')

        end_times_cascade = movie_cascade.iloc[len(movie_cascade)-1]['date']

        if end_times_cascade <= END_TIME:
            if len(movie_cascade) >= 1000:  # RESHARE MIN SIZE = 200 ( THE CASCADE MIN SIZE MAY BE < 200)
                movie_id_filter.append(mid)
                cnt_movies += 1
                print("Movie no: ", cnt_movies)

        if cnt_movies > 1000:
            break

    return movie_id_filter


if __name__ == '__main__':
    number_intervals = 500
    srcList_global = []
    tgtList_global = []
    rtList_global = []
    midList_global = []
    postTimeList_global = []
    isOn_global = []

    numProcessors = 5
    pool = multiprocessing.Pool(numProcessors)

    print("Loading cascade data...")

    cnt_mids = 0
    tasks = []
    # mids = create_cascades()
    #
    # pickle.dump(mids, open('data/moviesSample.pickle', 'wb'))

    mids = pickle.load(open('data/moviesSample.pickle', 'rb'))
    for mid in mids:
        tasks.append( (mid, cnt_mids) )
        cnt_mids += 1
        # if cnt_mids > 1:
        #     break

    results = pool.starmap_async(cascade_op, tasks)
    pool.close()
    pool.join()

    cascade_data = results.get()

    count_invalid = 0
    for idx in range(len(cascade_data)):
        try:
            srcList, tgtList, rtList, midList, postTimeList, isOnList = cascade_data[idx]

            srcList_global.extend(srcList)
            tgtList_global.extend(tgtList)
            rtList_global.extend(rtList)
            midList_global.extend(midList)
            postTimeList_global.extend(postTimeList)
            isOn_global.extend(isOnList)

        except:
            count_invalid += 1

    print('Invalid: ', count_invalid)

    cascade_df = pd.DataFrame()
    cascade_df['target'] = tgtList_global
    cascade_df['source'] = srcList_global
    cascade_df['rt_time'] = rtList_global
    cascade_df['mid'] = midList_global
    cascade_df['post_time'] = postTimeList_global
    cascade_df['isOn'] = isOn_global

    pickle.dump(cascade_df, open('data/cascade_df_flixster_v2+.pickle', 'wb'))

