import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import datetime
from dateutil.relativedelta import relativedelta



def create_cascades():
    '''

    Create the cascade files for all movie ids
    :return:
    '''

    ''' Some preprocessing checks '''

    # movie_stats_df = pickle.load(open('movie_stats_df.pickle', 'rb'))
    # movie_stats_df['month'] = movie_stats_df['month'].dt.date
    # movie_stats_df.plot(x='month', y='number_ratings', lw=3, kind='bar')
    # plt.xlabel('Time', size=40, )
    # plt.ylabel('# User Ratings', size=40, )
    # plt.tick_params(axis='x', labelsize=25)
    # plt.tick_params(axis='y', labelsize=25)
    # plt.grid(True)
    # # plt.subplots_adjust(left=0.16, bottom=0.16)
    # plt.show()

    ''' Create a mapping of cascade source times '''

    print('Create start mappings.....')
    movie_file = 'ratings_timed/Ratings_timed_v1.txt'
    movie_df = pd.read_csv(movie_file, sep="\s+", encoding="ISO-8859-1", index_col=False, )
    movie_df['date'] = pd.to_datetime(movie_df['date'])

    ''' These are the start and end tiems for analysis'''
    START_TIME = pd.to_datetime('2005-03-01')
    END_TIME = pd.to_datetime('2009-01-01')

    movie_df = movie_df[movie_df['date'] >= START_TIME]
    movie_df = movie_df[movie_df['date'] < END_TIME]

    movie_id_list = list(set(movie_df['movieid']))  # Set of movie ids

    ''' The following operations can be shortened by  the group by operation in pandas '''
    start_times_cascade = {}  # MOVIED ---- START DATE mapping
    for mid in movie_id_list:
        movie_cascade = movie_df[movie_df['movieid'] == mid]
        movie_cascade = movie_cascade.sort('date')

        start_times_cascade[mid] = movie_cascade.iloc[0]['date']

    ''' This is the cascade start and end times for consideration '''
    START_TIME_CAS = pd.to_datetime('2006-01-01')
    END_TIME_CAS = pd.to_datetime('2008-01-01')

    cascade_life_sizes = []
    movie_cascade_size_above_100 = 0
    movie_cascade_size_above_500 = 0
    movie_cascade_size_above_1000 = 0
    movie_cascade_size_above_5000 = 0

    movie_id_filter = []
    for mid in start_times_cascade:
        if start_times_cascade[mid] >= START_TIME_CAS and start_times_cascade[mid] < END_TIME_CAS: # checks if the cascades start in this range
            movie_cascade = movie_df[movie_df['movieid'] == mid]
            movie_cascade = movie_cascade[movie_cascade['date'] < END_TIME_CAS]

            # This stores the size of cascades starting in the defined time range
            # and extending till the end time of the range - NOT THE ACTUAL END TIME OF THAT CASCADE
            cascade_life_sizes.append(len(movie_cascade))

            if len(movie_cascade) > 100:
                movie_cascade_size_above_100 += 1
            if len(movie_cascade) > 500:
                movie_cascade_size_above_500 += 1
            if len(movie_cascade) > 1000:
                movie_cascade_size_above_1000 += 1
            else:
                movie_cascade_size_above_5000 += 1

            if len(movie_cascade) >= 200:  # CASCADE MIN SIZE = 300
                movie_id_filter.append(mid)

    # with open('movie_metadata_2007_03.txt', 'w') as fw:
    #     fw.write("Number of movies in data: %d" % len(cascade_life_sizes))
    #     fw.write("\nAverage cascade size: %d" % np.mean(np.array(cascade_life_sizes)))
    #     fw.write("\n# cascades above 100 ratings: %f and perc. of all: %f" % (
    #     movie_cascade_size_above_100, movie_cascade_size_above_100 / len(cascade_life_sizes)))
    #     fw.write("\n# cascades above 500 ratings: %f and perc. of all: %f" % (
    #     movie_cascade_size_above_500, movie_cascade_size_above_500 / len(cascade_life_sizes)))
    #     fw.write("\n# cascades above 1000 ratings: %f and perc. of all: %f" % (
    #     movie_cascade_size_above_1000, movie_cascade_size_above_1000 / len(cascade_life_sizes)))
    #     fw.write("\n# cascades above 5000 ratings: %f and perc. of all: %f" % (
    #     movie_cascade_size_above_5000, movie_cascade_size_above_5000 / len(cascade_life_sizes)))
    #
    # # CODE FOR HISTOGRAM PLOT OF MOVIE CASCADE LENGTHS
    # hfont = {'fontname': 'Arial'}
    #
    # plt.close()
    # plt.figure(figsize=(12, 8))
    # n, bins, patches = plt.hist(cascade_life_sizes, bins=20, facecolor='b')
    # plt.xlabel('Cascade lifetime', size=40, **hfont)
    # plt.ylabel('', size=40, **hfont)
    # plt.title('')
    # plt.tick_params(axis='x', labelsize=25)
    # plt.tick_params(axis='y', labelsize=25)
    # plt.subplots_adjust(left=0.16, bottom=0.16)
    # plt.grid(True)
    # plt.savefig('Cascade_size_hist_2007_03.png')
    # plt.close()

    print('Create the cascades file....')
    ''' Create the cascades rt_df file '''
    movie_file = 'ratings_timed/Ratings_timed_v1.txt'
    movie_df = pd.read_csv(movie_file, sep="\s+", encoding="ISO-8859-1", index_col=False, )
    movie_df['date'] = pd.to_datetime(movie_df['date'])

    sn_dict = pickle.load(open('flixster_sn_dict.pickle', 'rb'))

    ''' Storage structures for cascade dataframe'''
    srcList = []
    tgtList = []
    midList = []
    rtList = []
    postTimeList = []
    isOnList = []

    print('Total number of cascades to consider: ', len(movie_id_filter))
    for mid in movie_id_filter:
        movie_cascade = movie_df[movie_df['movieid'] == mid]
        movie_cascade = movie_cascade[movie_cascade['date'] < END_TIME_CAS]

        movie_cascade = movie_cascade.sort('date')

        ''' Create temporary lists for traversal - THINK OF A BETTER WAY OF DOING THIS '''
        uid_tempList = list(movie_cascade['userid'])
        time_tempList = list(movie_cascade['date'])

        ''' !. Create a link between current node and a future node if the time difference is less than 7 days '''
        ''' 2. If the link between the node and the future node is repeated, ignore it. '''
        for idx_cur in range(len(uid_tempList)):
            curr_time = time_tempList[idx_cur]
            cur_time_str = datetime.datetime.strptime(curr_time.strftime("%Y-%m-%d"), '%Y-%m-%d')

            for idx_link in range(idx_cur+1, len(uid_tempList)):
                fut_time = time_tempList[idx_link]
                fut_time_str = datetime.datetime.strptime(fut_time.strftime("%Y-%m-%d"), '%Y-%m-%d')

                diff = fut_time_str - cur_time_str

                if diff.days > 21:
                    break

                ''' Find the nodes which are friends - connect them'''
                try:
                    friends = sn_dict[uid_tempList[idx_cur]]
                    # print(len(friends))
                except:
                    continue
                if uid_tempList[idx_link] in friends: # if tgt in src's friends list
                    srcList.append(uid_tempList[idx_cur])
                    tgtList.append(uid_tempList[idx_link])

                    rtList.append(time_tempList[idx_link])
                    postTimeList.append(time_tempList[idx_cur])
                    midList.append(mid)
                    isOnList.append(1)

    cascade_df = pd.DataFrame()
    cascade_df['target'] = tgtList
    cascade_df['source'] = srcList
    cascade_df['rt_time'] = rtList
    cascade_df['mid'] = midList
    cascade_df['post_time'] = postTimeList
    cascade_df['isOn'] = isOnList

    pickle.dump(cascade_df, open('cascade_df_06-08_v1+.pickle', 'wb'))


def read_cascades():
    '''
    Test the cascade file to check for consistencies
    :return:
    '''

    cascades_df = pickle.load(open('../flixster_process/data/cascade_df_06-08_v1+.pickle', 'rb'))
    print(cascades_df[:10])

    # print('The number of cascades are: ', len(list(set(cascades_df['mid']))))


def create_network_cascades():
    '''
    Compute the average degree of the networks
    Plot the cascades
    :return:
    '''

    cascades_df = pickle.load(open('../flixster_process/data/cascade_df_06-08_v1+.pickle', 'rb'))

def main():
    '''

    :return:
    '''

    ''' This is the social network data file processing  '''
    # sn_dict = create_sn_dict()

    ''' This is for preprocessing the cascade files '''
    # movie_stats_df = cascade_stats()
    # pickle.dump(movie_stats_df, open('movie_stats_df.pickle', 'wb'))

    ''' This is for creating the cascades file '''
    # create_cascades()

    ''' This is for testing the cascades file '''
    read_cascades()

if __name__ == "__main__":
    main()