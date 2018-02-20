import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import datetime
from dateutil.relativedelta import relativedelta

def create_sn_dict():
    links_file = 'flixster/links.txt'
    sn_data = pd.read_csv(links_file, sep="\t", header=None)
    sn_data.columns = ['uid', 'frid']

    # print(sn_data.size)

    sn_dict = {}  # this stores the social network dictionary

    for idx, row in sn_data.iterrows():
        print(idx/sn_data.size*100)

        try:
            sn_dict[row['uid']].append(row['frid'])
        except:
            sn_dict[row['uid']] = []
            sn_dict[row['uid']].append(row['frid'])

    pickle.dump(sn_dict, open('flixster_sn_dict.pickle', 'wb'))

    return sn_dict


def cascade_stats():
    '''

    :return:
    '''
    movie_file = 'ratings_timed/Ratings_timed_v1.txt'
    movie_df = pd.read_csv(movie_file, sep="\s+", encoding = "ISO-8859-1", index_col=False, )
    print(movie_df[:10])
    movie_id_list = list(set(movie_df['movieid'])) # Set of movie ids


    ''' CASCADE SIZE STATS '''
    movie_cacade_size = []
    movie_cascade_size_above_100 = 0
    movie_cascade_size_above_500 = 0
    movie_cascade_size_above_1000 = 0
    movie_cascade_size_above_10000 = 0


    for mid in movie_id_list:
        movie_cascade = movie_df[movie_df['movieid'] == mid]

        movie_cacade_size.append(len(movie_cascade))

        if len(movie_cacade_size) > 100:
            movie_cascade_size_above_100 += 1
        if len(movie_cacade_size) > 500:
            movie_cascade_size_above_500 += 1
        if len(movie_cacade_size) > 1000:
            movie_cascade_size_above_1000 += 1
        if len(movie_cacade_size) > 10000:
            movie_cascade_size_above_10000 += 1

    with open('movie_metadata.txt', 'w') as fw:
        fw.write("Number of movies in data: %d" % len(movie_cacade_size))
        fw.write("\nAverage cascade size: %d" % np.mean(np.array(movie_cacade_size)))
        fw.write("\n# cascades above 100 ratings: %f and perc. of all: %f" %(movie_cascade_size_above_100, movie_cascade_size_above_100/len(movie_cacade_size)))
        fw.write("\n# cascades above 500 ratings: %f and perc. of all: %f" %(movie_cascade_size_above_500, movie_cascade_size_above_500/len(movie_cacade_size)))
        fw.write("\n# cascades above 1000 ratings: %f and perc. of all: %f" %(movie_cascade_size_above_1000, movie_cascade_size_above_1000/len(movie_cacade_size)))
        fw.write("\n# cascades above 10000 ratings: %f and perc. of all: %f" %(movie_cascade_size_above_10000, movie_cascade_size_above_10000/len(movie_cacade_size)))

    # CODE FOR HISTOGRAM PLOT OF MOVIE CASCADE LENGTHS
    hfont = {'fontname': 'Arial'}

    plt.close()
    plt.figure(figsize=(12, 8))
    n, bins, patches = plt.hist(movie_cacade_size, bins=20, facecolor='b')
    plt.xlabel('Cascade size ( # Users)', size=40, **hfont)
    plt.ylabel('', size=40, **hfont)
    plt.title('')
    plt.tick_params(axis='x', labelsize=25)
    plt.tick_params(axis='y', labelsize=25)
    plt.subplots_adjust(left=0.16, bottom=0.16)
    plt.grid(True)
    plt.savefig('Cascade_size_hist.png')
    plt.close()

    ''' CASCADE TEMPORAL STATS '''
    # movie_df['date'] = pd.to_datetime(movie_df['date'])

    time_diff_list = []
    min_time = 0
    max_time = 0

    for mid in movie_id_list:
        movie_cascade = movie_df[movie_df['movieid'] == mid]

        # Temporally order the movie ratings
        movie_cascade = movie_cascade.sort('date')

        first_time = movie_cascade.iloc[0]['date']
        last_time = movie_cascade.iloc[len(movie_cascade)-1]['date']

        if min_time == 0:
            min_time = first_time
            max_time = last_time
        else:
            if first_time < min_time:
                min_time = first_time
            if last_time > max_time:
                max_time = last_time

        first_time = datetime.datetime.strptime(first_time, '%Y-%m-%d')
        last_time = datetime.datetime.strptime(last_time, '%Y-%m-%d')
        diff = last_time - first_time

        time_diff_list.append(diff.days) # in minutes

    with open('movie_metadata.txt', 'a') as fw:
        fw.write("\n Ratings start time: " + min_time )
        fw.write("\n Ratings end time: " + max_time )
        fw.write("\nAverage cascade lifetimee: %f" % np.mean(np.array(time_diff_list)))

    # CODE FOR HISTOGRAM PLOT OF MOVIE CASCADE LENGTHS
    plt.close()
    plt.figure(figsize=(12, 8))
    n, bins, patches = plt.hist(time_diff_list, bins=20, facecolor='b')
    plt.xlabel('Cascade lifetime (days)', size=40, )
    plt.ylabel('', size=40, )
    plt.title('')
    plt.tick_params(axis='x', labelsize=25)
    plt.tick_params(axis='y', labelsize=25)
    plt.subplots_adjust(left=0.16, bottom=0.16)
    plt.grid(True)
    plt.savefig('Cascade_life_hist.png')
    plt.close()

    # Temporal Analysis of ratings by  months
    START_TIME = pd.to_datetime('2006-01-01')
    END_TIME = pd.to_datetime('2009-01-01')

    # FILTER FOR FASTER PROCESSING
    movie_df['date'] = pd.to_datetime(movie_df['date'])

    movie_df = movie_df[movie_df['date'] >= START_TIME]
    movie_df = movie_df[movie_df['date'] < END_TIME]

    start_pos = START_TIME
    end_pos = start_pos + relativedelta(months=1)

    start_dates_list = []
    num_ratings = []
    num_movies = []
    while start_pos <END_TIME:
        movie_df_curr = movie_df[movie_df['date'] >= start_pos]
        movie_df_curr = movie_df_curr[movie_df_curr['date'] < end_pos]

        num_ratings.append(len(movie_df_curr))
        num_movies.append(len(list(set(movie_df_curr['movieid']))))
        start_dates_list.append(start_pos)

        start_pos +=  relativedelta(months=1)
        end_pos += relativedelta(months=1)

    movie_stats_df = pd.DataFrame()
    movie_stats_df['month'] = start_dates_list
    movie_stats_df['number_ratings'] = num_ratings
    movie_stats_df['number_movies'] = num_movies

    return movie_stats_df



def create_cascades():
    '''

    Create the cascade files for all movie ids
    :return:
    '''

    ''' Some preprocessing checks '''

    movie_stats_df = pickle.load(open('movie_stats_df.pickle', 'rb'))
    movie_stats_df['month'] = movie_stats_df['month'].dt.date
    movie_stats_df.plot(x='month', y='number_ratings', lw=3, kind='bar')
    plt.xlabel('Time', size=40, )
    plt.ylabel('# User Ratings', size=40, )
    plt.tick_params(axis='x', labelsize=25)
    plt.tick_params(axis='y', labelsize=25)
    plt.grid(True)
    # plt.subplots_adjust(left=0.16, bottom=0.16)
    plt.show()



def main():
    '''

    :return:
    '''

    ''' This is the social network data file processing  '''
    # sn_dict = create_sn_dict()

    ''' This is for preprocessing the cascade files '''
    # movie_stats_df = cascade_stats()
    # pickle.dump(movie_stats_df, open('movie_stats_df.pickle', 'wb'))

    create_cascades()

if __name__ == "__main__":
    main()