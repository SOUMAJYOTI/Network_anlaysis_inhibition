import pickle
import numpy as np
import pandas as pd


def main():
    print('Loading diffusion file...')
    diff_file = 'cascade_df_06-08_v1+.pickle'
    df = pd.read_pickle(diff_file)
    # df = pd.read_csv(diff_file, names=['index', 'target', 'source', 'rt_time', 'mid', 'post_time', 'isOn'])
    df['mid'] = df['mid'].astype('str')
    steep_times, inhib_times = pickle.load(open('steep_inhib_times_dict_flixster.pickle', 'rb'))

    ''' Sample of cascades for testing '''
    mid_count = 0
    cascades_list_sample = pd.DataFrame()
    for mid in steep_times:
        cascade = df[df['mid'] == str(mid)]
        cascades_list_sample = pd.concat([cascades_list_sample, cascade])

        mid_count +=1
        if mid_count > 10:
            break

    ''' Cascades with inhib times '''
    mid_count = 0
    cascades_list = pd.DataFrame()
    for mid in steep_times:
        cascade = df[df['mid'] == str(mid)]
        cascades_list = pd.concat([cascades_list, cascade])

        mid_count += 1
        # if mid_count > 500:
        #     break

    pickle.dump(cascades_list_sample, open('cascades_flixster_sample.pickle', 'wb'))
    pickle.dump(cascades_list, open('cascades_flixster_train.pickle', 'wb'))

if __name__ == "__main__":
    main()