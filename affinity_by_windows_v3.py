import pandas as pd
import numpy as np
import argparse
from sklearn.cluster import AffinityPropagation

from warnings import filterwarnings
filterwarnings('ignore')


def get_target_region(by_windowd_db, chromosome, down_region, up_region, function=None, power_n=None):
    
    if function is not None:
        function_name = function.__name__
    function_name = "None"
    
    in_db = by_windowd_db
    target_region_df = in_db[\
                             (in_db['window'] >= down_region) & \
                             (in_db['window'] <= up_region) & \
                             (in_db['seqname'].str.contains(chromosome))\
                            ]
    header = target_region_df[['seqname','window']]
    df_values = target_region_df.iloc[:, 2:]
    
    if function_name == 'power':
        power_n = power_n        
        np_values = function(df_values, power_n)
        final_df = pd.concat([header, np_values], axis=1)
    
    elif function is not None:
        np_values = function(df_values).replace(-np.inf, 0)
        final_df = pd.concat([header, np_values], axis=1)
    else:
        final_df = target_region_df
    return final_df


def transform_data (function, power_n = None):
    
    if function is None:
        function_name = "None"
    else:
        function_name = function.__name__
    return function_name, function


def cluster_by_haplotype (df, damping):
    
    df_ = df.set_index(['seqname', 'window']).T
    df_ = df_.loc[:,~df_.columns.duplicated()]
    X = df_.to_numpy()

    dfs_ = []
    for damping_ in damping:
        af = AffinityPropagation(random_state=0, verbose=False, damping=damping_).fit(X)
        X_predicted = af.fit_predict(X)
        df_.insert(loc=0,  column=f'dmp_{damping_}', value=X_predicted, allow_duplicates=True)
    return df_

def header_prefix(items_list, prefix):
    prefix_ = prefix
    items_list = [prefix_ + str(x) for x in items_list]
    return items_list


def drop_samples(samples):
    samples_LC = pd.read_csv(samples, delimiter='\t')
    samples_l = samples_LC['genotype'].to_list()
    return samples_l

reference_rename = {\
'arinalrfor': 'Arina',
'jagger':'Jagger',
'stanley':'Stanley',
'julius':'Julius',
'lancer':'Lancer',
'landmark':'Landmark',
'mace':'Mace',
'norin61':'Norin61',
'chinese':'CS',
'sy_mattis':'Mattis',
'spelta':'Spelt'}

# data_path = '/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/combined_by_windows_JIC_AGIS_filtered'
# data_path = '/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/combined_by_windows'

# # def conversions_file(conversions_file,samples,chromosome,start_region,end_region,reference,assembly,window,num_of_windows,path_to_db):

def blocks_in_region(conversions_file, chromosome, start_region, end_region, reference, assembly):

    chromosome_sfx = chromosome+'__chi'

    file_df = \
        (pd.read_csv(conversions_file,
              delimiter='\t')
        )

    region_df = \
    file_df[\
           (file_df['start'] >= start_region) & \
           (file_df['end'] <= end_region) & \
           (file_df['chromosome'] == chromosome_sfx) & \
           (file_df['assembly'] == reference) \
            ]

    blocks = (region_df
            .sort_values(by=['start'],ascending=True)
            ['block_no'].unique()
            )

    return blocks, file_df


def regions_in_blocks(file_df, blocks, chromosome, reference, assembly):

    chromosome_sfx = chromosome+'__chi'

    target_region_df = \
        (file_df[file_df['block_no']
        .isin(blocks)]
        )

    target_region_df['reference'] = \
        (target_region_df['reference']
        .map(reference_rename)
        )

    target_region_df['assembly'] =  \
        (target_region_df['assembly']
        .map(reference_rename)
        )

    target_region_df.loc[target_region_df.chromosome == chromosome_sfx,'reference'] = assembly

    target_region_df['mapped_size'] = \
        (target_region_df['end'] - 
        target_region_df['start']
        )

    target_region_df['end'] = \
                        np.where(\
                        (target_region_df['mapped_size'] > 50000) & \
                        (target_region_df['mapped_size'] <= 100000), \
                        target_region_df['start'] + 100000, \
                        target_region_df['end'])

    target_region_df['end'] = \
                        np.where(\
                        target_region_df['mapped_size'] <= 50000, \
                        target_region_df['start'] + 50000, \
                        target_region_df['end'])

    target_region_df['end'] = \
                        np.where(\
                        target_region_df['mapped_size'] > 100000, \
                        (np.around(target_region_df['mapped_size'], decimals = -5) + target_region_df['start']), \
                        target_region_df['end'])

    target_region_df['mapped_size'] = \
            (target_region_df['end'] -
            target_region_df['start']
            )

    return target_region_df


def variations_by_windows(blocks, target_region_df, references, samples, chromosome, reference, window, num_of_windows, path_to_db, dampings,damping_list):
    
    # references = ['CS','Jagger','Arina']

    start_l = 0
    region_df = []

    while start_l <= len(blocks):
        
        block_no = blocks[start_l:(start_l + num_of_windows)]

        dfs_block = []
        for block in block_no:

            region_in_blocks = \
                (target_region_df
                    [target_region_df
                    ['block_no'] == block]
                )

            dfs = []
            for reference in references:

                if region_in_blocks [region_in_blocks
                        ['reference'] == reference].empty:
                    continue

                else:
                    region_in_blocks_ = (region_in_blocks
                        [region_in_blocks['reference'] == reference]
                    )

                if len(region_in_blocks_) > 1:

                    down_region = (region_in_blocks_
                        .groupby(by=['reference'])
                        .agg({'start':'min'})
                        .reset_index()['start'].values[0]
                        )

                    up_region = (region_in_blocks_
                        .groupby(by=['reference'])
                        .agg({'end':'max'})
                        .reset_index()['end'].values[0]
                    )

                else:
                    down_region = \
                    (int(region_in_blocks_[region_in_blocks_
                        ['reference'] == reference]['start'].values) 
                    )
                    up_region = \
                    (int(region_in_blocks_[region_in_blocks_
                        ['reference'] == reference]['end'].values)
                    )


                file_db = \
                pd.read_csv(f'{path_to_db}/{reference}_vs_all_variations_{window}_windows.tsv',
                    delimiter='\t')

                file_db = (file_db
                    .drop(samples,axis=1,errors='ignore')
                    )
                
                extracted_region = get_target_region(file_db,chromosome,down_region,up_region)

                dfs.append(extracted_region)

            dfs_concat = pd.concat(dfs, join="inner")
            dfs_block.append(dfs_concat)

        dfs_block_concat = pd.concat(dfs_block, join="inner")

        affinity_group = cluster_by_haplotype(dfs_block_concat, dampings)

        affinity_group['block_no_start'], affinity_group['block_no_end'] = block_no.min(),block_no.max()
        
        
        region_df.append(affinity_group)

        start_l += num_of_windows + 1
        # print(region_df)

    region_haplotype = pd.concat(region_df, axis=1, ignore_index=False)
    # print(region_haplotype)
    

    region_haplotype.index.names = ['query']
    region_haplotype.columns = ['_'.join(tuple(map(str, t))) for t in region_haplotype.columns.values]
    region_haplotype.rename(columns={\
                                     "dmp_0.9_": "dmp_0.9", \
                                     "dmp_0.8_": "dmp_0.8", \
                                     "dmp_0.7_": "dmp_0.7", \
                                     "dmp_0.6_": "dmp_0.6",\
                                     "dmp_0.5_": "dmp_0.5", \
                                     "block_no_end_": "block_no_end", \
                                     "block_no_start_": "block_no_start"}, inplace = True)
    return region_haplotype


def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p', '--path_to_db')
    parser.add_argument('-c', '--conversions_file')
    parser.add_argument('-k', '--chromosome')
    parser.add_argument('-s', '--start_region', type=int)  
    parser.add_argument('-e', '--end_region', type=int)
    parser.add_argument('-r', '--reference')
    parser.add_argument('-a', '--assembly')
    parser.add_argument('-n', '--num_of_windows', type=int)
    parser.add_argument('-w', '--window', type=int)
    parser.add_argument('-d', '--drop_samples_f')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    chromosome = args.chromosome
    start_region = args.start_region
    end_region = args.end_region
    reference = args.reference
    assembly = args.assembly
    window = args.window
    conversions_file = args.conversions_file

    blocks, file_df = blocks_in_region(conversions_file, chromosome, start_region, end_region, reference, assembly)
    target_region_df = regions_in_blocks(file_df, blocks, chromosome, reference, assembly)
    references = target_region_df.reference.unique()

    num_of_windows = args.num_of_windows
    path_to_db = args.path_to_db
    samples = drop_samples(args.drop_samples_f)

    dampings = [0.5,0.6,0.7,0.8,0.9]
    damping_list = header_prefix(dampings, 'dmp_')

    region_haplotype = variations_by_windows(blocks, target_region_df, references, samples, chromosome, reference, window, num_of_windows, path_to_db, dampings,damping_list)
    print(region_haplotype)
    region_haplotype.to_csv(args.output, sep='\t', index=True)

if __name__ == '__main__':
    main()
