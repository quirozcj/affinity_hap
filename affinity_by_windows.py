import pandas as pd
import numpy as np
import argparse
from sklearn.cluster import AffinityPropagation
from datetime import datetime
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from itertools import count

startTime = datetime.now()

<<<<<<< HEAD
def get_target_region(db_file, reference, down_region, up_region):
    target_region_df = db_file[(db_file['start'].between(down_region, up_region, inclusive='both')
        & (db_file['reference'] == reference))]
=======

def get_target_region(by_windowd_db, reference, chromosome, down_region, up_region, function=None, power_n=None):
    
    if function is not None:
        function_name = function.__name__
    function_name = "None"
    
    in_db = by_windowd_db
    target_region_df = in_db[(in_db['window'] >= down_region) &
                             (in_db['window'] <= up_region) &
                             (in_db['seqname'].str.contains(chromosome)) &
                             (in_db['reference'] == reference)
                             ]
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
    target_region_df.drop(columns=['reference'],inplace=True)
    return target_region_df

<<<<<<< HEAD
def cluster_by_haplotype (df, dampings):
    t_df = df.set_index(['seqname', 'start', 'end']).T
    t_df = t_df.loc[:,~t_df.columns.duplicated()]
    X = StandardScaler(with_mean=True).fit_transform(t_df)
    for dmp in dampings:
        af = AffinityPropagation(random_state=0, damping=dmp).fit(X)
=======
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


def transform_data (function):
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
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
        X_predicted = af.fit_predict(X)
        if len(af.cluster_centers_indices_) > 2:
            sc_score = '%0.2f'%metrics.silhouette_score(X, af.labels_, metric='sqeuclidean')
        else:
            sc_score = 0
        t_df.insert(loc=0,  column=f'sc_{dmp}', value=sc_score, allow_duplicates=True)
        t_df.insert(loc=0,  column=f'dmp_{dmp}', value=X_predicted, allow_duplicates=True)
    return t_df

def header_prefix(items_list, prefix):
    prefix_ = prefix
    items_list = [prefix_ + str(x) for x in items_list]
    return items_list

def drop_samples(samples):
    samples_LC = pd.read_csv(samples, delimiter='\t')
    samples_l = ['seqname', 'start', 'end']
    samples_l += samples_LC['genotype'].to_list()
    return samples_l

<<<<<<< HEAD
def assembly_reference_filter(df, assembly):
	filter_df = df[df['assembly'].isin([assembly])]
	return filter_df

def build_assembly_df(start_region, end_region, window, chromosome_sfx, assembly):
	blocks_list=[]
	for i in count(start_region, window):
		if i > end_region:
			break
		else:
			blocks_list.append(chromosome_sfx+':'+str(start_region+1)+'-'+str(start_region+100000))
			start_region += window
	blocks_df = pd.DataFrame(blocks_list, columns=['block_no'])
	blocks_df['assembly'] = assembly
	blocks_df['reference'] = assembly
	blocks_df['chromosome'] = chromosome_sfx
	blocks_df['start'] = blocks_df['block_no'].str.split(r':|-', expand=True)[1].astype(int)
	blocks_df['end'] = blocks_df['block_no'].str.split(r':|-', expand=True)[2].astype(int)
	blocks_df['orientation'] = '+'
	return blocks_df	

def blocks_in_region(conversions_file, start_region, end_region, assembly):
	region_df = assembly_reference_filter(pd.read_csv(conversions_file, delimiter='\t'), assembly)
	chromosome_sfx = region_df.iloc[0]['block_no'].split(':')[0]
	window = 50000
	assembly_df = build_assembly_df(start_region, end_region, window, chromosome_sfx, assembly)
	blocks = assembly_df['block_no'].unique()
	region_df = pd.concat([region_df, assembly_df], join="inner", ignore_index=True)
	return blocks, region_df
    
def regions_in_blocks(file_df, blocks, chromosome, assembly):
    chromosome_sfx = file_df.iloc[0]['block_no'].split(':')[0]
    target_region_df = file_df[file_df['block_no'].isin(blocks)]
    target_region_df.loc[target_region_df.chromosome == chromosome_sfx,'reference'] = assembly
    return target_region_df

def variations_by_windows(blocks, target_region_df, references, chromosome, num_of_windows, db_file, dampings):
    start_l = 0
    region_df = []
=======
# delete hard coded by a file
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


def blocks_in_region(conversions_file, chromosome, start_region, end_region, reference):
    chromosome_sfx = chromosome+'__chi'
    file_df =pd.read_csv(conversions_file, delimiter='\t')
    region_df = \
    file_df[\
           (file_df['start'] >= start_region) & 
           (file_df['end'] <= end_region) & 
           (file_df['chromosome'] == chromosome_sfx) & 
           (file_df['assembly'] == reference) 
            ]

    blocks = region_df.sort_values(by=['end','start'],ascending=True)['block_no'].unique()
    return blocks, file_df

def add_region_size(target_region_df):
    target_region_df['mapped_size'] = target_region_df['end'] - target_region_df['start'] 

def regions_in_blocks(file_df, blocks, chromosome, reference, assembly):
    chromosome_sfx = chromosome+'__chi'
    target_region_df = file_df[file_df['block_no'].isin(blocks)]
    target_region_df['reference'] = target_region_df['reference'].map(reference_rename)
    target_region_df['assembly']  = target_region_df['assembly'].map(reference_rename)
    target_region_df.loc[target_region_df.chromosome == chromosome_sfx,'reference'] = assembly
    add_region_size(target_region_df)

    #Test to remove this part. Remove training "\"
    target_region_df['end'] = np.where(
                        (target_region_df['mapped_size'] > 50000) & 
                        (target_region_df['mapped_size'] <= 100000), 
                        target_region_df['start'] + 100000, 
                        target_region_df['end'])
    target_region_df['end'] = np.where(
                        target_region_df['mapped_size'] <= 50000, 
                        target_region_df['start'] + 50000, 
                        target_region_df['end'])
    target_region_df['end'] = np.where(
                        target_region_df['mapped_size'] > 100000, 
                        (np.around(target_region_df['mapped_size'], decimals = -5) + target_region_df['start']), 
                        target_region_df['end'])

    assembly_start = target_region_df['block_no'].str.split(r':|-', expand=True)[1]
    assembly_end   = target_region_df['block_no'].str.split(r':|-', expand=True)[2]
    target_region_df.loc[target_region_df.assembly == target_region_df.reference,'start'] = assembly_start.astype(int)
    target_region_df.loc[target_region_df.assembly == target_region_df.reference,'end'] = assembly_end.astype(int)

    add_region_size(target_region_df)
        
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    return target_region_df


def variations_by_windows(blocks, target_region_df, references, chromosome, num_of_windows, db_file, dampings):
    
    # references = ['CS','Jagger','Arina']
    start_l = 0
    region_df = []
    print(references)
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
    while start_l < len(blocks):
        block_no = blocks[start_l:(start_l + num_of_windows)]
        dfs_block = []
        for block in block_no:
<<<<<<< HEAD
            region_in_blocks = target_region_df[target_region_df['block_no'] == block]
            dfs = []
            extract_position_for_reference(references, chromosome, db_file, region_in_blocks, dfs)
=======
            region_in_blocks = \
                (target_region_df
                    [target_region_df
                    ['block_no'] == block]
                )
            dfs = []
            extract_position_for_reference(references, chromosome, db_file, region_in_blocks, dfs)
            
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
            dfs_concat = pd.concat(dfs, join="inner")
            dfs_block.append(dfs_concat)
        dfs_block_concat = pd.concat(dfs_block, join="inner")
        dfs_block_concat = dfs_block_concat.sort_values(by=['seqname','start','end'], ascending=True)
        affinity_group = cluster_by_haplotype(dfs_block_concat, dampings)
        dmp_positions = (len(dampings) * 2)
        values = affinity_group.columns.get_level_values(0)[:dmp_positions].values
        affinity_group.columns = ['_'.join(tuple(map(str, t))) for t in affinity_group.columns.values]
        keys = affinity_group.columns.values[:dmp_positions]
        dictionary = dict(zip(keys, values))
        affinity_group.rename(columns=dictionary, inplace = True)
        affinity_group.index.names = ['query']
        affinity_group['block_start'], affinity_group['block_end'] = block_no[0],block_no[-1]
        affinity_group.reset_index(level=0, inplace=True)
        melt_list = affinity_group.columns.values[:dmp_positions+1].tolist() + affinity_group.columns.values[-2:].tolist()
        affinity_group = pd.melt(affinity_group,
                                id_vars = melt_list,
                                var_name = 'window',
                                value_name = 'variations',
                                value_vars = affinity_group.iloc[:, (dmp_positions):-2])
        start_l += num_of_windows + 1
        region_df.append(affinity_group)
    region_haplotype = pd.concat(region_df, axis=0, ignore_index=False)
    header = region_haplotype.columns[1:dmp_positions+1].values.tolist()
    region_haplotype['chr'] = region_haplotype['block_start'].str.split(r':|-', expand=True)[0]
    region_haplotype['start'] = region_haplotype['block_start'].str.split(r':|-', expand=True)[1].astype(int)
    region_haplotype['end'] = region_haplotype['block_end'].str.split(r':|-', expand=True)[2].astype(int)
    header = ['chr','start','end','query','window','variations'] + header
    region_haplotype = region_haplotype[header]
    region_haplotype = region_haplotype.sort_values(by=['query','start'], ascending=True)
    return region_haplotype

def extract_position_for_reference(references, chromosome, db_file, region_in_blocks, dfs):
    for reference in references:
<<<<<<< HEAD
        if region_in_blocks [region_in_blocks['reference'] == reference].empty:
            continue        
        
        region_by_ref = region_in_blocks[region_in_blocks['reference'] == reference]
        for index, row in region_by_ref.iterrows():
            if row['end'] - row['start'] < 50000:
                down_region = row['start'] - 50000
                up_region = row['end']
            else:
	            down_region = row['start']
	            up_region = row['end']

            extracted_region = get_target_region(db_file, reference, down_region, up_region)
            dfs.append(extracted_region)

=======
        if region_in_blocks [region_in_blocks
                        ['reference'] == reference].empty:
            continue
        
        region_in_blocks_ = region_in_blocks[region_in_blocks['reference'] == reference]
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

        extracted_region = get_target_region(db_file, reference, chromosome, down_region, up_region)
        dfs.append(extracted_region)
        

# merge with name conversion sabove
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
references_names = {\
        'WhJag':'jagger',
        'WhAri':'arinalrfor',
        'Whjul':'julius',
        'Whlan':'lancer',
        'WhLan':'landmark',
        'Whmac':'mace',
        'WhSYM':'sy_mattis',
        'WhNor':'norin61',
        'Whspe':'spelta',
        'WhSta':'stanley'}

def map_substring(s, dict_map):
    for key in dict_map.keys():
        if key in s: return dict_map[key]
    return 'chinese'

def get_highest_dmp(df, dmps):
    sc_list=[]
    for dmp in dmps:
        sc_list.append(f'sc_{dmp}')
    df['sc_max'] = df[sc_list].astype(float).idxmax(axis=1)
    df['sc_max'] = df['sc_max'].str.split('_', expand=True)[1]
    df['sc_max'] = 'dmp_' + df['sc_max'].astype(str)

    dmps_list=[]
    for dmp in dmps:
        dmps_list.append(f'dmp_{dmp}')
    idx = len(df.columns)-1
    df.insert(loc=idx,  column='sc_val', value=0)
    df.insert(loc=idx,  column='dmp_max', value=0)
    for dmp in dmps_list:
        df['dmp_max'] = np.where(df['sc_max'] == dmp, df[f'{dmp}'], df['dmp_max'])
        sc_dmp = dmp.replace('dmp_', 'sc_')
        df['sc_val'] = np.where(df['sc_max'] == dmp, df[f'{sc_dmp}'], df['sc_val'])
    return df

def chr_region(chr_lengths, assembly, chromosome):
    len_df = pd.read_csv(chr_lengths, delimiter='\t')
    ref_df = len_df[len_df['assembly'] == assembly]
    chr_df = ref_df[ref_df['chr'].str.contains(f'{chromosome}')]
    chr_start = chr_df.iloc[0]['start'].astype(int)
    chr_end = chr_df.iloc[0]['end'].astype(int)
    return chr_start, chr_end


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path_to_db')
    parser.add_argument('-c', '--conversions_file')
    parser.add_argument('-k', '--chromosome')
    parser.add_argument('-r', '--references', nargs='+')
    parser.add_argument('-a', '--assembly')
    parser.add_argument('-n', '--num_of_windows', type=int)
    parser.add_argument('-w', '--window', type=int)
    parser.add_argument('-d', '--drop_samples_f', default = None)
    parser.add_argument('-o', '--output')
    parser.add_argument('-l', '--chr_lengths')
    parser.add_argument('-g', '--dampings', nargs='+')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    chromosome = args.chromosome
    reference = args.assembly
    assembly = args.assembly
    window = args.window
    num_of_windows = args.num_of_windows
    conversions_file = args.conversions_file
    start_region, end_region = chr_region(args.chr_lengths, assembly, chromosome)

<<<<<<< HEAD
    blocks, file_df = blocks_in_region(conversions_file, start_region, end_region, assembly)
    target_region_df = regions_in_blocks(file_df, blocks, chromosome, assembly)
    target_blocks = target_region_df[target_region_df['reference'] == assembly].sort_values(by=['start'],ascending=True)
=======
    blocks, file_df = blocks_in_region(conversions_file, chromosome, start_region, end_region, reference)
    target_region_df = regions_in_blocks(file_df, blocks, chromosome, reference, assembly)
    target_blocks = target_region_df[target_region_df['reference']== assembly].sort_values(by=['start'],ascending=True)
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
    blocks = target_blocks['block_no'].unique()
    references = args.references
    path_to_db = args.path_to_db
    samples = drop_samples(args.drop_samples_f)
    dampings = list(map(float, args.dampings))
    damping_list = header_prefix(dampings, 'dmp_')
    sc_list = header_prefix(dampings, 'sc_')

    db_file = pd.read_csv(f'{path_to_db}/{chromosome}_variations_{window}w.tsv', delimiter='\t')
    db_file = db_file.filter(items=samples)

    idx = 2
    new_col = db_file['seqname'].apply(lambda x: map_substring(x, references_names))
    db_file.insert(loc=idx, column='reference', value=new_col)

    region_haplotype = variations_by_windows(blocks, target_region_df, references, chromosome, num_of_windows, db_file, dampings)

    Y = StandardScaler(with_mean=True).fit_transform(region_haplotype['variations'].to_numpy().reshape(-1, 1))
    region_haplotype.insert(loc=6,  column='variations_scl', value=Y)
    region_haplotype = get_highest_dmp(region_haplotype, dampings)
    region_haplotype = region_haplotype.drop(damping_list + sc_list, axis=1)

    grouped = region_haplotype.groupby(['start'])
    group_size = region_haplotype.groupby(['start','query'])

    region_haplotype['w_num'] = group_size['window'].transform(len)
    region_haplotype['dmp_num'] = grouped['dmp_max'].transform('max')+1
    region_haplotype['std'] = grouped['variations'].transform('std').round(1)
    region_haplotype['median'] = grouped['variations'].transform('median')
    region_haplotype['mean'] = grouped['variations'].transform('mean').round(1)
    # region_haplotype['variance'] = grouped['variations'].transform('var').round(1)
    region_haplotype['skew'] = grouped['variations'].transform('skew').round(2)
    # region_haplotype['kurt'] = grouped['variations'].transform(pd.DataFrame.kurt).round(2)

<<<<<<< HEAD
=======
    region_haplotype = variations_by_windows(blocks, target_region_df, references, chromosome, num_of_windows, db_file, dampings)
>>>>>>> b0b66ca6763e9a6cdfec35035e44ad2fce75bd22
    region_haplotype.to_csv(args.output, sep='\t', index=False)
    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)

    # print(region_haplotype)
    print(datetime.now() - startTime)

if __name__ == '__main__':
    main()
