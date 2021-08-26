import argparse
import pandas as pd

def combine_db(samples_id, reference):
    samples_id = pd.read_csv(samples_id, delimiter='\t')
    samples_by_reference = samples_id[samples_id['reference'] == reference]
    
    samples_combined_db = pd.DataFrame()
    
    for index, row in samples_by_reference.iterrows():
        variations_db = pd.read_csv(row['path']+row['file'], delimiter='\t')
        sample_name = row['query']
        samples_combined_db[sample_name] = variations_db['variations']
    row_names = variations_db[['seqname', 'start', 'end']]
    samples_combined_db = pd.concat([row_names, samples_combined_db], axis=1)
    return samples_combined_db

def count_by_windows(combined_samples, window_size):
    in_db = combined_samples
    # in_db = pd.read_csv(combined_samples, delimiter='\t')
    window_size = window_size
    chrLen = in_db['end'].max() # longest chromosome
    
    w_pos = 0
    db_byChr = pd.DataFrame()
    
    while w_pos <= chrLen:
        by_windows_df = in_db[(in_db['start'] >= w_pos) & (in_db['start'] < w_pos + window_size)]
        by_windows_df = by_windows_df.drop(['start','end'], axis=1)
        by_windows_df['window'] = w_pos + window_size
        w_pos += window_size
        db_byChr = db_byChr.append(by_windows_df)  
        
    by_windows_db = db_byChr.groupby(['seqname','window']).sum().reset_index()
    return by_windows_db

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--samples_id')
    parser.add_argument('-w', '--window_size', type=int)
    parser.add_argument('-r', '--reference')
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    variations_combined = combine_db(args.samples_id, args.reference)
    by_windowd_db = count_by_windows(variations_combined, args.window_size)
    by_windowd_db.to_csv(args.output, sep='\t', index=False)
    # variations_combined.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()