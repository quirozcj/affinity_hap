import argparse
import pandas as pd
import numpy as np


def combine_references(references, chromosome, window, path):
	dfs = []
	for reference in references:
		file_db = pd.read_csv(path+f'/{reference}_combined_queries_{window}w.tsv.gz', delimiter='\t')
		file_db = file_db[file_db['seqname'].str.contains(f'{chromosome}')]
		dfs.append(file_db)
	dfs_concat = pd.concat(dfs, join="inner")
	return dfs_concat

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument('-p', '--path')
	parser.add_argument('-w', '--window_size', type=int)
	parser.add_argument('-r', '--references', nargs='+') # this will be a list of references
	parser.add_argument('-k', '--chromosome')
	parser.add_argument('-o', '--output')
	args = parser.parse_args()
	return args

def main():
	args = parse_arguments()

	references = args.references
	path = args.path
	combine_references_df = combine_references(references, args.chromosome, args.window_size, path)
	combine_references_df.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
	main()