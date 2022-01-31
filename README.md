# Automatically call haplotypes by k-mers


### concatenate multiple outpus of IBSpy into one file by reference
- use the script ```run_combine.sh```. It requres a tab separated imput file (-i) argument with the name of "reference", "query", "file", & the patch to the individual files. See exaple (samples_metadata.tsv). (-r), reference, (-w), window size to combine. This can be > than 50,000 bp of the output of IBSpy. (-s), the score to extract form the original IBSpy output file. (-o), output file name (use .gz to compress).

```sh

test_data='./test_data'
out_dir='./out'
mkdir -p $out_dir

# reference="jagger"
reference="chinese"
window=50000
score='variations'

python3 combine_by_windows.py \
-i $test_data/samples_metadata.tsv \
-r ${reference} \
-w ${window} \
-s ${score} \
-o $out_dir/${reference}_combined_queries_${window}w.tsv.gz
```


### Merge all variations counts and query samples in a single file by chromosome
- use the script ```run_merge.sh```. (-p), path to the files to merge. (-w), windows size, this must match the windows used in the step before (combine_by_widnows). If a query sample is not present in one of the reference concatenated files, this will be dicarded. (-r) the references to include in the file. You can add as many references you want. (-k), chromosomes to merges. One file per chromosome will be created  (-o).

```sh
declare -a references=("chinese" "jagger")
chromosome='chr1A'
window=50000

db_dir='./out'
out_dir='./out'
mkdir -p $out_dir

python3 merge_references.py \
-p $db_dir \
-w ${window} \
-r ${references[@]} \
-k ${chromosome} \
-o $out_dir/${chromosome}_variations_${window}w.tsv
```

### Call haplotypes by windows using Affinity Propagation (AP) algorithm (https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AffinityPropagation.html)
- use the script ```run_affinity.sh```. This step uses the file of the previous step and the file with the corresponding windows physical position in each pengenome assembly. (-c), file with sinteny regions of each assembly. (-p), path to the database file with the variaitons counts crated in the step before. (-k) chromosome. (-r), references to include in the analysisis. (-a), assembly to use as a template, this will be used to map all the correspondin windows from the other pangenome assmblies. Any assembly can be used as a template. (-n), number of windows, this is the number of sub-windows to include in a block of windows; 19 windows corresponds (50,000 bp sub-windows) a 1 Mbp block of windows. (-g), damping parameter to tunne for the AP clustering algorithm. (-w), window size used when concatenated multiple individual files from the IBSpy outputs. (-l), chromosome lenghts of the genome asembly used as a template. (-d), genotypes to inlcude in the analysis. (-o), output file, compress by adding .gz.

```sh

declare -a references=("chinese" "jagger")
# declare -a dampings=(0.5 0.6 0.7 0.8 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99)
declare -a dampings=(0.5 0.6)

chromosome='chr1A'
assembly='chinese'
ibs_window=50000
windows_no=19
# sliding_w=3

test_data='./test_data'
db_dir='./out'
test_data='./test_data'
out_dir='./out'
mkdir -p $out_dir

python3 affinity_by_windows.py \
-c $test_data/${chromosome}_merged_in_windows__ws-100000_round-3_flank-50000.tsv.gz \
-p $db_dir \
-k $chromosome \
-r ${references[@]} \
-a $assembly \
-n $windows_no \
-g ${dampings[@]} \
-w $ibs_window \
-l $test_data/chr_lenghts_test.tsv \
-d $test_data/watseq_HC_test.tsv \
-o $out_dir/${chromosome}_${assembly}_${windows_no}w_${ibs_window}ibs_sliding${sliding_w}.tsv
# -s ${sliding_w} \
```