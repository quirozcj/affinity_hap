
declare -a references=("chinese" "arinalrfor")
declare -a dampings=(0.5 0.6)

chromosome='chr1A'
assembly='chinese'
ibs_window=100000
windows_no=19
sliding_w=3

test_data='./test_data'
conversion_dir=$base_dir/09_mummer_conversions/by_chromosome
out_dir='./out'
mkdir -p $out_dir

python3 affinity_by_windows.py \
-c $test_data/${chromosome}_merged_in_windows__ws-100000_round-3_flank-50000.tsv.gz \
-p $test_data \
-k $chromosome \
-r ${references[@]} \
-a $assembly \
-n $windows_no \
-g ${dampings[@]} \
-w $ibs_window \
-l chr_lenghts.tsv \
-d ${dataset}_HC.tsv \
-s ${sliding_w} \
-o $out_dir/${chromosome}_${assembly}_${windows_no}w_${ibs_window}ibs_sliding${sliding_w}.tsv.gz