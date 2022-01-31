
declare -a references=("chinese" "jagger")
# declare -a dampings=(0.5 0.6 0.7 0.8 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99)
declare -a dampings=(0.5 0.6)

chromosome='chr1A'
assembly='chinese'
ibs_window=50000
windows_no=19
sliding_w=1

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
-s ${sliding_w} \
-o $out_dir/${chromosome}_${assembly}_${windows_no}w_${ibs_window}ibs_sliding${sliding_w}.tsv
