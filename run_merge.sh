
declare -a references=("chinese" "arinalrfor")

chromosome='chr1A'
window=50000

python3 merge_references.py \
-p ./test_data \
-w ${window} \
-r ${references[@]} \
-k ${chromosome} \
-o $out_dir/${chromosome}_variations_${window}w.tsv