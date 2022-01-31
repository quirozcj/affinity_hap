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