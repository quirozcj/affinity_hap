
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