#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 5G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/scripts/slurmn/aff.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/scripts/slurmn/aff.%N.%j.err # STDERR
#SBATCH --job-name=aff


base_dir="."
script_dir="$base_dir"
db_dir="$base_dir/test_data"
conversion_dir="$base_dir/test_data"
out_dir="$base_dir/out"
mkdir -p $out_dir

chromosome=chr2B
start=0
end=1000000
ibs_window=100000
windows_no=9
assembly=CS
reference=chinese


python3 $script_dir/affinity_by_windows.py \
-c $conversion_dir/${chromosome}_merged_in_windows__ws-100000_round-3_flank-50000.tsv.gz \
-k $chromosome \
-s $start \
-e $end \
-r $reference \
-a $assembly \
-n $windows_no \
-w $ibs_window \
-d $db_dir/samples_LC.tsv \
-o $out_dir/${chromosome}_${windows_no}w_${ibs_window}bp_${start}_to_${end}_test.tsv \
-p $db_dir
