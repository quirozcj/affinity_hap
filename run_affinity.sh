#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 100G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/ibspy_analysis/scripts/log/aff_mt.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/ibspy_analysis/scripts/log/aff_mt.%N.%j.err # STDERR
#SBATCH --job-name=aff_mt
#SBATCH --array=0-20

i=$SLURM_ARRAY_TASK_ID

# chromosome list
declare -a chromosomes=()
for n in {1..7..1}; do
	for l in A B D; do
		chromosomes+=(chr${n}${l});
	done;
done

declare -a references=("chinese" "arinalrfor" "sy_mattis" "jagger" "stanley" "julius" "lancer" "landmark" "mace" "norin61" "spelta")
declare -a assemblies=("arinalrfor")

chr_i=$(($i%21))
chromosome=${chromosomes[$chr_i]}

ref_i=$(($i/21))
assembly=${assemblies[$ref_i]}

declare -a dampings=(0.5 0.6 0.7 0.8 0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99)
# declare -a dampings=(0.5 0.6)

ibs_window=50000
windows_no=19
score=variations_modern_tauschii

base_dir=/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq
script_dir=$base_dir/03_haplotypes/ibspy_analysis/scripts
db_dir=$base_dir/03_haplotypes/ibspy_analysis/ibspy_combined/variations_tauschii/references_combined
conversion_dir=$base_dir/09_mummer_conversions/by_chromosome
out_dir=/nbi/group-data/ifs/NBI/Cristobal-Uauy/Jesus/14_watseq_affinity/affinity_long_files/modern_tauschii/${assembly}/${score}
mkdir -p $out_dir

python3 $script_dir/affinity_by_windows.py \
-c $conversion_dir/${chromosome}_merged_in_windows__ws-100000_round-3_flank-50000.tsv.gz \
-p $db_dir \
-k $chromosome \
-r ${references[@]} \
-a $assembly \
-n $windows_no \
-g ${dampings[@]} \
-w $ibs_window \
-l chr_lenghts.tsv \
-d modern_tauschii_HC.tsv \
-o $out_dir/${chromosome}_${assembly}_${windows_no}w_${ibs_window}ibs.tsv.gz
