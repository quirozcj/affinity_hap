#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy,nbi-largemem,jic-short,nbi-short
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 15G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/scripts/slurmn/combine.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/scripts/slurmn/combine.%N.%j.err # STDERR
#SBATCH --array=0-23
#SBATCH --job-name=combine

source package /nbi/software/testing/bin/python-3.7.2

i=$SLURM_ARRAY_TASK_ID

base_dir=/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq
script_dir=$base_dir/03_haplotype_grouping/scripts
metadata_dir=$base_dir/03_haplotype_grouping/metadata
out_dir=$base_dir/03_haplotype_grouping/combined_by_windows
mkdir -p $out_dir

declare -a REFERENCES=("Arina" "Jagger" "Julius" "Lancer" "Landmark" "Mace" "Norin61" "Spelt" "Stanley" "Mattis" "CS" "Svevo")
declare -a WINDOWS=(100000 50000)

reference=$(($i%12))
reference=${REFERENCES[$reference]}

windows_i=$(($i/12))
window=${WINDOWS[$windows_i]}

python3 $script_dir/combine_by_windows.py \
-i $metadata_dir/samples_id_curated_path.tsv \
-r ${reference} \
-w ${window} \
-o $out_dir/${reference}_vs_all_variations_${window}_windows.tsv
