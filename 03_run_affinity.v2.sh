#!/bin/bash
#SBATCH --partition=nbi-medium,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 40G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/scripts/slurmn/aff.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/scripts/slurmn/aff.%N.%j.err # STDERR
#SBATCH --job-name=aff

# source package /nbi/software/testing/bin/python-3.7.2

# i=$SLURM_ARRAY_TASK_ID

base_dir=/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq
script_dir=$base_dir/03_haplotype_grouping/JIC_AGIS/scripts
db_dir=$base_dir/03_haplotype_grouping/JIC_AGIS/combined_by_windows_JIC_AGIS_filtered
out_dir=$base_dir/03_haplotype_grouping/JIC_AGIS/affinity
mkdir -p $out_dir

start=6000000
end=6500000
chromosome=chr5B

singularity exec --overlay \
~/overlays/python36_overlay.img \
/software/4493852d-9831-4b90-a949-b7ea7039ce16/4493852d-9831-4b90-a949-b7ea7039ce16.img \
python3 $script_dir/affinity_by_windows_v3.py \
-c $script_dir/chr5B_merged_in_windows__ws-100000_round-3_flank-50000.tsv \
-k $chromosome \
-s $start \
-e $end \
-r chinese \
-a CS \
-n 19 \
-w 100000 \
-d samples_LC.tsv \
-o $out_dir/${chromosome}_19w_100kb_${start}_to_${end}_cleaned_bash.tsv \
-p /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/combined_by_windows