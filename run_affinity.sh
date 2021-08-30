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

chromosome=chr2B
start=500000000
end=600000000
ibs_window=100000

singularity exec --overlay \
~/overlays/python36_overlay.img \
/software/4493852d-9831-4b90-a949-b7ea7039ce16/4493852d-9831-4b90-a949-b7ea7039ce16.img \
python3 $script_dir/affinity_by_windows_v4.py \
-c $script_dir/${chromosome}_merged_in_windows__ws-100000_round-3_flank-50000.tsv \
-k $chromosome \
-s $start \
-e $end \
-r chinese \
-a CS \
-n 19 \
-w $ibs_window \
-d samples_LC.tsv \
-o $out_dir/${chromosome}_19w_${ibs_window}kb_${start}_to_${end}.tsv \
-p /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/combined_by_windows_JIC_AGIS_filtered
# -p /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotype_grouping/JIC_AGIS/combined_by_windows_JIC_AGIS_filtered