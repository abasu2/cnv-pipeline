#!bin/bash

sample_name=$1
sample_path="/share/lanzarolab/seq/map/AgamP4/*samples/$sample_name/$sample_name.bam"
fasta_file="/share/lanzarolab/users/abasu/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
python ~/cnv_data_fixed.py $sample_name $fasta_file $sample_path /home/abasu/accessibility_masks /home/abasu/cnv_results/$sample_name