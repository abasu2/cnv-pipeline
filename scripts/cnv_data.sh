#!/bin/bash
trap handler SIGINT

shutdown=false
handler () {
    echo "Caught SIGINT. Finishing current iteration..."
    shutdown=true
}

fasta_file="/share/lanzarolab/users/abasu/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
accessibility_path=$(realpath ../accessibility_masks)
for samples_file in ../sample_names/saotome_samples.txt ../sample_names/principe_samples.txt ; do
    while read -u 5 sample_name; do # read lines from saotome_samples.txt and principe_samples.txt
        mkdir -p ../cnv_results/$sample_name
        # if any of the data files do not exist, call cnv_data on them.
        if [ ! -f ../cnv_results/$sample_name/2R.npz ] || \
           [ ! -f ../cnv_results/$sample_name/2L.npz ] || \
           [ ! -f ../cnv_results/$sample_name/3R.npz ] || \
           [ ! -f ../cnv_results/$sample_name/3L.npz ] || \
           [ ! -f ../cnv_results/$sample_name/X.npz ]
        then
            sample_path="/share/lanzarolab/seq/map/AgamP4/*samples/$sample_name/$sample_name.bam"
            results_path=$(realpath ../cnv_results/$sample_name)
            python cnv_data.py $sample_name $fasta_file $sample_path $accessibility_path $results_path &
            pypid=$!
            wait $pypid # wait for the python process to finish
        fi
        if [ $shutdown = true ] ; then
            echo "Iteration will finish in the background. Exiting."
            exit
        fi
    done 5<$samples_file
done

# python ~/cnv_calling.py $sample_name $sample_path /home/abasu/regions.txt /home/abasu/accessibility_masks --image /home/abasu/cnv_results/$sample_name | tee -a ~/cnv_results/$sample_name/summary.txt ~/all_cns.txt &