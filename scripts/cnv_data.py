"""
This script writes copy number data across an entire genome to a series of .npz files.

This script requires a FASTA reference genome file, the windowed accessibility, and a BAM alignment file.
"""

import argparse
import numpy as np
import numpy.lib.recfunctions
from numpy import lib
from scipy import stats
import sklearn
import pysam
import pysamstats
import hmmlearn
from hmmlearn import hmm

AUTOSOMES = ('2R', '2L', '3R', '3L')
GENOME = ('2R', '2L', '3R', '3L', 'X')
WINDOW_SIZE = 300  # use non-overlapping 300 bp windows
    
def append_reference_gc(whole_coverage, fasta):
    """
    Appends a field to each recarray in whole_coverage indicating the %GC in each window in the reference genome.
    whole_coverage: A dictionary of coverage recarrays across the genome. Modified such that each coverage array has an extra
    field representing the reference %GC.
    fasta: A pysam.FastaFile object representing the reference genome
    returns: None
    modifies: whole_coverage
    """
    for chromosome in GENOME:
        fasta_gc = np.empty(len(whole_coverage[chromosome]), dtype='u1') # %GC content in the reference sequence in each window.
        for window_num, pos in enumerate(whole_coverage[chromosome].pos):
            start = pos - WINDOW_SIZE // 2 # start of window
            end = pos + WINDOW_SIZE // 2  # end of window
            countC = fasta.fetch(reference=chromosome, start=start, end=end).count('C')
            countG = fasta.fetch(reference=chromosome, start=start, end=end).count('G')
            gc_content = round((countC + countG) / WINDOW_SIZE * 100) # calculate gc content in the window
            fasta_gc[window_num] = gc_content
        # add a field to the coverage indicating the window %GC in the reference genome.
        whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "reference_gc", fasta_gc, dtypes='u1', asrecarray=True, usemask=False)

def append_accessibility(whole_coverage, accessibility_dir):
    """
    Appends a field to each recarray in whole_coverage indicating whether or not each window is accessible.
    
    whole_coverage: A dictionary of coverage recarrays across the genome. Modified such that each coverage array has an extra
    field representing whether or not each window is accessible.
    accessibility_dir: The directory where the windowed accessibility is stored.
    returns: None
    modifies: whole_coverage
    """
    for chromosome in GENOME: # append field indicating accessibility to coverage recarray.
        with open(f"{accessibility_dir}/accessibility.{chromosome}_{WINDOW_SIZE}.mask", 'rb') as window_mask:
            whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "accessible", np.load(window_mask), dtypes='?', asrecarray=True, usemask=False)

def mean_reads_gc_normalized(whole_coverage):
    """
    Determines the mean number of reads, normalized by gc content in the reference genome.
    
    whole_coverage: A dictionary of coverage recarrays across the genome.
    returns: A recarray containing the mean number of reads for each GC percentage.
    """
    mean_reads_by_gc = np.recarray(101, dtype=[("gc", 'u1'), ("mean_reads", '<f8'), ("num_windows", '<i4')])
    mean_reads_by_gc.gc = np.arange(0, 101, dtype='u1')
    mean_reads_by_gc.mean_reads = np.zeros(101, dtype='<f8')
    mean_reads_by_gc.num_windows = np.zeros(101, dtype='<i4')
    for chromosome in AUTOSOMES:
        for window in whole_coverage[chromosome]:
            if window.accessible: # only include accessible windows
                # add the read counts in the window to the element with the specified reference GC content
                mean_reads_by_gc[window.reference_gc].mean_reads += window.reads_all
                mean_reads_by_gc[window.reference_gc].num_windows += 1
    for i in mean_reads_by_gc:
        i.mean_reads = i.mean_reads/i.num_windows if i.num_windows != 0 else np.nan # compute mean
    return mean_reads_by_gc

def append_depth_normed(whole_coverage, mean_reads_by_gc):
    """
    Appends a field to each recarray in whole_coverage indicating the normalized coverage in each window.
    
    whole_coverage: A dictionary of coverage recarrays across the genome. Modified such that each coverage array has an extra
    field representing the normalized coverage.
    mean_reads_by_gc: A recarray containing the mean number of reads for each GC percentage (in the reference genome). Should 
    contain fields for reference_gc and mean_reads.
    returns: None
    modifies: whole_coverage
    """
    for chromosome in GENOME:
        # depth_normed is calculated as 2 times the reads in each window divided by the mean number of reads (normalized by gc content)
        depth_normed = [0 if np.isnan(mean_reads_by_gc[i.reference_gc].mean_reads) else 2*i.reads_all/mean_reads_by_gc[i.reference_gc].mean_reads for i in whole_coverage[chromosome]]
        whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "depth_normed", depth_normed, dtypes='<f8', asrecarray=True, usemask=False)

def append_filter(whole_coverage, whole_quality, mean_reads_by_gc):
    """
    Appends a field to each recarray in whole_coverage indicating whether or not each window is filtered.
    
    whole_coverage: A dictionary of coverage recarrays across the genome. Modified such that each coverage array has an extra
    field representing a filter.
    whole_quality: A dictionary of quality recarrays across the genome.
    mean_reads_by_gc: A recarray containing the mean number of reads for each GC percentage (in the reference genome). Should 
    contain fields for reference_gc and mean_reads.
    returns: None
    modifies: whole_coverage
    """
    for chromosome in GENOME:
        quality_filter = np.array(whole_quality[chromosome].reads_mapq0 >= whole_coverage[chromosome].reads_all*0.02)
        gc_filter = np.array(mean_reads_by_gc[whole_coverage[chromosome].reference_gc].num_windows < 100)
        coverage_filter = np.logical_or(quality_filter, gc_filter)
        whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "filtered", coverage_filter, dtypes='?', asrecarray=True, usemask=False)

def fit_hmm(depth_normed, transition_probability, variance, variance_fixed, max_copy_number=12, n_iter=0, params='st', init_params=''):
    """
    Predicts copy number states using a HMM.
    
    depth_normed: A list of normalized coverage values.
    transition_probability: The probability of a state transition to any state other than the current.
    variance: The variance of the normalized coverage, per copy.
    variance_fixed: The variance to be used for the zero copy number state.
    n_iter: The number of iterations to perform when fitting the model
    params: Parameters to the model that can be changed through fitting
    init_params: Parameters to the model that are initialized from the data.
    returns: A list containing the predicted copy number state.
    """
    # convenience variable
    min_copy_number = 0  # minimum copy number to consider in the model
    n_states = max_copy_number - min_copy_number + 1
    # construct the transition matrix
    transmat = np.zeros((n_states, n_states))
    transmat[:] = transition_probability # fill transition matrix with identical transition probabilities
    transmat[np.diag_indices(n_states)] = 1-((n_states-1)*transition_probability) # set diagonal elements equal to appropriate values
    # construct means and covariance
    means_list = range(n_states)
    means = np.array([[n] for n in means_list])
    covars = np.array([[variance*n + variance_fixed] for n in means_list])
    # setup HMM 
    model = hmm.GaussianHMM(n_states, 
                        covariance_type="diag", 
                        n_iter=n_iter, 
                        params=params,
                        init_params=init_params)
    model.means_ = means
    model.covars_ = covars
    model.transmat_ = transmat
    # fit HMM
    obs = np.column_stack([depth_normed])
    model.fit(obs)
    # predict hidden states
    h = model.predict(obs)
    return h

def coverage_variance(whole_coverage):
    """
    Determines the variance in the coverage across the genome.
    
    whole_coverage: A dictionary of coverage recarrays across the genome.
    returns: The variance in coverage across the genome, excluding windows marked an inaccessible.
    """
    return np.var([i.depth_normed for chromosome in AUTOSOMES for i in whole_coverage[chromosome] if i.accessible])

def save_cnv_data(whole_coverage, data_dir):
    """
    Saves predicted copy number data to numpy npz files.
    
    whole_coverage: A dictionary of coverage recarrays across the genome.
    data_dir: The directory where the copy number data should be saved to.
    returns: None
    """
    for chromosome in GENOME:
        filtered_pos = np.array([i.pos for i in whole_coverage[chromosome] if not i.filtered], dtype='<i4')
        filtered_depth = np.array([i.depth_normed for i in whole_coverage[chromosome] if not i.filtered], dtype='<f8')
        variance = coverage_variance(whole_coverage)
        copy_number = np.array(fit_hmm(filtered_depth, transition_probability=0.0001, variance=variance/2, variance_fixed=0.01), dtype='u1') # apply hmm
        data = {"filtered_pos" : filtered_pos, "filtered_depth" : filtered_depth, "copy_number" : copy_number}
        np.savez(f"{data_dir}/{chromosome}", **data)

def main():
    parser = argparse.ArgumentParser(description="calls CNVs in a sample. writes data to a file.")
    parser.add_argument("name", help="name of the sample.")
    parser.add_argument("fastafile", help="path to the reference genome")
    parser.add_argument("bamfile", help="path to the sample.")
    parser.add_argument("accessibility", help="directory where the windowed accessibility is saved.")
    parser.add_argument("cnv_data", help="directory where the cnv data should be saved to.")
    args = parser.parse_args()
    sample_name = args.name
    fasta_filename = args.fastafile
    bam_filename = args.bamfile
    accessibility_dir = args.accessibility
    data_dir = args.cnv_data
    fasta = pysam.FastaFile(fasta_filename)
    bam = pysam.AlignmentFile(bam_filename)
    # create dictionaries for coverage and quality across the genome.
    whole_coverage = {c : pysamstats.load_coverage_ext_binned(bam, fasta, chrom=c, window_size=WINDOW_SIZE) for c in GENOME}
    whole_quality = {c : pysamstats.load_mapq_binned(bam, fasta, chrom=c, window_size=WINDOW_SIZE) for c in GENOME}
    # process coverage data
    append_reference_gc(whole_coverage, fasta)
    append_accessibility(whole_coverage, accessibility_dir)
    mean_reads_by_gc = mean_reads_gc_normalized(whole_coverage)
    append_depth_normed(whole_coverage, mean_reads_by_gc)
    append_filter(whole_coverage, whole_quality, mean_reads_by_gc)
    # predict copy number states, save results
    save_cnv_data(whole_coverage, data_dir)

if __name__ == "__main__":
    main()
    