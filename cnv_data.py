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

parser = argparse.ArgumentParser(description="calls CNVs in a sample. writes data to a file.")
parser.add_argument("name", help="name of the sample.")
parser.add_argument("bamfile", help="filename of the sample.")
parser.add_argument("accessibility", help="directory where the windowed accessibility is saved.")
parser.add_argument("cnv_data", help="directory where the cnv data should be saved to.")
args = parser.parse_args()
sample_name = args.name
bam_filename = args.bamfile
accessibility_dir = args.accessibility
data_dir = args.cnv_data

fasta_filename = "/share/lanzarolab/users/abasu/reference/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"  # path to the reference genome FASTA file
AUTOSOMES = ('2R', '2L', '3R', '3L')
GENOME = ('2R', '2L', '3R', '3L', 'X')
WINDOW_SIZE = 300  # use non-overlapping 300 bp windows
fasta = pysam.FastaFile(fasta_filename)
bam = pysam.AlignmentFile(bam_filename)
whole_coverage = {c : pysamstats.load_coverage_ext_binned(bam, fasta, chrom=c, window_size=WINDOW_SIZE) for c in genome}
whole_quality = {c : pysamstats.load_mapq_binned(bam, fasta, chrom=c, window_size=WINDOW_SIZE) for c in genome}

for chromosome in GENOME:
    fasta_gc = np.empty(len(whole_coverage[chromosome]), dtype='u1') # %GC content in the reference sequence in each window.
    for window_num, pos in enumerate(whole_coverage[chromosome].pos):
        start = pos - WINDOW_SIZE//2
        end = pos + WINDOW_SIZE//2 
        countC = fasta.fetch(reference=chromosome, start=start, end=end).count('C')
        countG = fasta.fetch(reference=chromosome, start=start, end=end).count('G')
        gc_content = round((countC + countG) / WINDOW_SIZE * 100)
        fasta_gc[window_num] = gc_content
    # add a field to the coverage indicating the window %GC in the reference genome.
    whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "reference_gc", fasta_gc, dtypes='u1', asrecarray=True, usemask=False)

for chromosome in GENOME: # append field indicating accessibility to coverage recarray.
    with open(f"{accessibility_dir}/accessibility.{chromosome}_{window_size}.mask", 'rb') as window_mask:
        whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "accessible", np.load(window_mask), dtypes='?', asrecarray=True, usemask=False)
        
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

for chromosome in GENOME:
    depth_normed = [0 if np.isnan(mean_reads_by_gc[i.reference_gc].mean_reads) else 2*i.reads_all/mean_reads_by_gc[i.reference_gc].mean_reads 
                    for i in whole_coverage[chromosome]]
    whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "depth_normed", depth_normed, dtypes='<f8', asrecarray=True, usemask=False)

for chromosome in GENOME:
    quality_filter = np.array(whole_quality[chromosome].reads_mapq0 >= whole_coverage[chromosome].reads_all*0.02)
    gc_filter = np.array(mean_reads_by_gc[whole_coverage[chromosome].reference_gc].num_windows < 100)
    coverage_filter = np.logical_or(quality_filter, gc_filter)
    whole_coverage[chromosome] = np.lib.recfunctions.append_fields(whole_coverage[chromosome], "filtered", coverage_filter, dtypes='?', asrecarray=True, usemask=False)

def fit_hmm(depth_normed,  # normalised coverage array 
            transition_probability,  # probability of state transition
            variance,  # variance per copy 
            variance_fixed,  # variance for the zero copy number state 
            max_copy_number=12,  # maximum copy number to consider in the model 
            n_iter=0,  # number of iterations to perform when fitting the model
            params='st',  # parameters that can be changed through fitting 
            init_params=''  # parameters that are initialised from the data
           ):
    
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

coverage_variance = np.var([i.depth_normed for chromosome in AUTOSOMES for i in whole_coverage[chromosome] if i.accessible])
                
for chromosome in GENOME:
    filtered_pos = np.array([i.pos for i in whole_coverage[chromosome] if not i.filtered], dtype='<i4')
    filtered_depth = np.array([i.depth_normed for i in whole_coverage[chromosome] if not i.filtered], dtype='<f8')
    copy_number = np.array(fit_hmm(filtered_depth, transition_probability=0.0001, variance=coverage_variance/2, variance_fixed=0.01), dtype='u1') # apply hmm
    data = {"filtered_pos" : filtered_pos, "filtered_depth" : filtered_depth, "copy_number" : copy_number}
    np.savez(f"{data_dir}/{chromosome}", **data)