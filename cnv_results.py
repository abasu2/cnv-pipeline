import argparse
import numpy as np
from numpy import lib
from scipy import stats
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="generates results from CNV data.")
parser.add_argument("name", help="name of the sample.")
parser.add_argument("regions", help="filepath to the regions to examine for CNVs.")
parser.add_argument("cnv_data", help="directory where the cnv data is saved.")
parser.add_argument("-i", "--image", help="directory where the graphs should be saved to. images are named by cluster.")
args = parser.parse_args()
sample_name = args.name
data_dir = args.cnv_data

class Gene:
    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
    @staticmethod
    def make_gene(args): # make a Gene object from an argv type argument array.
        return Gene(args[1], args[2], args[3], args[4])
    def __repr__(self):
        return f"name: {self.name}, chrom: {self.chrom}, start: {self.start}, end: {self.end}"

class Cluster:
    def __init__(self, name, chrom, start, end, genes):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.genes = [gene for gene in genes if gene.chrom == self.chrom and gene.start >= self.start and gene.end <= self.end]
    @staticmethod
    def make_cluster(args, genes): # make a Cluster object from an argv type argument array and a gene array.
        return Cluster(args[1], args[2], args[3], args[4], genes)
    def __repr__(self):
        rep = f"name: {self.name}, chrom: {self.chrom}, start: {self.start}, end: {self.end}, genes:\n"
        for gene in genes:
            rep += gene.__repr__() + "\n"
        return rep
    
with open(args.regions, 'r') as file: # read regions genes and clusters from regions file.
    next(file)
    lines = [s.split(", ") for s in file.read().splitlines()]
genes = [Gene.make_gene(line) for line in lines if line[0] == "gene"] # save genes in the file
clusters = [Cluster.make_cluster(line, genes) for line in lines if line[0] == "cluster"] # save gene clusters in the file

autosomes = ['2R', '2L', '3R', '3L']
genome = ['2R', '2L', '3R', '3L', 'X']
target_chroms = np.unique([c.chrom for c in genes]) # only the chromosomes which contain a gene of interest

for chromosome in target_chroms:
    npzfile = np.load(f"{data_dir}/{chromosome}.npz")
    filtered_pos = npzfile["filtered_pos"]
    filtered_depth = npzfile["filtered_depth"]
    copy_number = npzfile["copy_number"]
    on_chrom = lambda obj: obj.chrom == chromosome # checks if a Gene/Cluster is on chromosome
    for gene in filter(on_chrom, genes): # selects only genes on chromosome
        indices = [i for i, pos in enumerate(filtered_pos) if gene.start <= pos <= gene.end] # indices in filtered_pos corresponding to the gene
        duplication = stats.mode(np.take(copy_number, indices), axis=None).mode - 2 # copy number state of the gene. calculated as the mode minus 2.
        duplication = duplication.item()
        print(f"{repr(gene)}, duplication: {duplication}")
    if args.image is not None:
        for cluster in filter(on_chrom, clusters): # selects only clusters on chromosome
            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot(111)
            #           black      brown    darkgreen     red      darkblue      fuschia  springgreen    aqua      orange    palegold   hotpink   steelblue
            colors = ("#000000", "#8b4513", "#006400", "#ff0000", "#00008b", "#ff00ff", "#00ff7f", "#00ffff", "#ff8c00", "#eee8aa", "#ff69b4", "#4682b4")
            for i, color in zip(range(13), colors):
                ax.plot(filtered_pos[copy_number == i], 
                        filtered_depth[copy_number == i], 
                        color=color, marker='o', linestyle=' ', alpha=.2)
            ax.plot(filtered_pos, copy_number, linestyle='-', linewidth=2, color='k')
            ax.set_title(f"predicted copy number in {sample_name} around {cluster.name}", fontsize=20)
            ax.set_xlabel(f"{chromosome} position (bp)", fontsize=16)
            ax.set_ylabel("copy number", fontsize=16)
            ax.tick_params(labelsize=15)
            ax.set_yticks(range(13))
            ax.grid(axis='y')
            cluster_len = cluster.end - cluster.start
            normal_cnv_pos = filtered_pos[copy_number == 2] # marks the positions where there is no CNV
            # figure out how to make the lower xlim and upper xlim bound the CNV area + some no CNV area.
            lower_xlim = min(cluster.start - int(cluster_len*0.15), cluster.start)
            ax.set_xlim(left=cluster.start - int(cluster_len*0.15), right=cluster.end + int(cluster_len*0.15)) # sets limits to be slightly larger than cluster length
            gene_fontsize = min(min([700*(g.end - g.start) / cluster_len for g in cluster.genes]), 16) # sets the fontsize for all gene labels in the cluster
            for i, gene in enumerate(cluster.genes):
                plt.axvspan(gene.start, gene.end, color="#808080", alpha=0.5) # highlights individual genes.
                gene_midpoint = (gene.end + gene.start) // 2
                # label individual genes
                plt.text(gene_midpoint, 10, gene.name, 
                         verticalalignment = "center", horizontalalignment="center", rotation="vertical", fontsize=gene_fontsize)
                plt.savefig(f"{args.image}/{cluster.name}.png")