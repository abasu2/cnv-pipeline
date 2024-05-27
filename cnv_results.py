import argparse
import numpy as np
from numpy import lib
from scipy import stats
import csv
import matplotlib.pyplot as plt

class Gene:
    def __init__(self, name, chrom, start, end):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
    @classmethod
    def make_gene(cls, args): # make a Gene object from an argv type argument array.
        return cls(args[1], args[2], args[3], args[4])
    def __repr__(self):
        return f"name: {self.name}, chrom: {self.chrom}, start: {self.start}, end: {self.end}"
class Cluster:
    def __init__(self, name, chrom, start, end, genes):
        self.name = name
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.genes = [gene for gene in genes if gene.chrom == self.chrom and gene.start >= self.start and gene.end <= self.end]
    @classmethod
    def make_cluster(cls, args, genes): # make a Cluster object from an argv type argument array and a gene array.
        return cls(args[1], args[2], args[3], args[4], genes)
    def __repr__(self):
        rep = f"name: {self.name}, chrom: {self.chrom}, start: {self.start}, end: {self.end}, genes:\n"
        for gene in genes:
            rep += gene.__repr__() + "\n"
        return rep

NPZFILE = np.load(f"{data_dir}/{chromosome}.npz")
FILTERED_POS = npzfile["filtered_pos"]
FILTERED_DEPTH = npzfile["filtered_depth"]
COPY_NUMBER = npzfile["copy_number"]
COLORS = ("#000000", # black 
          "#8b4513", # brown
          "#006400", # darkgreen
          "#ff0000", # red
          "#00008b", # darkblue
          "#ff00ff", # fuschia
          "#00ff7f", # springgreen
          "#00ffff", # aqua
          "#ff8c00", # orange
          "#eee8aa", # palegold
          "#ff69b4", # hotpink
          "#4682b4") # steelblue
    
def read_regions(region_file):
    """
    Reads and parses a region file.
    region_file: A text file containing the regions to be examined.
    returns: A list of the lines in the regions file, containing type, name, chromosome, start, and end fields.
    """
    with open(region_file, 'r') as file: # read regions genes and clusters from regions file.
        next(file)
        lines = [s.split(", ") for s in file.read().splitlines()]
    return lines

def create_images(image_dir, genes, clusters, sample_name):
    """
    Creates images of coverage and copy number for a sample.
    image_dir: The directory where the images should be saved.
    genes: A list of gene objects of interest.
    clusters: A list of gene clusters. One image is made per cluster.
    sample_name: The name of the sample being analyzed.
    """
    TARGET_CHROMS = np.unique([c.chrom for c in genes]) # only the chromosomes which contain a gene of interest
    for chromosome in TARGET_CHROMS:
        on_chrom = lambda obj: obj.chrom == chromosome # predicate that checks if a Gene/Cluster is on chromosome
        for cluster in filter(on_chrom, clusters): # selects only clusters on chromosome
            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot(111)
            for i, COLORS in zip(range(13), COLORS):
                ax.plot(FILTERED_POS[COPY_NUMBER == i], 
                        FILTERED_DEPTH[COPY_NUMBER == i], 
                        color=color, marker='o', linestyle=' ', alpha=.2)
            ax.plot(FILTERED_POS, COPY_NUMBER, linestyle='-', linewidth=2, color='k')
            ax.set_title(f"predicted copy number in {sample_name} around {cluster.name}", fontsize=20)
            ax.set_xlabel(f"{chromosome} position (bp)", fontsize=16)
            ax.set_ylabel("copy number", fontsize=16)
            ax.tick_params(labelsize=15)
            ax.set_yticks(range(13))
            ax.grid(axis='y')
            cluster_len = cluster.end - cluster.start
            lower_xlim = min(cluster.start - int(cluster_len*0.15), cluster.start)
            # sets limits to be slightly larger than cluster length
            ax.set_xlim(left=cluster.start - int(cluster_len*0.15), right=cluster.end + int(cluster_len*0.15))
            # sets the fontsize for all gene labels in the cluster
            gene_fontsize = min(min([700*(g.end - g.start) / cluster_len for g in cluster.genes]), 16) 
            for i, gene in enumerate(cluster.genes):
                plt.axvspan(gene.start, gene.end, color="#808080", alpha=0.5) # highlights individual genes.
                gene_midpoint = (gene.end + gene.start) // 2
                # label individual genes
                plt.text(gene_midpoint, 10, gene.name, 
                         verticalalignment = "center", horizontalalignment="center", rotation="vertical", fontsize=gene_fontsize)
                plt.savefig(f"{image_dir}/{cluster.name}.png")

def serialize_results(csv_file, genes, sample_name):
    """
    Serializes the copy number results in a CSV file.
    csv_file: The name of the csv file where the results should be serialized.
    genes: A list of gene objects of interest.
    sample_name: The name of the sample.
    """
    # create the csv file with its fields if it does not exist.
    HEADER = ["Sample Name"] + [g.name for g in genes]
    try:
        with open(csv_file, 'x') as file:
            writer = csv.DictWriter(file, fieldnames=HEADER)
            writer.writeheader(HEADER)
    except FileExistsError:
        pass
    # dictionary containing gene names and duplications
    gene_duplications = {g.name : 0 for g in genes}
    TARGET_CHROMS = np.unique([c.chrom for c in genes]) # only the chromosomes which contain a gene of interest
    for chromosome in TARGET_CHROMS:
        on_chrom = lambda obj: obj.chrom == chromosome # checks if a Gene/Cluster is on chromosome
        for gene in filter(on_chrom, genes): # selects only genes on chromosome
            # indices in filtered_pos corresponding to the gene
            indices = [i for i, pos in enumerate(FILTERED_POS) if gene.start <= pos <= gene.end]
            # copy number state of the gene. calculated as the mode minus 2.
            duplication = stats.mode(np.take(copy_number, indices), axis=None).mode - 2
            duplication = duplication.item()
            gene_duplications[gene.name] = duplication
            print(f"{repr(gene)}, duplication: {duplication}")
    with open(csv_file, 'w+') as file:
        writer = csv.DictWriter(file, fieldnames=HEADER, restval=sample_name)
        writer.writerows(gene_duplications)
    
def main():
    parser = argparse.ArgumentParser(description="generates results from CNV data.")
    parser.add_argument("name", help="name of the sample.")
    parser.add_argument("regions", help="filepath to the regions to examine for CNVs.")
    parser.add_argument("cnv_data", help="directory where the CNV data is saved.")
    parser.add_argument("result_csv", help="path to the csv file where the results should be saved to.")
    parser.add_argument("-i", "--image", help="directory where the graphs should be saved to. images are named by cluster.")
    args = parser.parse_args()
    sample_name = args.name
    region_file = args.regions
    csv_file = args.result_csv
    data_dir = args.cnv_data
    region_lines = read_regions(region_file)
    genes = [Gene.make_gene(line) for line in region_lines if line[0] == "gene"] # save genes in the file
    clusters = [Cluster.make_cluster(line, genes) for line in region_lines if line[0] == "cluster"] # save gene clusters in the file
    if args.image is not None:
        create_images(args.image, genes, clusters, sample_name)
    serialize_results(csv_file, genes, sample_name)

if __name__ == "__main__":
    main()