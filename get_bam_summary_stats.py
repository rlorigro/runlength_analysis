from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from collections import defaultdict
import numpy
import math


class IterativeHistogram:
    def __init__(self, start, stop, n_bins, unbounded_upper_bin=False, unbounded_lower_bin=False, include_upper_edge=True):
        self.start = start
        self.stop = stop
        self.n_bins = n_bins
        self.histogram = numpy.zeros(n_bins)
        self.bin_size = (stop - start)/n_bins
        self.edges = [start + self.bin_size*i for i in range(n_bins+1)]

        self.unbounded_upper_bin = unbounded_upper_bin
        self.unbounded_lower_bin = unbounded_lower_bin
        self.include_upper_edge = include_upper_edge

    def get_bin(self, x):
        # find index of bin by normalizing and centering the value w.r.t. bin edges and add value to that bin
        bin_index = int(math.floor((x - self.start)/self.bin_size))

        if x == self.stop and self.include_upper_edge:
            bin_index = self.n_bins - 1

        if self.unbounded_lower_bin and x < self.start:
            bin_index = 0

        if self.unbounded_upper_bin and x > self.stop:
            bin_index = self.n_bins - 1

        return bin_index

    def update(self, x):
        bin_index = self.get_bin(x)

        if 0 <= bin_index <= (self.n_bins-1):
            self.histogram[bin_index] += 1
        else:
            print("WARNING: value not within user-defined bins:", x, "\tBin index:", bin_index)

    def get_histogram(self):
        return self.histogram

    def get_normalized_histogram(self):
        total = sum(self.histogram)
        normalized_histogram = self.histogram/numpy.sum(self.histogram)

        return normalized_histogram


def get_cigar_name_from_code(cigar_code):
    if cigar_code == 0:
        name = "MATCH"
    elif cigar_code == 1:
        name = "INS"
    elif cigar_code == 2:
        name = "DEL"
    elif cigar_code == 3:
        name = "REF_SKIP"
    elif cigar_code == 4:
        name = "SOFT_CLIP"
    elif cigar_code == 5:
        name = "HARD_CLIP"
    elif cigar_code == 6:
        name = "PAD"
    elif cigar_code == 7:
        name = "EQUAL"
    elif cigar_code == 8:
        name = "DIFF"
    elif cigar_code == 9:
        name = "BACK"
    elif cigar_code == 10:
        name = "NM"
    else:
        exit("ERROR: invalid cigar_code: %d" % cigar_code)

    return name


def count_cigar_operations(reads, chromosome_name, chromosomal_cigar_counts, n_alignments, n_primary, n_supplementary, n_secondary, map_qualities):
    """
    Given a set of pysam read objects, generate data from each read
    :param reads:
    :param chromosome_name:
    :param chromosomal_cigar_counts:
    :return:
    """

    for read in reads:
        n_alignments += 1
        map_qualities.update(read.mapping_quality)

        if read.is_secondary:
            n_secondary += 1

        if read.is_read1:
            n_primary += 1

        if read.is_supplementary:
            n_supplementary += 1

        if read.mapping_quality > 0:
            cigar_tuples = read.cigartuples

            for c, cigar in enumerate(cigar_tuples):
                cigar_code = cigar[0]
                length = cigar[1]

                chromosomal_cigar_counts[chromosome_name][cigar_code] += 1

    return chromosomal_cigar_counts, n_alignments, n_primary, n_supplementary, n_secondary, map_qualities


def parse_bam(bam_path, reference_path):
    """
    Iterate a BAM file and count summary stats from that file
    :param bam_path:
    :param reference_path:
    :return:
    """
    fasta_handler = FastaHandler(reference_path)
    chromosome_names = fasta_handler.get_contig_names()

    chromosomal_cigar_counts = defaultdict(lambda: defaultdict(int))

    n_alignments = 0
    n_primary = 0
    n_supplementary = 0
    n_secondary = 0

    map_qualities = IterativeHistogram(start=0, stop=60, n_bins=6)

    for chromosome_name in chromosome_names:
        bam_handler = BamHandler(bam_path)

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=0, stop=chromosome_length)

        chromosomal_cigar_counts, \
        n_alignments, \
        n_primary, \
        n_supplementary, \
        n_secondary, \
        map_qualities = count_cigar_operations(reads=reads,
                                               chromosome_name=chromosome_name,
                                               chromosomal_cigar_counts=chromosomal_cigar_counts,
                                               n_alignments=n_alignments,
                                               n_primary=n_primary,
                                               n_supplementary=n_supplementary,
                                               n_secondary=n_secondary,
                                               map_qualities=map_qualities)

    return chromosomal_cigar_counts, n_alignments, n_primary, n_supplementary, n_secondary, map_qualities


def main():
    bam_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_runnie_output_2019_4_5_11_16_8_642422/runnie_subset_train_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"
    reference_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"

    chromosomal_cigar_counts, \
    n_alignments, \
    n_primary, \
    n_supplementary, \
    n_secondary, \
    map_qualities = parse_bam(bam_path=bam_path, reference_path=reference_path)

    for key in chromosomal_cigar_counts:
        print(key)
        for cigar_code in chromosomal_cigar_counts[key]:
            print(get_cigar_name_from_code(cigar_code), "\t", chromosomal_cigar_counts[key][cigar_code])

    print()
    print("n_alignments:\t\t", n_alignments)
    print("n_primary:\t\t", n_primary)
    print("n_supplementary:\t", n_supplementary)
    print("n_secondary:\t\t", n_secondary)
    print("map_qualities:\t\t", map_qualities.get_normalized_histogram())


if __name__ == "__main__":
    main()
