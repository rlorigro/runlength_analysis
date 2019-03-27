from matplotlib import pyplot
import pysam
import math
import sys


def get_coverage(bam_file_path, chromosome_name, start, stop):
    sam_file = pysam.AlignmentFile(bam_file_path, "rb")

    pileup_columns = sam_file.pileup(chromosome_name, start, stop)
    coverages = list()

    for pileup_column in pileup_columns:
        position = pileup_column.pos
        coverage = pileup_column.nsegments

        if start < position < stop:
            coverages.append(coverage)

    sam_file.close()

    return coverages


def main():
    # BAM of reads aligned to some sequence
    bam_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_3_27_13_13_54_735626/sequence_subset_train_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"

    # Fasta containing the sequences that are the reference in the BAM
    fasta_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_3_27_11_15_11_736904/refEcoli_rle.fasta"

    fasta_handler = pysam.FastaFile(fasta_path)
    chromosome_names = fasta_handler.references
    print(chromosome_names)

    n_samples = 100

    coverages = list()
    for name in chromosome_names:
        length = fasta_handler.get_reference_length(name)

        start = 0
        stop = length
        step_size = int(round(length/n_samples))

        steps = range(start, stop, step_size)

        for c,coord in enumerate(steps):
            coord = int(math.floor(coord))
            coverage = get_coverage(bam_file_path=bam_path, chromosome_name=name, start=coord, stop=coord+2)[0]
            coverages.append(coverage)
            sys.stderr.write("\r %.2f%%" % (c+1/n_samples*100))

        sys.stderr.write("\n")

    pyplot.plot(coverages)
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
