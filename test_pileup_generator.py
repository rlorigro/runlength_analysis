from modules.PileupGenerator import INSERT_CHAR, DELETE_CHAR
from modules.PileupGenerator import PileupGenerator
from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from matplotlib import pyplot
import numpy


def get_encoding(character):
    if character == "A":
        encoding = 0
    elif character == "C":
        encoding = 1
    elif character == "G":
        encoding = 2
    elif character == "T":
        encoding = 3
    elif character == INSERT_CHAR:
        encoding = 5
    elif character == DELETE_CHAR:
        encoding = 5
    else:
        raise("ERROR: unrecognized character: %s" % character)

    return encoding


def get_aligned_segments(fasta_handler, bam_handler, chromosome_name, pileup_start, pileup_end, include_ref=False):
    """
    Get read segments from a pair of coordinates given that each read has an aligned match at the start and end
    coordinate
    :param fasta_handler:
    :param bam_handler:
    :param chromosome_name:
    :param pileup_start:
    :param pileup_end:
    :return:
    """
    ref_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                              start=pileup_start,
                                              stop=pileup_end + 1)

    reads = bam_handler.get_reads(chromosome_name=chromosome_name,
                                  start=pileup_start,
                                  stop=pileup_end)

    segment_grabber = PileupGenerator(chromosome_name=chromosome_name,
                                      start_position=pileup_start,
                                      end_position=pileup_end,
                                      ref_sequence=ref_sequence,
                                      reads=reads)

    # if a reference sequence is intended to be added to the pileup, then leave a space for it
    if include_ref:
        segment_grabber.max_coverage -= 1

    sequence_dictionary = segment_grabber.get_aligned_read_segments()

    return sequence_dictionary


def main():
    # bam_file_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_3_27_14_59_24_409353/sequence_subset_test_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"
    # ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"

    bam_file_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_runnie_output_2019_4_8_17_33_14_191911/runnie_subset_test_60x_10kb_rle_VS_refEcoli_rle.sorted.bam"
    ref_fasta_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_runnie_output_2019_4_8_17_33_14_191911/refEcoli_rle.fasta"
    # -------------------------------------------------------------------------

    fasta_handler = FastaHandler(ref_fasta_path)
    contig_names = fasta_handler.get_contig_names()
    chromosome_name = contig_names[0]

    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    bam_handler = BamHandler(bam_file_path)
    fasta_handler = FastaHandler(ref_fasta_path)

    pileup_start = 0
    pileup_end = pileup_start+1000      # add random variation here ?

    aligned_segments = get_aligned_segments(fasta_handler=fasta_handler,
                                            bam_handler=bam_handler,
                                            chromosome_name=chromosome_name,
                                            pileup_start=pileup_start,
                                            pileup_end=pileup_end,
                                            include_ref=True)

    encoding = list()
    for alignment in aligned_segments.values():
        encoding.append(list(map(get_encoding, alignment)))

    encoding = -numpy.array(encoding, dtype=numpy.float)

    pyplot.imshow(encoding)
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
