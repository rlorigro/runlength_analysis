from measure_runlength_distribution_from_runnie import runlength_encode_fasta, align_as_RLE
from modules.RunniePileupGenerator import PileupGenerator
from handlers.RunlengthHandler_v2 import RunlengthHandler
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from handlers.BamHandler import BamHandler
import random
import sys
import os


# Indexing RLE tuples
SEQUENCE = 0
LENGTHS = 1

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


def get_read_segments(fasta_handler, bam_handler, chromosome_name, pileup_start, pileup_end, read_data, runlength_ref_sequences):
    """
    Get raw read segments from a pair of coordinates given that each read has an aligned match at the start and end
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

    ref_lengths = runlength_ref_sequences[chromosome_name][LENGTHS][pileup_start:pileup_end]

    reads = bam_handler.get_reads(chromosome_name=chromosome_name,
                                  start=pileup_start,
                                  stop=pileup_end)

    segment_grabber = PileupGenerator(chromosome_name=chromosome_name,
                                      start_position=pileup_start,
                                      end_position=pileup_end,
                                      ref_sequence=ref_sequence,
                                      ref_lengths=ref_lengths,
                                      read_data=read_data,
                                      reads=reads)

    sequences, scales, shapes = segment_grabber.get_read_segments()

    return sequences, scales, shapes


def chunk_chromosome_coordinates(chromosome_length, chunk_size):
    chunks = list()

    index = 0
    while index < chromosome_length:
        start = index
        stop = min(chromosome_length, start + chunk_size)

        if stop - start > 2:
            chunks.append([start, stop])

        index = stop

    return chunks


def fetch_random_windows(n_windows,window_size, chromosome_length, fasta_handler, bam_handler, chromosome_name, runlength_ref_sequences, read_data):
    """
    Fetch sequences from n random windows of length window_size
    :param n_windows:
    :param chromosome_length:
    :param window_size:
    :param fasta_handler:
    :param bam_handler:
    :param chromosome_name:
    :param runlength_ref_sequences:
    :param read_data:
    :return:
    """

    i = 0
    while i < n_windows:
        sys.stdout.write("\r %.3f" % (100 * (float(i) / n_windows)))

        start = random.randint(0, (chromosome_length - window_size))
        stop = start + window_size - 1

        window = (start, stop)

        print("\nGENERATING WINDOW: ", window)

        sequences, scales, shapes = get_read_segments(fasta_handler=fasta_handler,
                                                      bam_handler=bam_handler,
                                                      chromosome_name=chromosome_name,
                                                      pileup_start=start,
                                                      pileup_end=stop,
                                                      runlength_ref_sequences=runlength_ref_sequences,
                                                      read_data=read_data)

        for k,key in enumerate(sequences):
            print(sequences[key][:10])
            print(scales[key][:10])
            print(shapes[key][:10])

        if len(sequences) > 0:
            i += 1


def main():
    # ref_fasta_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_ref.fasta"
    # runlength_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_runnie.out"

    ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    runlength_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/test/rad2_pass_runnie_4_5_6_7.out"

    output_parent_dir = "output/"
    output_dir = "runlength_matrix_from_runnie_output_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    assembly_fasta_filename_prefix = ".".join(os.path.basename(runlength_path).split(".")[:-1])
    runlength_read_fasta_filename = assembly_fasta_filename_prefix + "_rle.fasta"
    runlength_read_fasta_path = os.path.join(output_dir, runlength_read_fasta_filename)

    handler = RunlengthHandler(runlength_path)

    reads = handler.iterate_file(sequence_cutoff=sys.maxsize, print_status=True)
    read_data = dict()

    for r, read in enumerate(reads):
        read_data[read.id] = read

    print("\nRLE encoding reference sequence...")

    runlength_ref_sequences = runlength_encode_fasta(fasta_sequence_path=ref_fasta_path)

    assembly_vs_ref_bam_path = align_as_RLE(runlength_reference_path=runlength_ref_fasta_path,
                                            runlength_ref_sequences=runlength_ref_sequences,
                                            runlength_read_path=runlength_read_fasta_path,
                                            runlength_read_sequences=read_data,
                                            output_dir=output_dir)

    bam_handler = BamHandler(assembly_vs_ref_bam_path)
    fasta_handler = FastaHandler(runlength_ref_fasta_path)

    contig_names = fasta_handler.get_contig_names()
    chromosome_name = contig_names[0]
    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    sequences, scales, shapes = get_read_segments(fasta_handler=fasta_handler,
                                                  bam_handler=bam_handler,
                                                  chromosome_name=chromosome_name,
                                                  pileup_start=100000,
                                                  pileup_end=100000+100,
                                                  runlength_ref_sequences=runlength_ref_sequences,
                                                  read_data=read_data)

    for k, key in enumerate(sequences):
        print(key)
        print(sequences[key][:10])
        print(scales[key][:10])
        print(shapes[key][:10])

    # n_windows = 200
    # window_size = 50
    #
    # fetch_random_windows(n_windows=n_windows,
    #                      window_size=window_size,
    #                      chromosome_length=chromosome_length,
    #                      fasta_handler=fasta_handler,
    #                      bam_handler=bam_handler,
    #                      chromosome_name=chromosome_name,
    #                      runlength_ref_sequences=runlength_ref_sequences,
    #                      read_data=read_data)


if __name__ == "__main__":
    main()
