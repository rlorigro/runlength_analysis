from measure_runlength_distribution_from_fasta import runlength_encode_fasta, align_as_RLE
from modules.RunlengthClassifier import RunlengthClassifier
from modules.RunlengthPileupGenerator import INSERT_CHAR, DELETE_CHAR
from modules.RunlengthPileupGenerator import PileupGenerator
from handlers.RunlengthHandler_v2 import RunlengthHandler
from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from handlers.BamHandler import BamHandler
from matplotlib import pyplot
import numpy
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


def get_accuracy_from_confusion_matrix(matrix):
    diagonal = numpy.eye(matrix.shape[0], matrix.shape[1])

    diagonal_mask = (diagonal == 1)
    off_diagonal_mask = (diagonal == 0)

    diagonal_confusion = matrix[diagonal_mask]
    off_diagonal_confusion = matrix[off_diagonal_mask]

    diagonal_sum = numpy.sum(diagonal_confusion)
    off_diagonal_sum = numpy.sum(off_diagonal_confusion)

    accuracy = diagonal_sum/(diagonal_sum + off_diagonal_sum)

    return accuracy


def get_runlength_confusion(true_lengths, predicted_lengths, max_length=50):
    confusion = numpy.zeros([max_length+1, max_length+1])

    for true_length,predicted_length in zip(true_lengths, predicted_lengths):
        true_length = min(max_length, true_length)
        predicted_length = min(max_length, predicted_length)

        confusion[true_length,predicted_length] += 1

    return confusion


def plot_runlength_pileup(sequences, lengths, ref_sequence, ref_lengths):
    """
    Given 2d numpy float matrices, make an organized plot of pileup data
    :param sequences:
    :param scales:
    :param shapes:
    :param lengths:
    :param ref_sequence:
    :param ref_lengths:
    :return:
    """
    figure, axes = pyplot.subplots(nrows=6, sharex=True)

    axes[0].imshow(ref_sequence)
    axes[1].imshow(sequences)
    axes[2].imshow(ref_lengths)
    axes[3].imshow(lengths)

    axes[0].set_ylabel("Ref Nucleotide")
    axes[1].set_ylabel("Read Nucleotide")
    axes[2].set_ylabel("Ref length")
    axes[3].set_ylabel("Read Length")

    pyplot.show()
    pyplot.close()


def mode(x):
    """
    Given an array of repeat observations e.g.: [1, 2, 2, 1, 2], find the most frequent observation
    :param x:
    :return:
    """
    unique, inverse = numpy.unique(x, return_inverse=True)
    bincount = numpy.bincount(inverse)
    mode_index = bincount.argmax()
    mode = unique[mode_index]

    # print()
    # print(x)
    # print(unique)
    # print(inverse)
    # print(bincount)
    # print(mode)

    return mode


def get_consensus_from_runlength_pileup_encoding(length_classifier, sequence_encoding, length_encoding, reversal_encoding, bayesian=True):
    # print(sequence_encoding)
    # print(length_encoding)

    window_length = sequence_encoding.shape[1]

    consensus_sequence = list()
    consensus_lengths = list()

    for i in range(window_length):
        sequence_observations = sequence_encoding[:,i]
        length_observations = length_encoding[:,i]

        base_consensus = mode(sequence_observations)

        # print("sequence_observations\t", sequence_observations)
        # print("length_observations\t", length_observations)
        # print("reversal_encoding\t", reversal_encoding)
        # print("base consensus\t\t", base_consensus)

        if base_consensus == BASE_TO_INDEX["-"]:
            length_consensus = 0

        else:
            if bayesian:
                normalized_y_log_likelihoods, length_consensus = length_classifier.predict(x=length_observations,
                                                                                           character_index=base_consensus,
                                                                                           reversal=reversal_encoding)
            else:
                length_consensus = int(mode(length_observations))

        consensus_sequence.append(base_consensus)
        consensus_lengths.append(length_consensus)

        # print("length normalized_y_log_likelihoods\n", 10**normalized_y_log_likelihoods[:10])
        # print("length consensus\t", length_consensus)
        # print()

    return consensus_sequence, consensus_lengths


def map_parameters_to_mode(scale_shape_tuple, max_runlength=50):
    x = numpy.arange(0,max_runlength)

    mode = RunlengthHandler.calculate_numerical_mode(scale=scale_shape_tuple[0], shape=scale_shape_tuple[1], x_range=x)

    return mode


def convert_parameters_to_modes(scales, shapes):
    modes = list()
    for read_id in scales.keys():
        modes.append(list(map(map_parameters_to_mode, (scales[read_id], shapes[read_id]))))

    return modes


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
        encoding = 4
    elif character == DELETE_CHAR:
        encoding = 4
    else:
        raise("ERROR: unrecognized character: %s" % character)

    return encoding


def get_aligned_segments(fasta_handler, bam_handler, chromosome_name, pileup_start, pileup_end, read_data, runlength_ref_sequences):
    """
    Get aligned read segments from a pair of coordinates given that each read has an aligned match at the start and end
    coordinate.
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

    sequence_data = segment_grabber.get_aligned_read_segments()

    return sequence_data


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


def main():
    # ref_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_runlength_test_2019_3_25_13_8_0_341509_ref.fasta"
    # read_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_runlength_test_2019_3_25_13_8_0_341509_reads.fasta"
    # matrix_path = "/home/ryan/code/runnie_parser/output/runlength_matrix_from_assembly_contigs_2019_3_19_13_29_14_657613/probability_matrices_2019_3_19_13_29_19_362916.csv"

    ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    read_fasta_path = "/home/ryan/code/runlength_analysis/data/sequence_subset_ecoli_guppy-runnie_60x_test.fastq"
    matrix_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_sequence_2019_4_5_15_29_28_403950/probability_matrices_2019_4_5_15_35_57_920301.csv"

    output_parent_dir = "output/"
    output_dir = "runlength_matrix_from_sequence_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    read_fasta_filename_prefix = ".".join(os.path.basename(read_fasta_path).split(".")[:-1])
    runlength_read_fasta_filename = read_fasta_filename_prefix + "_rle.fasta"
    runlength_read_fasta_path = os.path.join(output_dir, runlength_read_fasta_filename)

    runlength_ref_sequences = runlength_encode_fasta(fasta_sequence_path=ref_fasta_path)
    runlength_read_sequences = runlength_encode_fasta(fasta_sequence_path=read_fasta_path)

    read_vs_ref_bam_path = align_as_RLE(runlength_reference_path=runlength_ref_fasta_path,
                                        runlength_ref_sequences=runlength_ref_sequences,
                                        runlength_read_path=runlength_read_fasta_path,
                                        runlength_read_sequences=runlength_read_sequences,
                                        output_dir=output_dir)

    bam_handler = BamHandler(read_vs_ref_bam_path)
    fasta_handler = FastaHandler(runlength_ref_fasta_path)

    contig_names = fasta_handler.get_contig_names()
    chromosome_name = contig_names[0]
    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    windows = chunk_chromosome_coordinates(chromosome_length=chromosome_length, chunk_size=1000)

    # Initialize empty confusion matrices
    total_confusion = get_runlength_confusion([],[],10)
    total_modal_confusion = get_runlength_confusion([],[],10)

    length_classifier = RunlengthClassifier(matrix_path)

    print("reading BAM")
    for pileup_start, pileup_end in windows[:10]:
        print("window", pileup_start, pileup_end)

        sys.stderr.write("\r%s" % pileup_start)
        aligned_ref_sequence, aligned_ref_lengths, aligned_sequences, aligned_lengths, reversal_statuses = \
            get_aligned_segments(fasta_handler=fasta_handler,
                                 bam_handler=bam_handler,
                                 chromosome_name=chromosome_name,
                                 pileup_start=pileup_start,
                                 pileup_end=pileup_end,
                                 runlength_ref_sequences=runlength_ref_sequences,
                                 read_data=runlength_read_sequences)

        sequence_encoding = list()
        length_encoding = list()
        reversal_encoding = list()

        # No reads here?
        if len(aligned_sequences) == 0:
            continue

        # print("REF\t", "".join(aligned_ref_sequence))
        for read_id in aligned_sequences.keys():
            # print("READ\t","".join(aligned_sequences[read_id]))
            sequence_encoding.append(list(map(get_encoding, aligned_sequences[read_id])))
            length_encoding.append(aligned_lengths[read_id])
            reversal_encoding.append(reversal_statuses[read_id])

        ref_sequence_encoding = [list(map(get_encoding, aligned_ref_sequence))]
        ref_lengths_encoding = [aligned_ref_lengths]

        ref_sequence_encoding = numpy.array(ref_sequence_encoding, dtype=numpy.int)
        ref_length_encoding = numpy.array(ref_lengths_encoding, dtype=numpy.int)
        sequence_encoding = numpy.array(sequence_encoding, dtype=numpy.int)
        length_encoding = numpy.array(length_encoding, dtype=numpy.float)
        reversal_encoding = numpy.array(reversal_encoding, dtype=numpy.bool)

        ref_sequence_encoding = numpy.atleast_2d(ref_sequence_encoding)
        ref_length_encoding = numpy.atleast_2d(ref_length_encoding)
        sequence_encoding = numpy.atleast_2d(sequence_encoding)
        length_encoding = numpy.atleast_2d(length_encoding)

        # plot_runlength_pileup(sequences=-sequence_encoding,
        #                       lengths=length_encoding,
        #                       ref_sequence=-ref_sequence_encoding,
        #                       ref_lengths=ref_length_encoding)

        consensus_sequence, consensus_lengths = \
            get_consensus_from_runlength_pileup_encoding(length_classifier=length_classifier,
                                                         sequence_encoding=sequence_encoding,
                                                         length_encoding=length_encoding,
                                                         reversal_encoding=reversal_encoding)

        modal_consensus_sequence, modal_consensus_lengths = \
            get_consensus_from_runlength_pileup_encoding(length_classifier=length_classifier,
                                                         sequence_encoding=sequence_encoding,
                                                         length_encoding=length_encoding,
                                                         reversal_encoding=reversal_encoding,
                                                         bayesian=False)

        print()
        print("PREDICTED\t",consensus_lengths[:10])
        print("TRUE\t\t",aligned_ref_lengths[:10])

        confusion = get_runlength_confusion(true_lengths=aligned_ref_lengths,
                                            predicted_lengths=consensus_lengths,
                                            max_length=10)

        total_confusion += confusion

        modal_confusion = get_runlength_confusion(true_lengths=aligned_ref_lengths,
                                                  predicted_lengths=modal_consensus_lengths,
                                                  max_length=10)

        total_modal_confusion += modal_confusion

        # except Exception as e:
        #     print(e)
        #     continue
    print()

    accuracy = get_accuracy_from_confusion_matrix(total_confusion)

    print("Bayes:", accuracy)

    accuracy = get_accuracy_from_confusion_matrix(total_modal_confusion)

    print("No Bayes", accuracy)

    plot_filename = "confusion.png"
    plot_path = os.path.join(output_dir, plot_filename)

    figure = pyplot.figure()
    axes = pyplot.axes()
    axes.set_xlabel("Predicted")
    axes.set_ylabel("True")

    pyplot.imshow(numpy.log10(total_confusion))
    pyplot.show()
    figure.savefig(plot_path)

    pyplot.close()

    plot_filename = "modal_confusion.png"
    plot_path = os.path.join(output_dir, plot_filename)

    figure = pyplot.figure()
    axes = pyplot.axes()
    axes.set_xlabel("Predicted")
    axes.set_ylabel("True")

    pyplot.imshow(numpy.log10(total_modal_confusion))
    pyplot.show()
    figure.savefig(plot_path)

    pyplot.close()


if __name__ == "__main__":
    main()
