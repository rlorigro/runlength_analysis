from measure_runlength_distribution_from_runnie import runlength_encode_fasta, align_as_RLE
from modules.WeibullRunlengthClassifier import RunlengthClassifier as WeibullRunlengthClassifier
from modules.RunlengthClassifier import RunlengthClassifier
from modules.RunniePileupGenerator import INSERT_CHAR, DELETE_CHAR
from modules.RunniePileupGenerator import PileupGenerator
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


def plot_runlength_pileup(sequences, scales, shapes, modes, ref_sequence, predicted_sequence, ref_lengths, predicted_lengths):
    """
    Given 2d numpy float matrices, make an organized plot of pileup data
    :param sequences:
    :param scales:
    :param shapes:
    :param modes:
    :param ref_sequence:
    :param ref_lengths:
    :return:
    """
    figure, axes = pyplot.subplots(nrows=8, sharex=True)

    axes[0].imshow(ref_sequence)
    axes[1].imshow(predicted_sequence)
    axes[2].imshow(sequences)
    axes[3].imshow(scales)
    axes[4].imshow(shapes)
    axes[5].imshow(ref_lengths)
    axes[6].imshow(predicted_lengths)
    axes[7].imshow(modes)

    axes[0].set_ylabel("Ref Nucleotide")
    axes[1].set_ylabel("Predicted nt")
    axes[2].set_ylabel("Read Nucleotide")
    axes[3].set_ylabel("Scale (alpha)")
    axes[4].set_ylabel("Shape (beta)")
    axes[5].set_ylabel("Ref length")
    axes[6].set_ylabel("Predicted length")
    axes[7].set_ylabel("Mode")

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


def get_consensus_from_weibull_pileup_encoding(length_classifier, sequence_encoding, scale_encoding, shape_encoding, reversal_encoding):
    # print(sequence_encoding)
    # print(length_encoding)

    window_length = sequence_encoding.shape[1]

    consensus_sequence = list()
    consensus_lengths = list()

    for i in range(window_length):
        sequence_observations = sequence_encoding[:,i]
        scale_observations = scale_encoding[:,i]
        shape_observations = shape_encoding[:,i]

        base_consensus = mode(sequence_observations)

        # print("sequence_observations\t", sequence_observations)
        # print("scale_observations\t", scale_observations)
        # print("shape_observations\t", shape_observations)
        # print("reversal_encoding\t", reversal_encoding)
        # print("base consensus\t\t", base_consensus)

        if base_consensus == BASE_TO_INDEX["-"]:
            length_consensus = 0

        else:
            normalized_y_log_likelihoods, length_consensus = length_classifier.predict(x_scales=scale_observations,
                                                                                       x_shapes=shape_observations,
                                                                                       character_index=base_consensus,
                                                                                       reversal=reversal_encoding)

        # if length_consensus == 1:
        #     print(base_consensus)
        #     print(scale_observations)
        #     print(shape_observations)

        consensus_sequence.append(base_consensus)
        consensus_lengths.append(length_consensus)

        # print("length normalized_y_log_likelihoods\n", 10**normalized_y_log_likelihoods[:10])
        # print("length consensus\t", length_consensus)
        # print()

    return consensus_sequence, consensus_lengths


def get_consensus_from_modal_pileup_encoding(length_classifier, sequence_encoding, length_encoding, reversal_encoding):
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
            normalized_y_log_likelihoods, length_consensus = length_classifier.predict(x=length_observations,
                                                                                       character_index=base_consensus,
                                                                                       reversal=reversal_encoding)

        consensus_sequence.append(base_consensus)
        consensus_lengths.append(length_consensus)

        # if length_consensus == 0:
        #     print()
        #     print(base_consensus, INDEX_TO_BASE[base_consensus])
        #     print(length_observations)
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
    # ref_fasta_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_ref.fasta"
    # runlength_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_runnie.out"

    # ref_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_runnie_test_2019_4_8_14_33_30_333396_ref.fasta"
    # runlength_path = "/home/ryan/code/runlength_analysis/data/synthetic_runnie_test_2019_4_8_14_33_30_333396_runnie.out"

    ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    runlength_path = "/home/ryan/code/runlength_analysis/data/runnie_subset_test_flipflop_regional_0to10k.out"

    # WG ecoli 60x
    matrix_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_runnie_WG_train_60x_guppy_2019_4_23/probability_matrices_2019_4_23_15_9_14_837893.csv"
    raw_matrix_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_from_runnie_WG_train_60x_guppy_2019_4_23/frequency_matrices_2019_4_23_15_9_14_833128.csv"

    output_parent_dir = "output/"
    output_dir = "runlength_prediction_from_runnie_output_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    assembly_fasta_filename_prefix = ".".join(os.path.basename(runlength_path).split(".")[:-1])
    runlength_assembly_fasta_filename = assembly_fasta_filename_prefix + "_rle.fasta"
    runlength_assembly_fasta_path = os.path.join(output_dir, runlength_assembly_fasta_filename)

    handler = RunlengthHandler(runlength_path)

    reads = handler.iterate_file(sequence_cutoff=sys.maxsize, print_status=True)
    read_data = dict()

    for r, read in enumerate(reads):
        read_data[read.id] = read

    print("\nRLE encoding reference sequence...")

    runlength_ref_sequences = runlength_encode_fasta(fasta_sequence_path=ref_fasta_path)

    assembly_vs_ref_bam_path = align_as_RLE(runlength_reference_path=runlength_ref_fasta_path,
                                            runlength_ref_sequences=runlength_ref_sequences,
                                            runlength_read_path=runlength_assembly_fasta_path,
                                            runlength_read_sequences=read_data,
                                            output_dir=output_dir)

    bam_handler = BamHandler(assembly_vs_ref_bam_path)
    fasta_handler = FastaHandler(runlength_ref_fasta_path)

    contig_names = fasta_handler.get_contig_names()
    chromosome_name = contig_names[0]
    chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

    windows = chunk_chromosome_coordinates(chromosome_length=chromosome_length, chunk_size=1000)

    total_confusion = get_runlength_confusion([],[],10)
    total_confusion_weibull = get_runlength_confusion([],[],10)

    length_classifier = RunlengthClassifier(matrix_path)
    # length_classifier_weibull = WeibullRunlengthClassifier(matrix_path)
    length_classifier_weibull = WeibullRunlengthClassifier(raw_matrix_path, normalize_matrix=True, pseudocount=0.05)

    print("reading BAM")
    for pileup_start, pileup_end in windows[10:20]:
        sys.stderr.write("\r%s" % pileup_start)
        aligned_ref_sequence, aligned_ref_lengths, aligned_sequences, aligned_scales, aligned_shapes, reversal_statuses = \
            get_aligned_segments(fasta_handler=fasta_handler,
                                 bam_handler=bam_handler,
                                 chromosome_name=chromosome_name,
                                 pileup_start=pileup_start,
                                 pileup_end=pileup_end,
                                 runlength_ref_sequences=runlength_ref_sequences,
                                 read_data=read_data)

        sequence_encoding = list()
        scale_encoding = list()
        shape_encoding = list()
        modes_encoding = list()
        reversal_encoding = list()

        # No reads here?
        if len(aligned_sequences) == 0:
            continue

        try:
            # print("REF\t", "".join(aligned_ref_sequence))
            for read_id in aligned_sequences.keys():
                # print("READ\t%s\t%s" % (read_id,"".join(aligned_sequences[read_id])))
                sequence_encoding.append(list(map(get_encoding, aligned_sequences[read_id])))
                scale_encoding.append(aligned_scales[read_id])
                shape_encoding.append(aligned_shapes[read_id])
                modes_encoding.append(list(map(map_parameters_to_mode, zip(aligned_scales[read_id], aligned_shapes[read_id]))))
                reversal_encoding.append(reversal_statuses[read_id])

            ref_sequence_encoding = [list(map(get_encoding, aligned_ref_sequence))]
            ref_lengths_encoding = [aligned_ref_lengths]

            ref_sequence_encoding = numpy.atleast_2d(numpy.array(ref_sequence_encoding, dtype=numpy.int))
            ref_length_encoding = numpy.atleast_2d(numpy.array(ref_lengths_encoding, dtype=numpy.int))
            sequence_encoding = numpy.atleast_2d(numpy.array(sequence_encoding, dtype=numpy.int))
            scale_encoding = numpy.atleast_2d(numpy.array(scale_encoding, dtype=numpy.float))
            shape_encoding = numpy.atleast_2d(numpy.array(shape_encoding, dtype=numpy.float))
            modes_encoding = numpy.atleast_2d(numpy.array(modes_encoding, dtype=numpy.int))
            reversal_encoding = numpy.array(reversal_encoding, dtype=numpy.bool)

            consensus_sequence, consensus_lengths = \
                get_consensus_from_modal_pileup_encoding(length_classifier=length_classifier,
                                                         sequence_encoding=sequence_encoding,
                                                         length_encoding=modes_encoding,
                                                         reversal_encoding=reversal_encoding)

            weibull_consensus_sequence, weibull_consensus_lengths = \
                get_consensus_from_weibull_pileup_encoding(length_classifier=length_classifier_weibull,
                                                           sequence_encoding=sequence_encoding,
                                                           scale_encoding=scale_encoding,
                                                           shape_encoding=shape_encoding,
                                                           reversal_encoding=reversal_encoding)

            plot_runlength_pileup(sequences=-sequence_encoding,
                                  scales=scale_encoding,
                                  shapes=shape_encoding,
                                  modes=modes_encoding,
                                  ref_sequence=-ref_sequence_encoding,
                                  ref_lengths=ref_length_encoding,
                                  predicted_sequence=-numpy.atleast_2d(numpy.array(weibull_consensus_sequence, dtype=numpy.int)),
                                  predicted_lengths=numpy.atleast_2d(numpy.array(weibull_consensus_lengths, dtype=numpy.int)))

            print()
            print("PREDICTED\t",weibull_consensus_lengths[:10])
            print("TRUE\t\t",aligned_ref_lengths[:10])

            confusion = get_runlength_confusion(true_lengths=aligned_ref_lengths,
                                                predicted_lengths=consensus_lengths,
                                                max_length=10)

            confusion_weibull = get_runlength_confusion(true_lengths=aligned_ref_lengths,
                                                        predicted_lengths=weibull_consensus_lengths,
                                                        max_length=10)

            total_confusion += confusion
            total_confusion_weibull += confusion_weibull

        except Exception as e:
            print(e)
            continue
    print()

    accuracy = get_accuracy_from_confusion_matrix(total_confusion)

    print("Modal: ", accuracy)

    accuracy = get_accuracy_from_confusion_matrix(total_confusion_weibull)

    print("Full: ", accuracy)

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

    plot_filename = "confusion_weibull.png"
    plot_path = os.path.join(output_dir, plot_filename)

    figure = pyplot.figure()

    axes = pyplot.axes()
    axes.set_xlabel("Predicted")
    axes.set_ylabel("True")

    pyplot.imshow(numpy.log10(total_confusion_weibull))
    pyplot.show()
    figure.savefig(plot_path)

    pyplot.close()


if __name__ == "__main__":
    main()
