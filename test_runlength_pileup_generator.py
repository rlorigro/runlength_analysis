from measure_runlength_distribution_from_runnie import runlength_encode_fasta, align_as_RLE
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


def map_parameters_to_mode(scale_shape_tuple, max_runlength=50):
    x = numpy.arange(0,max_runlength)

    mode = RunlengthHandler.calculate_numerical_mode(scale=scale_shape_tuple[0], shape=scale_shape_tuple[1], x_range=x)

    return mode


def convert_parameters_to_modes(scales, shapes):
    modes = list()
    for read_id in scales.keys():
        modes.append(list(map(map_parameters_to_mode, (scales[read_id], shapes[read_id]))))

    return modes


def plot_runlength_pileup(sequences, scales, shapes, modes):
    figure, axes = pyplot.subplots(nrows=4)

    axes[0].imshow(sequences)
    axes[1].imshow(scales)
    axes[2].imshow(shapes)
    axes[3].imshow(modes)

    axes[0].set_ylabel("Nucleotide")
    axes[1].set_ylabel("Scale (alpha)")
    axes[2].set_ylabel("Shape (beta)")
    axes[3].set_ylabel("Mode")

    pyplot.show()
    pyplot.close()


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


def get_aligned_segments(fasta_handler, bam_handler, chromosome_name, pileup_start, pileup_end, read_data, include_ref=False):
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
                                      read_data=read_data,
                                      reads=reads)

    # if a reference sequence is intended to be added to the pileup, then leave a space for it
    if include_ref:
        segment_grabber.max_coverage -= 1

    sequence_dictionary = segment_grabber.get_aligned_read_segments()

    return sequence_dictionary


def main():
    ref_fasta_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_ref.fasta"
    runlength_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_18_11_56_2_830712_runnie.out"

    pileup_start = 0
    pileup_end = 127

    output_parent_dir = "output/"
    output_dir = "runlength_matrix_from_runnie_output_" + FileManager.get_datetime_string()
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

    aligned_sequences, aligned_scales, aligned_shapes = get_aligned_segments(fasta_handler=fasta_handler,
                                                                             bam_handler=bam_handler,
                                                                             chromosome_name=chromosome_name,
                                                                             pileup_start=pileup_start,
                                                                             pileup_end=pileup_end,
                                                                             include_ref=True,
                                                                             read_data=read_data)

    sequence_encoding = list()
    scale_encoding = list()
    shape_encoding = list()
    modes_encoding = list()
    for read_id in aligned_sequences.keys():
        sequence_encoding.append(list(map(get_encoding, aligned_sequences[read_id])))
        scale_encoding.append(aligned_scales[read_id])
        shape_encoding.append(aligned_shapes[read_id])
        modes_encoding.append(list(map(map_parameters_to_mode, zip(aligned_scales[read_id], aligned_shapes[read_id]))))

    sequence_encoding = -numpy.array(sequence_encoding, dtype=numpy.float)
    scale_encoding = numpy.array(scale_encoding, dtype=numpy.float)
    shape_encoding = numpy.array(shape_encoding, dtype=numpy.float)
    modes_encoding = numpy.array(modes_encoding, dtype=numpy.float)

    plot_runlength_pileup(sequences=sequence_encoding,
                          scales=scale_encoding,
                          shapes=shape_encoding,
                          modes=modes_encoding)


if __name__ == "__main__":
    main()
