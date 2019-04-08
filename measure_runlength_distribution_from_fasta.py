from discrete_weibull_distribution import evaluate_discrete_weibull, calculate_mode
from handlers.RunlengthHandler import RunlengthHandler
from handlers.FileManager import FileManager
from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from modules.align import align_minimap
from modules.matrix import *
import numpy
import os

MAX_RUNLENGTH = 50

# Indexes for "observation" data tuples
OBS_BASE = 0
OBS_REVERSAL = 1
OBS_LENGTH = 2
OBS_WEIGHT = 3

# Dimension labels for matrix of runlengths
MATRIX_REVERSAL = 0
MATRIX_BASE = 1
MATRIX_TRUE_LENGTH = 2
MATRIX_OBSERVED_LENGTH = 3

# Indexing RLE tuples
SEQUENCE = 0
LENGTHS = 1


def update_frequency_matrix(observed_length, true_base, true_length, alignment_reversal, matrix):
    # print(alignment_reversal, int(alignment_reversal))

    # did alignment reverse complement the sequence (via BAM) and ref (via pysam)? if so, revert to forward direction
    # if alignment_reversal:
    #     true_base = complement_base(true_base)

    # observed_base_index = BASE_TO_INDEX[observed_base]
    true_base_index = BASE_TO_INDEX[true_base]

    observed_length = min(MAX_RUNLENGTH, observed_length)
    true_length = min(MAX_RUNLENGTH, true_length)

    matrix[int(alignment_reversal), true_base_index, true_length, observed_length] += 1

    return


def complement_base(base):
    complement = None

    if base == "A":
        complement = "T"
    elif base == "T":
        complement = "A"
    elif base == "C":
        complement = "G"
    elif base == "G":
        complement = "C"
    else:
        exit("ERROR: invalid base has no complement: %s" % base)

    return complement


def reverse_complement_runlength_read(runlength_read):
    sequence = runlength_read[SEQUENCE]
    lengths = runlength_read[LENGTHS]

    reverse_complemented_sequence = list()
    reverse_complemented_lengths = list()

    for i in reversed(range(len(sequence))):
        reverse_complemented_sequence.append(complement_base(sequence[i]))
        reverse_complemented_lengths.append(lengths[i])

    reverse_complemented_sequence = "".join(reverse_complemented_sequence)

    return reverse_complemented_sequence, reverse_complemented_lengths


def parse_match(alignment_position, length, read_sequence, observed_lengths, ref_sequence, ref_runlengths, reversal_status, matrix):
    """
    Process a cigar operation that is a match
    :param alignment_position: Position where this match happened
    :param read_sequence: Read sequence
    :param ref_sequence: Reference sequence
    :param length: Length of the operation
    :return:

    This method updates the candidates dictionary.
    """
    start = alignment_position
    stop = start + length

    # print(len(ref_sequence))
    # print(len(read_sequence))
    # print(len(observed_lengths))

    for i in range(start, stop):
        allele = read_sequence[i-alignment_position]
        ref_base = ref_sequence[i-alignment_position]
        ref_runlength = ref_runlengths[i-alignment_position]

        try:
            observed_runlength = observed_lengths[i-alignment_position]

            update_frequency_matrix(observed_length=observed_runlength,
                                    true_base=ref_base,
                                    true_length=ref_runlength,
                                    alignment_reversal=reversal_status,
                                    matrix=matrix)

            # print(len(ref_sequence))
            # print(len(read_sequence))
            # print(len(observed_lengths))
            # print(i-alignment_position)
            # print(observed_runlength, ref_runlength)

        except IndexError as e:
            print()
            print(e)
            print(len(ref_sequence))
            print(len(read_sequence))
            print(len(observed_lengths))
            print(i-alignment_position)
            exit()

        if allele != ref_base:
            # mismatch
            pass

        else:
            # match
            pass

    return


def parse_cigar_tuple(cigar_code, length, alignment_position, read_sequence_segment, ref_sequence_segment,
                      ref_runlengths_segment, observed_lengths_segment, reversal_status, matrix):
    """
    Read an RLE reference segment, a segment of an assembly, and all the aligned reads that contributed to that assembly
    in order to generate a matrix of true vs observed runlengths. FOR NOW using only matches.

    cigar key map based on operation.
    details: http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    0: "MATCH",
    1: "INSERT",
    2: "DELETE",
    3: "REFSKIP",
    4: "SOFTCLIP",
    5: "HARDCLIP",
    6: "PAD"

    :param cigar_code:
    :param length:
    :param alignment_position:
    :param read_sequence_segment:
    :param ref_sequence_segment:
    :param ref_runlengths_segment:
    :param observed_lengths_segment:
    :param matrix:
    :return:
    """
    ref_index_increment = length
    read_index_increment = length

    if cigar_code == 0:
        # match
        parse_match(alignment_position=alignment_position,
                    length=length,
                    read_sequence=read_sequence_segment,
                    observed_lengths=observed_lengths_segment,
                    ref_sequence=ref_sequence_segment,
                    ref_runlengths=ref_runlengths_segment,
                    reversal_status=reversal_status,
                    matrix=matrix)

    elif cigar_code == 1:
        # insert
        ref_index_increment = 0

    elif cigar_code == 2 or cigar_code == 3:
        # delete or ref_skip
        read_index_increment = 0

    elif cigar_code == 4:
        # soft clip
        ref_index_increment = 0

    elif cigar_code == 5:
        # hard clip
        ref_index_increment = 0
        read_index_increment = 0

    elif cigar_code == 6:
        # pad
        ref_index_increment = 0
        read_index_increment = 0

    else:
        raise ("INVALID CIGAR CODE: %s" % cigar_code)

    return ref_index_increment, read_index_increment


def get_read_stop_position(read):
    """
    Returns the stop position of the reference to where the read stops aligning
    :param read: The read (pysam object)
    :return: stop position of the reference where the read last aligned
    """
    ref_alignment_stop = read.reference_end

    # only find the position if the reference end is fetched as none from pysam API
    # From pysam docs:
    #   "None values will be included for any soft-clipped or unaligned positions within the read."
    if ref_alignment_stop is None:
        positions = read.get_reference_positions()

        # find last entry that isn't None
        i = len(positions) - 1
        ref_alignment_stop = positions[-1]
        while i > 0 and ref_alignment_stop is None:
            i -= 1
            ref_alignment_stop = positions[i]

    return ref_alignment_stop


def parse_reads(reads, chromosome_name, fasta_handler, runlength_read_sequences, complete_ref_runlengths, matrix):
    """
    Given a set of pysam read objects, generate data for matches/mismatches/inserts/deletes and contig size/position for
    each read
    :param reads: pysam aligned segment objects
    :param chromosome_name:
    :param fasta_handler: fasta_handler object that can retrieve substrings from the reference sequence
    :return:
    """
    read_data = list()

    for read in reads:
        if read.mapping_quality > 0:
            read_id = read.query_name
            ref_alignment_start = read.reference_start
            ref_alignment_stop = get_read_stop_position(read)
            ref_length = ref_alignment_stop - ref_alignment_start
            cigar_tuples = read.cigartuples
            reversal_status = read.is_reverse
            runlength_read = runlength_read_sequences[read_id]

            if reversal_status:
                runlength_read = reverse_complement_runlength_read(runlength_read)

            observed_lengths = runlength_read[LENGTHS]
            read_sequence = runlength_read[SEQUENCE]

            n_initial_hard_clipped_bases = 0
            if cigar_tuples[0][0] == 5:
                n_initial_hard_clipped_bases = cigar_tuples[0][1]

            n_final_hard_clipped_bases = 0
            if cigar_tuples[-1][0] == 5:
                n_final_hard_clipped_bases = cigar_tuples[-1][1]

            clipped_start = n_initial_hard_clipped_bases
            clipped_stop = len(observed_lengths) - n_final_hard_clipped_bases

            read_sequence = read_sequence[clipped_start:clipped_stop]
            observed_lengths = observed_lengths[clipped_start:clipped_stop]

            ref_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                      start=ref_alignment_start,
                                                      stop=ref_alignment_stop + 10)

            ref_runlengths = complete_ref_runlengths[ref_alignment_start:ref_alignment_stop + 10]

            # read_quality = read.query_qualities
            # read_index: index of read sequence
            # ref_index: index of reference sequence
            read_index = 0
            ref_index = 0
            found_valid_cigar = False

            for c, cigar in enumerate(cigar_tuples):
                cigar_code = cigar[0]
                length = cigar[1]

                # Stop looping if encountered terminal hard/softclips
                if c == len(cigar_tuples) - 1 and (cigar_code == 5 or cigar_code == 4):
                    break

                # get the sequence segments that are affected by this operation
                read_sequence_segment = read_sequence[read_index:read_index+length]
                observed_lengths_segment = observed_lengths[read_index:read_index+length]
                ref_sequence_segment = ref_sequence[ref_index:ref_index+length]
                ref_runlengths_segment = ref_runlengths[ref_index:ref_index+length]

                if len(observed_lengths_segment) == 0:
                    print(len(read_sequence), read_index, read_index + length, length)

                # skip parsing the first segment if it is not a match
                if cigar_code != 0 and found_valid_cigar is False:
                    # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                    if cigar_code == 1 or cigar_code == 4:
                        read_index += length
                    continue

                found_valid_cigar = True

                # send the cigar tuple to get attributes we got by this operation
                ref_index_increment, read_index_increment = parse_cigar_tuple(cigar_code=cigar_code,
                                                                              length=length,
                                                                              alignment_position=ref_alignment_start + ref_index,
                                                                              read_sequence_segment=read_sequence_segment,
                                                                              ref_sequence_segment=ref_sequence_segment,
                                                                              ref_runlengths_segment=ref_runlengths_segment,
                                                                              observed_lengths_segment=observed_lengths_segment,
                                                                              reversal_status=reversal_status,
                                                                              matrix=matrix)

                # increase the read index iterator
                read_index += read_index_increment
                ref_index += ref_index_increment

    return read_data


def generate_runlength_frequency_matrix(runlength_ref_sequence_path, read_vs_ref_bam_path,
                                        runlength_ref_sequences, runlength_read_sequences):
    """
    Take an alignment of RLE sequences (in BAM format, using minimap as an aligner) in combination with the series of
    lengths (which have been excluded from the BAM) and aligned observations from Benedicts' model to generate a matrix
    of true vs observed lengths.
    :param runlength_ref_sequence_path:
    :param read_vs_ref_bam_path:
    :return:
    """

    for chromosome_name in runlength_ref_sequences:
        shape = [2,4,MAX_RUNLENGTH+1,MAX_RUNLENGTH+1]
        matrix = numpy.zeros(shape, dtype=numpy.float64)

        # print(read_vs_ref_bam_path)

        bam_handler = BamHandler(read_vs_ref_bam_path)
        fasta_handler = FastaHandler(runlength_ref_sequence_path)

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=0, stop=chromosome_length)

        parse_reads(chromosome_name=chromosome_name,
                    fasta_handler=fasta_handler,
                    reads=reads,
                    complete_ref_runlengths=runlength_ref_sequences[chromosome_name][LENGTHS],
                    runlength_read_sequences=runlength_read_sequences,
                    matrix=matrix)

        # plot_base_matrices(matrix=matrix, cutoff=40)

        yield matrix


def runlength_encode(sequence):
    character_sequence = list()
    character_counts = list()
    current_character = None

    for character in sequence:
        if character != current_character:
            character_sequence.append(character)
            character_counts.append(1)
        else:
            character_counts[-1] += 1

        current_character = character

    character_sequence = ''.join(character_sequence)

    return character_sequence, character_counts


def runlength_encode_fasta(fasta_sequence_path):
    fasta_handler = FastaHandler(fasta_sequence_path)

    contig_names = fasta_handler.get_contig_names()

    runlength_sequences = dict()

    for contig_name in contig_names:
        chromosome_length = fasta_handler.get_chr_sequence_length(contig_name)

        sequence = fasta_handler.get_sequence(chromosome_name=contig_name, start=None, stop=None)

        bases, lengths = runlength_encode(sequence)

        runlength_sequences[contig_name] = (bases, lengths)

        # print(contig_name, len(bases), len(lengths))

    return runlength_sequences


def iteratively_align_as_RLE(ref_fasta_path, read_fasta_path, output_dir):
    """
    Given 2 fasta files for reads and reference, iterate them, runlength encode their sequences, and write the RLE
    sequences to a new file, then align them with minimap2
    :param ref_fasta_path:
    :param read_fasta_path:
    :param output_dir:
    :return:
    """
    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    read_fasta_filename_prefix = ".".join(os.path.basename(read_fasta_path).split(".")[:-1])
    runlength_read_fasta_filename = read_fasta_filename_prefix + "_rle.fasta"
    runlength_read_fasta_path = os.path.join(output_dir, runlength_read_fasta_filename)

    print("SAVING run length fasta file:", runlength_ref_fasta_path)
    print("SAVING run length fasta file:", runlength_read_fasta_path)

    with open(runlength_ref_fasta_path, "w") as file:
        fasta_handler = FastaHandler(ref_fasta_path)
        names = fasta_handler.get_contig_names()

        for name in names:
            sequence = fasta_handler.get_sequence(chromosome_name=name, start=None, stop=None)
            sequence, lengths = runlength_encode(sequence)

            file.write(">"+name+" RLE\n")
            file.write(sequence + "\n")

    with open(runlength_read_fasta_path, "w") as file:
        fasta_handler = FastaHandler(read_fasta_path)
        names = fasta_handler.get_contig_names()

        for name in names:
            sequence = fasta_handler.get_sequence(chromosome_name=name, start=None, stop=None)
            sequence, lengths = runlength_encode(sequence)

            file.write(">"+name+" RLE\n")
            file.write(sequence + "\n")

    output_sam_file_path, output_bam_file_path = align_minimap(output_dir=output_dir,
                                                               ref_sequence_path=runlength_ref_fasta_path,
                                                               reads_sequence_path=runlength_read_fasta_path)

    return output_bam_file_path


def align_as_RLE(runlength_reference_path, runlength_ref_sequences, runlength_read_path, runlength_read_sequences, output_dir):
    print("SAVING run length fasta file:", runlength_reference_path)
    print("SAVING run length fasta file:", runlength_read_path)

    with open(runlength_reference_path, "w") as file:
        for contig_name in runlength_ref_sequences:
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_ref_sequences[contig_name][SEQUENCE] + "\n")

    with open(runlength_read_path, "w") as file:
        for contig_name in runlength_read_sequences:
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_read_sequences[contig_name][SEQUENCE] + "\n")

    output_sam_file_path, output_bam_file_path = align_minimap(output_dir=output_dir,
                                                               ref_sequence_path=runlength_reference_path,
                                                               reads_sequence_path=runlength_read_path,
                                                               k=18)

    return output_bam_file_path


def main():
    ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    read_fasta_path = "/home/ryan/code/runlength_analysis/output/guppy_vs_runnie_ecoli_rad2_train_test_sequences/sequence_subset_train_60x_10kb.fastq"
    # read_fasta_path = "/home/ryan/software/shasta/output/run_2019_4_2_9_28_TEST_BAYES/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/run_2019_4_2_9_29_TEST_NO_BAYES/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/run_2019_4_2_13_34_BAYES_NO_PRIOR/Assembly.fasta"

    # read_fasta_path = "/home/ryan/software/shasta/output/NO_BAYES_4_4_2019/run_2019_4_4_16_22_2_217698/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_DEFAULT_4_4_2019/run_2019_4_4_16_23_46_269190/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_NEW_4_4_2019/run_2019_4_4_16_26_22_687430/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_STRAND_4_4_2019/run_2019_4_5_12_58_42_872166/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_STRAND_NO_PRIOR_4_4_2019/run_2019_4_5_13_22_31_842825/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_STRAND_4_5_2019/run_2019_4_5_13_46_19_628/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/BAYES_STRAND_NO_PRIOR_4_5_2019/run_2019_4_5_13_58_30_740520/Assembly.fasta"

    # read_fasta_path = "/home/ryan/software/shasta/output/run_2019_3_28_15_28_11_844359/Assembly.fasta"
    # read_fasta_path = "/home/ryan/software/shasta/output/run_2019_3_28_12_55_NO_BAYES/Assembly.fasta"

    # ---- TEST DATA ----
    # ref_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_runlength_test_2019_3_27_16_34_11_810671_ref.fasta"
    # read_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_runlength_test_2019_3_27_16_34_11_810671_reads.fasta"
    # -------------------

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

    chromosomal_matrices = generate_runlength_frequency_matrix(runlength_ref_sequence_path=runlength_ref_fasta_path,
                                                               read_vs_ref_bam_path=read_vs_ref_bam_path,
                                                               runlength_ref_sequences=runlength_ref_sequences,
                                                               runlength_read_sequences=runlength_read_sequences)

    for matrix in chromosomal_matrices:
        save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                              frequency_matrices=matrix,
                                                              log_normalize=False,
                                                              plot=False)

        save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                              frequency_matrices=matrix,
                                                              log_normalize=True,
                                                              plot=False)

        nondirectional_matrix = sum_complementary_matrices(matrix, max_runlength=MAX_RUNLENGTH)

        save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                                 frequency_matrices=nondirectional_matrix,
                                                                 log_normalize=False,
                                                                 plot=False)

        save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                                 frequency_matrices=nondirectional_matrix,
                                                                 log_normalize=True,
                                                                 plot=False)

        # zero_mask = (matrix == 0)
        # nonzero_mask = numpy.invert(zero_mask)
        # matrix[zero_mask] += numpy.min(matrix[nonzero_mask])

        plot_directional_complementary_diff(matrix, max_runlength=MAX_RUNLENGTH)
        plot_base_matrices(matrix, test_spot=False, normalize_matrices=False)
        plot_base_matrices(matrix, test_spot=False, normalize_matrices=True)


if __name__ == "__main__":
    main()
