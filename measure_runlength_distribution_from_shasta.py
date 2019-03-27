from discrete_weibull_distribution import evaluate_discrete_weibull, calculate_mode
from handlers.ShastaRunlengthHandler import ShastaRunlengthHandler
from handlers.FileManager import FileManager
from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from modules.align import align_minimap
from matplotlib import pyplot
from matplotlib import colors
import numpy
import math
import sys
import os

MAX_RUNLENGTH = 50
X_RANGE = numpy.arange(0, MAX_RUNLENGTH)

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

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


def print_frequency_matrices(frequency_matrices, cutoff=None):
    if cutoff is None:
        cutoff = frequency_matrices.shape[-1]

    if frequency_matrices.ndim == 4:
        for i in range(frequency_matrices.shape[1]):
            print(numpy.array2string(frequency_matrices[0,i,:cutoff,:cutoff], max_line_width=sys.maxsize, separator='\t', formatter={'int': lambda x: "%d"%x}))
    elif frequency_matrices.ndim == 3:
        for i in range(frequency_matrices.shape[0]):
            print(numpy.array2string(frequency_matrices[i,:cutoff,:cutoff], max_line_width=sys.maxsize, separator='\t', formatter={'int': lambda x: "%d"%x}))


def plot_directional_residuals(frequency_matrices, cutoff=20):
    frequency_matrices = diff_complementary_matrices(frequency_matrices)
    # print_frequency_matrices(frequency_matrices=frequency_matrices, cutoff=10)

    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    vmax = numpy.max(frequency_matrices)
    vmin = -vmax

    colormap = pyplot.cm.RdBu
    norm = colors.SymLogNorm(linthresh=0.03, vmin=vmin, vmax=vmax)

    c_axes = axes[0][0].imshow(frequency_matrices[0, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][0].set_title(INDEX_TO_BASE[0])
    figure.colorbar(c_axes, ax=axes[0][0])

    c_axes = axes[1][0].imshow(frequency_matrices[3, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][0].set_title(INDEX_TO_BASE[3])
    figure.colorbar(c_axes, ax=axes[1][0])

    c_axes = axes[0][1].imshow(frequency_matrices[2, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][1].set_title(INDEX_TO_BASE[2])
    figure.colorbar(c_axes, ax=axes[0][1])

    c_axes = axes[1][1].imshow(frequency_matrices[1, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][1].set_title(INDEX_TO_BASE[1])
    figure.colorbar(c_axes, ax=axes[1][1])

    pyplot.show()
    pyplot.close()


def diff_complementary_matrices(directional_matrices):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by subtracting the counts of any reverse bases from their base complement's matrix
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, MAX_RUNLENGTH+1, MAX_RUNLENGTH+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        if r == 0:
            reversal = False
        else:
            reversal = True

        for b in range(n_bases):
            if reversal:
                destination_index = complement_base_index(b)
                print(nondirectional_matrices[destination_index,:,:].shape, directional_matrices[r,b,:,:].shape)
                nondirectional_matrices[destination_index,:,:] -= directional_matrices[r,b,:,:]
            else:
                destination_index = b
                print(nondirectional_matrices[destination_index,:,:].shape, directional_matrices[r,b,:,:].shape)
                nondirectional_matrices[destination_index,:,:] = directional_matrices[r,b,:,:]

    return nondirectional_matrices


def plot_base_matrices(matrix, test_spot=False, cutoff=20, normalize_matrices=False):
    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    if matrix.ndim == 4:
        matrix_A = matrix[0, 0, :cutoff, :cutoff]

        if test_spot:
            matrix_A[cutoff-1,0] = numpy.max(numpy.log10(matrix_A)) + 1

        if normalize_matrices:
            matrix_A = normalize(matrix_A, pseudocount=0)
            axes[0][0].imshow(matrix_A)
        else:
            axes[0][0].imshow(numpy.log10(matrix_A))

        axes[0][0].set_title(INDEX_TO_BASE[0])

        matrix_T = matrix[0, 3, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_T  = normalize(matrix_T, pseudocount=0)
            axes[1][0].imshow(matrix_T)
        else:
            axes[1][0].imshow(numpy.log10(matrix_T))

        axes[1][0].set_title(INDEX_TO_BASE[3])

        matrix_G = matrix[0, 2, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_G  = normalize(matrix_G, pseudocount=0)
            axes[0][1].imshow(matrix_G)
        else:
            axes[0][1].imshow(numpy.log10(matrix_G))

        axes[0][1].set_title(INDEX_TO_BASE[2])

        matrix_C = matrix[0, 1, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_C  = normalize(matrix_C, pseudocount=0)
            axes[1][1].imshow(matrix_C)
        else:
            axes[1][1].imshow(numpy.log10(matrix_C))

        axes[1][1].set_title(INDEX_TO_BASE[1])

    elif matrix.ndim == 3:
        matrix_A = matrix[0, :cutoff, :cutoff]

        if test_spot:
            matrix_A[cutoff, 0] = numpy.max(matrix_A)

        if normalize_matrices:
            matrix_A = normalize(matrix_A, pseudocount=0)
            axes[0][0].imshow(matrix_A)
        else:
            axes[0][0].imshow(numpy.log10(matrix_A))

        axes[0][0].set_title(INDEX_TO_BASE[0])

        matrix_T = matrix[3, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_T  = normalize(matrix_T, pseudocount=0)
            axes[1][0].imshow(matrix_T)
        else:
            axes[1][0].imshow(numpy.log10(matrix_T))
        axes[1][0].set_title(INDEX_TO_BASE[3])

        matrix_G = matrix[2, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_G  = normalize(matrix_G, pseudocount=0)
            axes[0][1].imshow(matrix_G)
        else:
            axes[0][1].imshow(numpy.log10(matrix_G))
        axes[0][1].set_title(INDEX_TO_BASE[2])

        matrix_C = matrix[1, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_C  = normalize(matrix_C, pseudocount=0)
            axes[1][1].imshow(matrix_C)
        else:
            axes[1][1].imshow(numpy.log10(matrix_C))
        axes[1][1].set_title(INDEX_TO_BASE[1])

    axes[1][1].set_xlabel("Observed length")
    axes[1][0].set_xlabel("Observed length")
    axes[1][0].set_ylabel("True length")
    axes[0][0].set_ylabel("True length")

    axes[1][1].set_yticks(numpy.arange(0,cutoff, 2))
    axes[1][0].set_yticks(numpy.arange(0,cutoff, 2))
    axes[0][0].set_yticks(numpy.arange(0,cutoff, 2))
    axes[0][1].set_yticks(numpy.arange(0,cutoff, 2))

    axes[1][1].set_xticks(numpy.arange(0, cutoff, 2))
    axes[1][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][1].set_xticks(numpy.arange(0, cutoff, 2))

    pyplot.show()
    pyplot.close()


def sum_complementary_matrices(directional_matrices):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by adding the counts of any reverse bases into their base complement's matrix
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, MAX_RUNLENGTH+1, MAX_RUNLENGTH+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        if r == 0:
            reversal = False
        else:
            reversal = True

        for b in range(n_bases):
            destination_index = b
            if reversal:
                destination_index = complement_base_index(b)

            nondirectional_matrices[destination_index,:,:] += directional_matrices[r,b,:,:].squeeze()

    return nondirectional_matrices


def complement_base_index(base_index):
    return 3 - base_index


def normalize(frequency_matrix, pseudocount, diagonal_bias):
    frequency_matrix = frequency_matrix.astype(numpy.float32)

    frequency_matrix += pseudocount

    diagonal_pseudocount = diagonal_bias
    diagonal_mask = numpy.eye(frequency_matrix.shape[0], frequency_matrix.shape[1], dtype=numpy.bool)
    frequency_matrix[diagonal_mask] += diagonal_pseudocount

    sum_y = numpy.sum(frequency_matrix, axis=1)

    probability_matrix = frequency_matrix / sum_y[:, numpy.newaxis]
    probability_matrix = numpy.log10(probability_matrix)

    return probability_matrix


def save_directional_frequency_matrices_as_delimited_text(output_dir, frequency_matrices, delimiter=",", log_normalize=False, plot=False, pseudocount=0, diagonal_bias=0, default_type=int):
    if log_normalize:
        filename = "probability_matrices_directional_" + FileManager.get_datetime_string() + ".csv"
    else:
        filename = "frequency_matrices_directional_" + FileManager.get_datetime_string() + ".csv"

    reversal_suffixes = ["F", "R"]
    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    for reversal in [0,1]:
        for base_index in range(4):
            base = INDEX_TO_BASE[base_index]
            suffix = reversal_suffixes[reversal]

            matrix = numpy.squeeze(frequency_matrices[reversal,base_index,:,:])

            type = default_type
            if log_normalize:
                matrix = normalize(matrix, pseudocount=pseudocount, diagonal_bias=diagonal_bias)
                type = float

            if plot:
                pyplot.imshow(matrix)
                pyplot.show()
                pyplot.close()

            matrix_name = "_".join([base, suffix])
            header = ">" + matrix_name + "\n"

            file.write(header)
            for r in range(matrix.shape[0]):
                row = [str(type(x)) for x in matrix[r]]

                row = delimiter.join(row) + "\n"

                file.write(row)

            file.write("\n")

    file.close()


def save_nondirectional_frequency_matrices_as_delimited_text(output_dir, frequency_matrices, delimiter=",", log_normalize=False, pseudocount=0, diagonal_bias=0, plot=False):
    if log_normalize:
        filename = "probability_matrices_" + FileManager.get_datetime_string() + ".csv"
    else:
        filename = "frequency_matrices_" + FileManager.get_datetime_string() + ".csv"

    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    for base_index in range(4):
        base = INDEX_TO_BASE[base_index]

        matrix = numpy.squeeze(frequency_matrices[base_index,:,:])

        type = int
        if log_normalize:
            matrix = normalize(matrix, pseudocount=15)
            type = float

        if plot:
            pyplot.imshow(matrix)
            pyplot.show()
            pyplot.close()

        matrix_name = base
        header = ">" + matrix_name + "\n"

        file.write(header)

        for r in range(matrix.shape[0]):
            row = [str(type(x)) for x in matrix[r]]
            row = delimiter.join(row) + "\n"

            file.write(row)

        file.write("\n")

    file.close()


def update_frequency_matrix(observed_pileup, consensus_base, true_base, true_length, alignment_reversal, matrix):
    # Did alignment reverse complement the sequence (via BAM) and ref (via pysam)? if so, revert to forward direction
    if alignment_reversal:
        true_base = complement_base(true_base)

    true_base_index = BASE_TO_INDEX[true_base]

    for item in observed_pileup:
        observed_length = item.length
        observed_base = item.base
        count = item.count

        if alignment_reversal:
            observed_base = complement_base(observed_base)

        if observed_base == consensus_base:
            matrix[int(alignment_reversal), true_base_index, true_length, observed_length] += count

        # matrix[int(alignment_reversal), true_base_index, true_length, observed_length] += count

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
    elif base == "_":
        complement = "_"
    else:
        exit("ERROR: invalid base has no complement: %s" % base)

    return complement


def reverse_complement_runlength_read(runlength_read):
    sequence = runlength_read.sequence
    pileup = runlength_read.pileup

    reverse_complemented_sequence = list()
    reverse_complemented_pileup = list()

    for i in reversed(range(len(sequence))):
        reverse_complemented_sequence.append(complement_base(sequence[i]))
        reverse_complemented_pileup.append(pileup[i])

    reverse_complemented_sequence = "".join(reverse_complemented_sequence)

    runlength_read.sequence = reverse_complemented_sequence
    runlength_read.pileup = reverse_complemented_pileup

    return runlength_read


def parse_match(alignment_position, length, read_sequence, observed_pileup_segment,
                ref_sequence, ref_runlengths, reversal_status, matrix):
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

    for i in range(start, stop):
        ref_base = ref_sequence[i-alignment_position]
        read_base = read_sequence[i-alignment_position]

        ref_runlength = ref_runlengths[i-alignment_position]

        observed_pileup = observed_pileup_segment[i-alignment_position]

        update_frequency_matrix(observed_pileup=observed_pileup,
                                consensus_base=read_base,
                                true_base=ref_base,
                                true_length=ref_runlength,
                                alignment_reversal=reversal_status,
                                matrix=matrix)

    return


def parse_cigar_tuple(cigar_code, length, alignment_position, read_sequence_segment, ref_sequence_segment,
                      ref_runlengths_segment, observed_pileup_segment, reversal_status, matrix):
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
                    observed_pileup_segment=observed_pileup_segment,
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


def parse_reads(reads, chromosome_name, fasta_handler, runlength_read_data, complete_ref_runlengths, matrix):
    """
    Given a set of pysam read objects, generate data for matches/mismatches/inserts/deletes and contig size/position for
    each read
    :param reads: pysam aligned segment objects
    :param chromosome_name:
    :param fasta_handler: fasta_handler object that can retrieve substrings from the reference sequence
    :return:
    """
    read_data = list()

    for r,read in enumerate(reads):
        if read.mapping_quality > 0 and not read.is_secondary:
            read_id = read.query_name
            ref_alignment_start = read.reference_start
            ref_alignment_stop = get_read_stop_position(read)
            ref_length = ref_alignment_stop - ref_alignment_start

            reversal_status = read.is_reverse

            runlength_read = runlength_read_data[read_id]

            if reversal_status:
                runlength_read = reverse_complement_runlength_read(runlength_read)

            observed_pileup = runlength_read.pileup

            ref_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                      start=ref_alignment_start,
                                                      stop=ref_alignment_stop + 10)

            ref_runlengths = complete_ref_runlengths[ref_alignment_start:ref_alignment_stop + 10]

            cigar_tuples = read.cigartuples
            read_sequence = read.query_sequence

            # read_length = len(read_sequence)
            # contig_length = read.infer_read_length()
            # read_quality = read.query_qualities

            # read_index: index of read sequence
            # ref_index: index of reference sequence
            read_index = 0
            ref_index = 0
            found_valid_cigar = False

            n_initial_clipped_bases = 0

            for c, cigar in enumerate(cigar_tuples):
                cigar_code = cigar[0]
                length = cigar[1]

                # get the sequence segments that are affected by this operation
                read_sequence_segment = read_sequence[read_index:read_index+length]
                observed_pileup_segment = observed_pileup[read_index:read_index+length]
                ref_sequence_segment = ref_sequence[ref_index:ref_index+length]
                ref_runlengths_segment = ref_runlengths[ref_index:ref_index+length]

                if len(observed_pileup_segment) == 0:
                    print(len(read_sequence), read_index, read_index + length, length)

                # skip parsing the first segment if it is not a match
                if cigar_code != 0 and found_valid_cigar is False:
                    # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                    if cigar_code == 1 or cigar_code == 4:
                        read_index += length
                    if cigar_code == 5 or cigar_code == 4:
                        n_initial_clipped_bases = length
                    continue

                found_valid_cigar = True

                ref_index_increment, read_index_increment = parse_cigar_tuple(cigar_code=cigar_code,
                                                                              length=length,
                                                                              alignment_position=ref_alignment_start + ref_index,
                                                                              read_sequence_segment=read_sequence_segment,
                                                                              ref_sequence_segment=ref_sequence_segment,
                                                                              ref_runlengths_segment=ref_runlengths_segment,
                                                                              observed_pileup_segment=observed_pileup_segment,
                                                                              reversal_status=reversal_status,
                                                                              matrix=matrix)

                # increase the read index iterator
                read_index += read_index_increment
                ref_index += ref_index_increment

    return read_data


def generate_runlength_frequency_matrix(runlength_ref_sequence_path, reads_vs_ref_bam_path,
                                        runlength_ref_sequences, runlength_read_data):
    """
    Take an alignment of RLE sequences (in BAM format, using minimap as an aligner) in combination with the series of
    lengths (which have been excluded from the BAM) and aligned observations from Benedicts' model to generate a matrix
    of true vs observed lengths.
    :param runlength_ref_sequence_path:
    :param reads_vs_ref_bam_path:
    :param runlength_ref_sequences:
    :param runlength_read_data:
    :return:
    """
    for chromosome_name in runlength_ref_sequences:
        shape = [2,4,MAX_RUNLENGTH+1,MAX_RUNLENGTH+1]
        matrix = numpy.zeros(shape, dtype=numpy.float64)

        bam_handler = BamHandler(reads_vs_ref_bam_path)
        fasta_handler = FastaHandler(runlength_ref_sequence_path)

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=0, stop=chromosome_length)

        parse_reads(chromosome_name=chromosome_name,
                    fasta_handler=fasta_handler,
                    reads=reads,
                    complete_ref_runlengths=runlength_ref_sequences[chromosome_name][LENGTHS],
                    runlength_read_data=runlength_read_data,
                    matrix=matrix)

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

        sequence = fasta_handler.get_sequence(chromosome_name=contig_name, start=0, stop=chromosome_length)

        bases, lengths = runlength_encode(sequence)

        runlength_sequences[contig_name] = (bases, lengths)

        print(contig_name, len(bases), len(lengths))

    return runlength_sequences


def align_as_RLE(runlength_reference_path, runlength_ref_sequences, runlength_reads_path, runlength_reads_sequences, output_dir):
    print("SAVING run length fasta file:", runlength_reference_path)
    print("SAVING run length fasta file:", runlength_reads_path)

    with open(runlength_reference_path, "w") as file:
        for contig_name in runlength_ref_sequences:
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_ref_sequences[contig_name][SEQUENCE] + "\n")

    with open(runlength_reads_path, "w") as file:
        for contig_name in runlength_reads_sequences:
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_reads_sequences[contig_name].sequence + "\n")

    output_sam_file_path, output_bam_file_path = align_minimap(output_dir=output_dir,
                                                               ref_sequence_path=runlength_reference_path,
                                                               reads_sequence_path=runlength_reads_path,
                                                               k=18)

    return output_bam_file_path


def main():
    # ref_fasta_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_14_18_19_4_215702_ref.fasta"
    # runlength_path = "/home/ryan/code/runnie_parser/data/synthetic_runnie_test_2019_3_14_18_19_4_215702_runnie.out"

    ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    shasta_path = "/home/ryan/code/runlength_analysis/data/ecoli_shasta_pileup_data.csv"

    output_parent_dir = "output/"
    output_dir = "runlength_matrix_from_runnie_output_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    reads_fasta_filename_prefix = ".".join(os.path.basename(shasta_path).split(".")[:-1])
    runlength_reads_fasta_filename = reads_fasta_filename_prefix + "_rle.fasta"
    runlength_reads_fasta_path = os.path.join(output_dir, runlength_reads_fasta_filename)

    shasta_handler = ShastaRunlengthHandler(shasta_path)

    print("Parsing shasta file...")

    name, read = shasta_handler.get_pileup_data(cutoff=sys.maxsize)

    read_data = {name:read}

    print("\nRLE encoding reference sequence...")

    runlength_ref_sequences = runlength_encode_fasta(fasta_sequence_path=ref_fasta_path)

    reads_vs_ref_bam_path = align_as_RLE(runlength_reference_path=runlength_ref_fasta_path,
                                         runlength_ref_sequences=runlength_ref_sequences,
                                         runlength_reads_path=runlength_reads_fasta_path,
                                         runlength_reads_sequences=read_data,
                                         output_dir=output_dir)

    print("\nGenerating matrices...")

    chromosomal_matrices = generate_runlength_frequency_matrix(runlength_ref_sequence_path=runlength_ref_fasta_path,
                                                               reads_vs_ref_bam_path=reads_vs_ref_bam_path,
                                                               runlength_ref_sequences=runlength_ref_sequences,
                                                               runlength_read_data=read_data)

    for matrix in chromosomal_matrices:
        save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                              frequency_matrices=matrix,
                                                              log_normalize=False,
                                                              plot=False,
                                                              default_type=float)

        save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                              frequency_matrices=matrix,
                                                              log_normalize=True,
                                                              pseudocount=15,
                                                              diagonal_bias=0,
                                                              plot=False,
                                                              default_type=float)

        plot_directional_residuals(matrix)

        # zero_mask = (matrix == 0)
        # nonzero_mask = numpy.invert(zero_mask)
        # matrix[zero_mask] += numpy.min(matrix[nonzero_mask])

        plot_base_matrices(sum_complementary_matrices(matrix), test_spot=False, normalize_matrices=False)
        plot_base_matrices(sum_complementary_matrices(matrix), test_spot=False, normalize_matrices=True)

        frequency_matrices = sum_complementary_matrices(matrix)

        save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                                 frequency_matrices=frequency_matrices,
                                                                 log_normalize=False,
                                                                 plot=False)

        save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                                 frequency_matrices=frequency_matrices,
                                                                 log_normalize=True,
                                                                 plot=False)


if __name__ == "__main__":
    main()
