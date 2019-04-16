from handlers.MarginpolishRunlengthHandler import MarginpolishRunlengthHandler
from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from multiprocessing import Manager, Pool, cpu_count
from modules.align import align_minimap
from modules.matrix import *
import argparse
import numpy
import copy
import sys
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

DEBUG_INDEX = 0


def update_frequency_matrix(observed_pileup, read_consensus_base, true_base, true_length, alignment_reversal, matrix):
    # Did alignment reverse complement the sequence (via BAM) and ref (via pysam)? if so, revert to forward direction
    # if alignment_reversal:
    #     true_base = complement_base(true_base)

    true_base_index = BASE_TO_INDEX[true_base]

    # global DEBUG_INDEX
    # DEBUG_INDEX += 1
    #
    # if DEBUG_INDEX < 20:
    #     # if observed_length == 1 and true_length > 2:
    #     bases = [p.base for p in observed_pileup]
    #     lengths = [p.length for p in observed_pileup]
    #     counts = [p.count for p in observed_pileup]
    #     reversals = [p.reversal for p in observed_pileup]
    #     observations = list(zip(lengths, counts, bases, reversals))
    #     observations = sorted(observations, key=lambda x: x[0])
    #     print(true_base, read_consensus_base, true_length, read_consensus_length, [[l] * c for l, c, b, r in observations])
    #     print(true_base, read_consensus_base, true_length, read_consensus_length, [[b] * c for l, c, b, r in observations])
    #     # print(true_base, read_consensus_base, true_length, read_consensus_length, [[int(r)] * c for l, c, b, r in observations])

    for item in observed_pileup:
        observed_length = item.length
        observed_base = item.base
        weight = item.weight

        if observed_length > MAX_RUNLENGTH:
            observed_length = MAX_RUNLENGTH

        if true_length > MAX_RUNLENGTH:
            true_length = MAX_RUNLENGTH

        if observed_base == read_consensus_base:
            matrix[int(item.reversal), true_base_index, true_length, observed_length] += weight

        # matrix[int(item.reversal), base_index, true_length, observed_length] += count

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


def reverse_runlength_read(runlength_read):
    sequence = runlength_read.sequence
    pileup = runlength_read.pileup

    reversed_sequence = list()
    reversed_pileup = list()

    for i in reversed(range(len(sequence))):
        reversed_sequence.append(sequence[i])
        reversed_pileup.append(pileup[i])

    reversed_sequence = "".join(reversed_sequence)

    runlength_read.sequence = reversed_sequence
    runlength_read.pileup = reversed_pileup

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
        # read_consensus_length = read_consensus_lengths_segment[i-alignment_position]

        ref_runlength = ref_runlengths[i-alignment_position]

        observed_pileup = observed_pileup_segment[i-alignment_position]

        update_frequency_matrix(observed_pileup=observed_pileup,
                                read_consensus_base=read_base,
                                # read_consensus_length=read_consensus_length,
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
                    # read_consensus_lengths_segment=read_consensus_lengths_segment,
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
    r = 0
    for read in reads:
        if read.mapping_quality > 0 and not read.is_secondary:
            r += 1

            read_id = read.query_name
            runlength_read = runlength_read_data[read_id]

            if read.is_reverse:
                runlength_read = reverse_runlength_read(copy.copy(runlength_read))

            ref_alignment_start = read.reference_start
            ref_alignment_stop = get_read_stop_position(read)
            ref_length = ref_alignment_stop - ref_alignment_start
            cigar_tuples = read.cigartuples
            read_sequence = runlength_read.sequence
            # read_sequence = read.query_sequence

            # read_consensus_lengths = runlength_read.lengths
            observed_pileup = runlength_read.pileup

            ref_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                      start=ref_alignment_start,
                                                      stop=ref_alignment_stop + 10)

            ref_runlengths = complete_ref_runlengths[ref_alignment_start:ref_alignment_stop + 10]

            # print("\nREAD BEFORE CLIPPING")
            # print("name:\t", read_id)
            # print("is_reversed:\t", read.is_reverse)
            # print("read_consensus_lengths:\t\t", len(read_consensus_lengths))
            # print("read_sequence:\t\t\t\t", len(read_sequence))
            # print("observed_pileup:\t\t\t", len(observed_pileup))
            # print("ref_sequence:\t\t\t\t", len(ref_sequence))
            # print("ref_runlengths:\t\t\t\t", len(ref_runlengths))
            # print("first 10 cigars:\t\t", cigar_tuples[:10])
            # print("last 10 cigars:\t\t", cigar_tuples[-10:])
            #
            # print(ref_runlengths[0:35])
            # print(read_consensus_lengths[0:35])
            # print(ref_sequence[0:35])
            # print(read_sequence[0:35])
            # print()

            n_initial_hard_clipped_bases = 0
            if cigar_tuples[0][0] == 5:
                n_initial_hard_clipped_bases = cigar_tuples[0][1]

            n_final_hard_clipped_bases = 0
            if cigar_tuples[-1][0] == 5:
                n_final_hard_clipped_bases = cigar_tuples[-1][1]

            clipped_start = n_initial_hard_clipped_bases
            clipped_stop = len(observed_pileup) - n_final_hard_clipped_bases

            read_sequence = read_sequence[clipped_start:clipped_stop]
            # read_consensus_lengths = read_consensus_lengths[clipped_start:clipped_stop]
            observed_pileup = observed_pileup[clipped_start:clipped_stop]

            # print("\nREAD AFTER CLIPPING")
            # print("name:\t", read_id)
            # print("is_reversed:\t", read.is_reverse)
            # print("read_consensus_lengths:\t\t", len(read_consensus_lengths))
            # print("read_sequence:\t\t\t\t", len(read_sequence))
            # print("observed_pileup:\t\t\t", len(observed_pileup))
            # print("ref_sequence:\t\t\t\t", len(ref_sequence))
            # print("ref_runlengths:\t\t\t\t", len(ref_runlengths))
            # print("first 10 cigars:\t\t", cigar_tuples[:10])
            # print("last 10 cigars:\t\t", cigar_tuples[-10:])
            #
            # print(ref_runlengths[0:35])
            # print(read_consensus_lengths[0:35])
            # print(ref_sequence[0:35])
            # print(read_sequence[0:35])
            # print()

            # read_index: index of read sequence
            # ref_index: index of reference sequence
            read_index = 0
            ref_index = 0
            found_valid_cigar = False

            for c,cigar in enumerate(cigar_tuples):
                cigar_code = cigar[0]
                length = cigar[1]

                # skip parsing the first segment if it is not a match
                if cigar_code != 0 and found_valid_cigar is False:
                    # only increment the read index if the non-match cigar code is INS or SOFTCLIP
                    if cigar_code == 1 or cigar_code == 4:
                        read_index += length
                    continue

                found_valid_cigar = True

                # get the sequence segments that are affected by this operation
                read_sequence_segment = read_sequence[read_index:read_index+length]
                ref_sequence_segment = ref_sequence[ref_index:ref_index+length]

                ref_runlengths_segment = ref_runlengths[ref_index:ref_index+length]

                # read_consensus_lengths_segment = read_consensus_lengths[read_index:read_index+length]
                observed_pileup_segment = observed_pileup[read_index:read_index+length]

                ref_index_increment, read_index_increment = parse_cigar_tuple(cigar_code=cigar_code,
                                                                              length=length,
                                                                              alignment_position=ref_alignment_start + ref_index,
                                                                              read_sequence_segment=read_sequence_segment,
                                                                              # read_consensus_lengths_segment=read_consensus_lengths_segment,
                                                                              ref_sequence_segment=ref_sequence_segment,
                                                                              ref_runlengths_segment=ref_runlengths_segment,
                                                                              observed_pileup_segment=observed_pileup_segment,
                                                                              reversal_status=read.is_reverse,
                                                                              matrix=matrix)

                # increase the read index iterator
                read_index += read_index_increment
                ref_index += ref_index_increment

    return r


def generate_runlength_frequency_matrix(runlength_ref_sequence_path, reads_vs_ref_bam_path,
                                        runlength_ref_sequences, runlength_read_data, max_runlength):
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
        shape = [2,4,max_runlength+1,max_runlength+1]
        matrix = numpy.zeros(shape, dtype=numpy.float64)

        try:
            bam_handler = BamHandler(reads_vs_ref_bam_path)
            fasta_handler = FastaHandler(runlength_ref_sequence_path)

            chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

            reads = bam_handler.get_reads(chromosome_name=chromosome_name, start=0, stop=chromosome_length)

            n_reads = parse_reads(chromosome_name=chromosome_name,
                                  fasta_handler=fasta_handler,
                                  reads=reads,
                                  complete_ref_runlengths=runlength_ref_sequences[chromosome_name][LENGTHS],
                                  runlength_read_data=runlength_read_data,
                                  matrix=matrix)

        except Exception as e:
            print("ERROR FILE NOT FOUND\n %s \n %s" % (runlength_ref_sequence_path, reads_vs_ref_bam_path))
            print(e)
            continue

        if n_reads > 0:
            yield (chromosome_name, matrix)
        else:
            sys.stderr.write("No reads found for chromosome: %s\n" % chromosome_name)


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
        sequence = fasta_handler.get_sequence(chromosome_name=contig_name, start=None, stop=None)

        bases, lengths = runlength_encode(sequence)

        runlength_sequences[contig_name] = (bases, lengths)

        sys.stderr.write("\rRun length encoded %s            " % contig_name)

    sys.stderr.write("\n")

    return runlength_sequences


def align_as_RLE(runlength_reference_path, runlength_ref_sequences, runlength_reads_path, runlength_reads_sequences, output_dir, max_threads):
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
                                                               preset="asm20",
                                                               k=18,
                                                               max_threads=max_threads)

    return output_bam_file_path


def parse_coverage_data(output_dir, marginpolish_path, ref_fasta_path, runlength_ref_fasta_path):  #, all_matrices, counter):
    reads_fasta_filename_prefix = ".".join(os.path.basename(marginpolish_path).split(".")[:-1])
    runlength_reads_fasta_filename = reads_fasta_filename_prefix + "_rle.fasta"
    runlength_reads_fasta_path = os.path.join(output_dir, runlength_reads_fasta_filename)

    shasta_handler = MarginpolishRunlengthHandler(marginpolish_path)

    # print("Parsing tsv file...")

    name, read = shasta_handler.get_pileup_data(cutoff=sys.maxsize)

    read_data = {name: read}

    # print("\nRLE encoding reference sequence...")

    runlength_ref_sequences = runlength_encode_fasta(fasta_sequence_path=ref_fasta_path)

    reads_vs_ref_bam_path = align_as_RLE(runlength_reference_path=runlength_ref_fasta_path,
                                         runlength_ref_sequences=runlength_ref_sequences,
                                         runlength_reads_path=runlength_reads_fasta_path,
                                         runlength_reads_sequences=read_data,
                                         output_dir=output_dir,
                                         max_threads=1)

    # print("\nGenerating matrices...")

    chromosomal_matrices = generate_runlength_frequency_matrix(runlength_ref_sequence_path=runlength_ref_fasta_path,
                                                               reads_vs_ref_bam_path=reads_vs_ref_bam_path,
                                                               runlength_ref_sequences=runlength_ref_sequences,
                                                               runlength_read_data=read_data,
                                                               max_runlength=MAX_RUNLENGTH)

    numpy_filename_prefix = ".".join(os.path.basename(reads_vs_ref_bam_path).split(".")[:-1])

    for chromosome_name, matrix in chromosomal_matrices:
        numpy_filename = numpy_filename_prefix + "_" + chromosome_name + "_matrix.pkl"
        numpy_file_path = os.path.join(output_dir, numpy_filename)
        matrix.dump(numpy_file_path)

        #     all_matrices.append(matrix)

    # counter.value += 1
    # sys.stderr.write("\r%d" % counter.value)


def main(ref_fasta_path, marginpolish_parent_dir, max_threads):
    # ref_fasta_path = "/home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta"
    # marginpolish_parent_dir = "/home/ryan/data/Nanopore/ecoli/benedict/flappy/ecoli.repeatCount.flippity_flappity"

    marginpolish_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=marginpolish_parent_dir,
                                                                file_extension=".tsv")

    # Test
    # ref_fasta_path = "/home/ryan/code/runlength_analysis/data/synthetic_coverage_data_marginpolish_2019_4_12_13_31_6_689541_ref.fasta"
    # marginpolish_paths = ["/home/ryan/code/runlength_analysis/data/synthetic_coverage_data_marginpolish_2019_4_12_13_31_6_689541.tsv"]

    if max_threads is None:
        max_threads = max(1, cpu_count() - 2)

    output_parent_dir = "output/"
    output_dir = "runlength_matrix_from_marginpolish_output_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    # manager = Manager()
    # all_matrices = manager.list()
    # counter = manager.Value('i', 0)

    process_arguments = list()
    for marginpolish_path in marginpolish_paths:
        process_arguments.append([output_dir, marginpolish_path, ref_fasta_path, runlength_ref_fasta_path]) # , all_matrices, counter])

    if max_threads > len(process_arguments):
        max_threads = len(process_arguments)

    with Pool(processes=max_threads, maxtasksperchild=5) as pool:
        pool.starmap(parse_coverage_data, process_arguments, chunksize=1)

    matrix_file_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=output_dir,
                                                               file_extension=".pkl")

    all_matrices = list()
    for path in matrix_file_paths:
        all_matrices.append(numpy.load(path))

    matrix = numpy.stack(all_matrices, axis=0)
    matrix = numpy.sum(matrix, axis=0)

    save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                          frequency_matrices=matrix,
                                                          chromosome_name="genomic",
                                                          log_normalize=False,
                                                          diagonal_bias=0,
                                                          plot=False,
                                                          default_type=float)

    save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                          frequency_matrices=matrix,
                                                          chromosome_name="genomic",
                                                          log_normalize=True,
                                                          pseudocount=15,
                                                          diagonal_bias=0,
                                                          plot=False,
                                                          default_type=float)

    nondirectional_matrix = sum_complementary_matrices(matrix, max_runlength=MAX_RUNLENGTH)

    save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                             frequency_matrices=nondirectional_matrix,
                                                             chromosome_name="genomic",
                                                             log_normalize=False,
                                                             plot=False)

    save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                             frequency_matrices=nondirectional_matrix,
                                                             chromosome_name="genomic",
                                                             log_normalize=True,
                                                             plot=False)

    # zero_mask = (matrix == 0)
    # nonzero_mask = numpy.invert(zero_mask)
    # matrix[zero_mask] += numpy.min(matrix[nonzero_mask])

    plot_directional_complementary_diff(matrix, max_runlength=MAX_RUNLENGTH)
    plot_base_matrices(matrix, test_spot=False, normalize_matrices=False)
    plot_base_matrices(matrix, test_spot=False, normalize_matrices=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_dir", "-i",
        type=str,
        required=True,
        help="path to directory containing marginpolish TSV files"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="file path of FASTA reference sequence file"
    )
    parser.add_argument(
        "--max_threads", "-t",
        type=int,
        required=False,
        help="total number of vCPU to allocate to this job"
    )
    args = parser.parse_args()

    main(ref_fasta_path=args.ref, marginpolish_parent_dir=args.input_dir, max_threads=args.max_threads)


"""
python3 measure_runlength_distribution_from_marginpolish.py -i /home/ryan/data/Nanopore/ecoli/benedict/flappy/ecoli.repeatCount.flippity_flappity --ref /home/ryan/data/Nanopore/ecoli/miten/refEcoli.fasta -t 12
"""
