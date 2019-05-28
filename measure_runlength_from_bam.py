from handlers.BamHandler import BamHandler
from collections import defaultdict, Counter
import argparse
import sys
import os


def count_runlength_per_character(sequence):
    """
    For each observed character, append observed runlengths of that character to a list in a dictionary with key=char
    :param sequence:
    :return:
    """
    character_count_distributions = defaultdict(lambda: defaultdict(int))
    current_length = 0
    current_character = None

    for character in sequence:
        if character != current_character and current_character is not None:
            character_count_distributions[current_character][current_length] += 1
            current_length = 1
        else:
            current_length += 1

        current_character = character

    character_count_distributions[current_character][current_length] += 1

    return character_count_distributions


def print_all_counts(all_counts):
    for character in all_counts:
        print(">Character %s" % character)

        for item in sorted(all_counts[character].items(), key=lambda x: x[0]):
            print("%d %d" % (item[0], item[1]))


def main(bam_file_path, cutoff, contig_name):
    # ---- GIAB E. Coli - (dev machine) ---------------------------------------
    # bam_file_path = "/home/ryan/data/GIAB/GRCh38_WG.fa"
    # bam_file_path = "/home/ryan/data/Nanopore/ecoli/flapppie/03_22_19_R941_gEcoli_first_410k_VS_refEcoli.sorted.bam"
    # -------------------------------------------------------------------------

    bam_handler = BamHandler(bam_file_path)
    reads = bam_handler.get_reads(chromosome_name=contig_name, start=None, stop=None)

    all_counts = defaultdict(lambda: Counter())

    sys.stderr.write("reading file...\n")
    sys.stderr.flush()

    c = 0
    for read in reads:
        if read.mapping_quality <= 5 or read.is_secondary or read.is_unmapped \
                or read.is_qcfail:
            continue

        c += 1

        if c % 100 == 0:
            sys.stderr.write("\rParsed %d reads" % c)

        if c > cutoff:
            break

        sequence = read.query_sequence

        # print(read.query_name)
        # print(len(sequence))
        # print(sequence[:10])

        character_counts = count_runlength_per_character(sequence)

        for character in character_counts:
            all_counts[character] += character_counts[character]

    sys.stderr.write("\n")

    for character in sorted(all_counts):
        print(">%s" % character)
        for length in sorted(all_counts[character].keys()):
            print(length, all_counts[character][length])


def test_runlength_counter():
    sequence = "ACCGGGTTTTAAAAACCCCCCGGGGGGGTTTTTTTTAAAAAAAAACCCCCCCCCC"

    character_counts = count_runlength_per_character(sequence)

    for character in character_counts:
        print(character)

        print(character_counts[character].keys())
        for length in sorted(character_counts[character].keys()):
            print(length, character_counts[character][length])


if __name__ == "__main__":
    # test_runlength_counter()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="path to directory containing marginpolish TSV files"
    )
    parser.add_argument(
        "--contig",
        type=str,
        required=True,
        help="name of contig from which to fetch reads from BAM file"
    )
    parser.add_argument(
        "--cutoff",
        type=int,
        default=sys.maxsize,
        required=False,
        help="limit on how many reads to parse (default = %d)" % sys.maxsize
    )

    args = parser.parse_args()
    main(bam_file_path=args.bam, cutoff=args.cutoff, contig_name=args.contig)
