from handlers.FastaHandler import FastaHandler
from collections import defaultdict, Counter
import sys


def count_runlength_per_character(sequence, threshold, chromosome_name):
    """
    For each observed character, append observed runlengths of that character to a list in a dictionary with key=char
    :param sequence:
    :return:
    """
    character_counts = defaultdict(list)
    current_character = None
    
    for c,character in enumerate(sequence):
        if character != current_character:
            if character in {"A","C","G","T"}:
                if len(character_counts[character]) > 0 and character_counts[character][-1] > threshold:
                    print('\t'.join([chromosome_name, str(c), character, str(character_counts[character][-1])]))

            character_counts[character].append(1)
        else:
            character_counts[character][-1] += 1

        current_character = character

    return character_counts


def main():
    output_dir = "output/ref_run_lengths/"
    filename_prefix = "ref_runlength_distribution"

    reference_file_path = "/home/ryan/data/Nanopore/Human/paolo/LC2019/kishwar/shasta_assembly_GM24385_chr20.fasta"

    # ---- GIAB E. Coli - (dev machine) -------------------------
    # reference_file_path = "/home/ryan/data/GIAB/GRCh38_WG.fa"
    # reference_file_path = "/home/ryan/data/Nanopore/ecoli/refEcoli.fasta"
    # -------------------------------------------------------------------------

    threshold = 5

    fasta_handler = FastaHandler(reference_file_path)
    contig_names = fasta_handler.get_contig_names()

    all_counts = defaultdict(lambda: Counter())

    sys.stderr.write("reading fasta file...\n")
    sys.stderr.flush()

    c = 0
    for chromosome_name in contig_names:
        if len(contig_names) > 1:
            if not chromosome_name != "chr1":
                continue
        c += 1

        # sys.stderr.write("Parsing chromosome %s\n" % chromosome_name)
        # sys.stderr.flush()

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reference_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name,
                                                        stop=chromosome_length,
                                                        start=0)

        character_counts = count_runlength_per_character(sequence=reference_sequence,
                                                         threshold=threshold,
                                                         chromosome_name=chromosome_name)


if __name__ == "__main__":
    # test_runlength_counter()
    main()
