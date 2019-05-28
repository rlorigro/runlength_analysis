from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from collections import defaultdict, Counter
from matplotlib import pyplot
import numpy
import math
import sys
import os


def count_runlength_per_character(sequence):
    """
    For each observed character, append observed runlengths of that character to a list in a dictionary with key=char
    :param sequence:
    :return:
    """
    character_counts = defaultdict(list)
    current_character = None

    for character in sequence:
        if character != current_character:
            character_counts[character].append(1)
        else:
            character_counts[character][-1] += 1

        current_character = character

    return character_counts


def print_all_counts(all_counts):
    for character in all_counts:
        print(">Character %s" % character)

        for item in sorted(all_counts[character].items(), key=lambda x: x[0]):
            print("%d %d" % (item[0], item[1]))


def print_all_counts_as_shasta_matrix(all_counts, max_count=50, pseudocount=1):
    """
    Format output as a row with index 0-50 (max_count) with normalized log probabilities of each length. Each comma
    separated row is preceded with a header describing which bases the prior corresponds to. Since observed reads are
    not directional, counts for complementary bases are summed together (A+T, G+C)
    :param all_counts:
    :param max_count:
    :return:
    """
    a_t_counts = all_counts["A"] + all_counts["T"]
    g_c_counts = all_counts["G"] + all_counts["C"]

    total = 0
    for i in range(max_count + 1):
        total += max(pseudocount, a_t_counts[i])

    line = list()
    for i in range(max_count + 1):
        count = max(pseudocount, a_t_counts[i])
        line.append("%.9f" % math.log((count/total),10))

    print(">AT prior")
    print(",".join(line))
    print()

    total = 0
    for i in range(max_count + 1):
        total += max(pseudocount, g_c_counts[i])

    line = list()
    for i in range(max_count + 1):
        count = max(pseudocount, g_c_counts[i])
        line.append("%.9f" % math.log((count/total),10))

    print(">GC prior")
    print(",".join(line))
    print()


def main():
    output_dir = "output/ref_run_lengths/bac/"
    filename_prefix = "ref_runlength_distribution"

    # ---- GIAB E. Coli - (dev machine) -------------------------
    # reference_file_path = "/home/ryan/data/GIAB/GRCh38_WG.fa"
    # reference_file_path = "/home/ryan/data/Nanopore/ecoli/refEcoli.fasta"
    reference_file_path = "/home/ryan/data/Nanopore/Human/paolo/LC2019/HG00733_BAC_VMRC62_79_clones.fasta"
    # -------------------------------------------------------------------------

    FileManager.ensure_directory_exists(output_dir)

    fasta_handler = FastaHandler(reference_file_path)
    contig_names = fasta_handler.get_contig_names()

    print(contig_names)
    print(sorted([(x,fasta_handler.get_chr_sequence_length(x)) for x in contig_names],key=lambda x: x[1]))

    all_counts = defaultdict(lambda: Counter())

    sys.stderr.write("reading fasta file...\n")
    sys.stderr.flush()

    c = 0
    for chromosome_name in contig_names:
        # if len(contig_names) > 1:
            # if not chromosome_name.startswith("chr"):
            #     print("WARNING: SKIPPING CHROMOSOME %s" % chromosome_name)
            #     continue

        # if c > 1:
        #     break
        c += 1

        sys.stderr.write("Parsing chromosome %s\n" % chromosome_name)
        sys.stderr.flush()

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reference_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name, start=0, stop=chromosome_length)
        character_counts = count_runlength_per_character(reference_sequence)

        # figure, axes = pyplot.subplots(nrows=len(character_counts.keys()), sharex=True)
        # figure.set_size_inches(6,12)

        for k,key in enumerate(character_counts.keys()):
            counts = character_counts[key]
            counter = Counter(counts)
            all_counts[key] += counter

            max_count = 50

            # step = 1
            # bins = numpy.arange(1, max_count + step, step=step)
            # frequencies, bins = numpy.histogram(counts, bins=bins, normed=False)
            #
            # frequencies = numpy.log10(frequencies+0.01) + 0.01
            #
            # center = (bins[:-1] + bins[1:])/2 - step/2 - 1
            #
            # axes[k].bar(center, frequencies, width=step, align="center")
            # axes[k].set_ylabel(str(key))
            # axes[k].set_ylim([-0.5,10])

        # axes[0].set_title(chromosome_name)

        # filename = filename_prefix + "_" + chromosome_name + ".png"
        # file_path = os.path.join(output_dir, filename)
        # figure.savefig(file_path)
        # pyplot.show()
        # pyplot.close()

    print_all_counts_as_shasta_matrix(all_counts, max_count=50)

    print_all_counts(all_counts)


def test_runlength_counter():
    sequence = "ACCGGGTTTTAAAAACCCCCCGGGGGGGTTTTTTTTAAAAAAAAACCCCCCCCCC"

    character_counts = count_runlength_per_character(sequence)

    print(character_counts)


if __name__ == "__main__":
    # test_runlength_counter()
    main()
