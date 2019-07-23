from handlers.FastaHandler import FastaHandler
from handlers.FileManager import FileManager
from collections import defaultdict, Counter
from matplotlib import pyplot
import argparse
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


def print_all_counts(all_counts, output_dir):
    file_path = os.path.join(output_dir, "raw_counts_genomic.txt")

    with open(file_path, "w") as file:
        for character in all_counts:
            file.write(">Character %s" % character + "\n")

            for item in sorted(all_counts[character].items(), key=lambda x: x[0]):
                file.write("%d %d" % (item[0], item[1]) + "\n")


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


def plot_counts_as_histogram(axes, counts, max_count, step):
    bins = numpy.arange(1, max_count + step, step=step)
    frequencies, bins = numpy.histogram(counts, bins=bins, normed=False)

    frequencies = numpy.log10(frequencies + 0.01) + 0.01

    center = (bins[:-1] + bins[1:]) / 2 - step / 2 - 1

    axes.bar(center, frequencies, width=step, align="center")
    axes.set_ylabel("Log10 Frequency")
    axes.set_xlabel("Run Length (bp)")
    axes.set_ylim([-0.5, 10])


def main(reference_file_path):
    input_prefix_name = os.path.basename(reference_file_path).split("/")[-1].split(".")[0]
    output_dir = os.path.join("output/ref_run_lengths/", input_prefix_name)
    filename_prefix = "ref_runlength_distribution"

    FileManager.ensure_directory_exists(output_dir)

    fasta_handler = FastaHandler(reference_file_path)
    contig_names = fasta_handler.get_contig_names()

    print(contig_names)
    print(sorted([(x,fasta_handler.get_chr_sequence_length(x)) for x in contig_names],key=lambda x: x[1]))

    all_counts = defaultdict(lambda: Counter())
    raw_counts_AT = list()
    raw_counts_GC = list()

    sys.stderr.write("reading fasta file...\n")
    sys.stderr.flush()

    max_count = 100
    step = 1
    c = 0
    for chromosome_name in contig_names:
        # if len(contig_names) > 1:
        #     if not chromosome_name.startswith("chr") or "alt" in chromosome_name or "v" in chromosome_name:
        #         print("WARNING: SKIPPING CHROMOSOME %s" % chromosome_name)
        #         continue

        # if c == 1:
        #     break
        c += 1

        sys.stderr.write("Parsing chromosome %s\n" % chromosome_name)
        sys.stderr.flush()

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reference_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name, start=0, stop=chromosome_length)
        character_counts = count_runlength_per_character(reference_sequence)

        figure, axes = pyplot.subplots(nrows=len(character_counts.keys()), sharex=True)
        figure.set_size_inches(6,12)

        for k,key in enumerate(character_counts.keys()):
            counts = character_counts[key]
            counter = Counter(counts)
            all_counts[key] += counter

            if key in {"C","G"}:
                raw_counts_GC += counts

            if key in {"A","T"}:
                raw_counts_AT += counts

            plot_counts_as_histogram(axes=axes[k], counts=counts, max_count=max_count, step=step)

            axes[k].set_ylabel(str(key))
            axes[k].set_ylim([-0.5,10])

        axes[0].set_title(chromosome_name)

        filename = filename_prefix + "_" + chromosome_name + ".png"
        file_path = os.path.join(output_dir, filename)
        figure.savefig(file_path)
        # pyplot.show()
        pyplot.close()

    figure, axes = pyplot.subplots(nrows=2)

    filename = filename_prefix + "_genomic.png"
    file_path = os.path.join(output_dir, filename)

    plot_counts_as_histogram(axes=axes[0], counts=raw_counts_AT, max_count=max_count, step=step)
    plot_counts_as_histogram(axes=axes[1], counts=raw_counts_GC, max_count=max_count, step=step)
    axes[0].set_ylabel("AT Log10 Frequency")
    axes[1].set_ylabel("GC Log10 Frequency")

    figure.savefig(file_path)
    # pyplot.show()
    pyplot.close()

    print_all_counts_as_shasta_matrix(all_counts, max_count=50)
    print_all_counts(all_counts, output_dir)


def test_runlength_counter():
    sequence = "ACCGGGTTTTAAAAACCCCCCGGGGGGGTTTTTTTTAAAAAAAAACCCCCCCCCC"

    character_counts = count_runlength_per_character(sequence)

    print(character_counts)


if __name__ == "__main__":
    # test_runlength_counter()

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences", "-i",
        type=str,
        required=True,
        help="path to fastq file"
    )
    args = parser.parse_args()
    main(args.sequences)

