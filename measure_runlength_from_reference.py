from handlers.FastaHandler import FastaHandler
from collections import defaultdict
from matplotlib import pyplot
import numpy


def count_runlength_per_character(sequence):
    character_counts = defaultdict(list)
    current_character = None

    for character in sequence:
        if character != current_character:
            character_counts[character].append(1)
        else:
            character_counts[character][-1] += 1

        current_character = character

    return character_counts


def main():
    output_dir = "output/ref_run_lengths/"
    filename_prefix = "ref_runlength_distribution"

    # ---- GIAB E. Coli - (dev machine) -------------------------
    reference_file_path = "/home/ryan/data/GIAB/GRCh38_WG.fa"

    # -------------------------------------------------------------------------

    fasta_handler = FastaHandler(reference_file_path)
    contig_names = fasta_handler.get_contig_names()

    print(contig_names)
    print(sorted([(x,fasta_handler.get_chr_sequence_length(x)) for x in contig_names],key=lambda x: x[1]))

    for chromosome_name in contig_names:
        if not chromosome_name.startswith("chr"):
            continue

        chromosome_length = fasta_handler.get_chr_sequence_length(chromosome_name)

        reference_sequence = fasta_handler.get_sequence(chromosome_name=chromosome_name, start=0, stop=chromosome_length)

        # print(reference_sequence[380021-12:380021+2])

        character_counts = count_runlength_per_character(reference_sequence)

        figure, axes = pyplot.subplots(nrows=len(character_counts.keys()), sharex=True)
        figure.set_size_inches(4,12)

        for k,key in enumerate(character_counts.keys()):
            print("Reading %s" % chromosome_name)

            counts = character_counts[key]
            max_count = 50

            step = 1
            bins = numpy.arange(1, max_count + step, step=step)
            frequencies, bins = numpy.histogram(counts, bins=bins, normed=False)

            print(bins)
            print(frequencies)

            frequencies = numpy.log10(frequencies+0.01) + 0.01

            print(bins.shape)
            center = (bins[:-1] + bins[1:])/2 - step/2 - 1

            axes[k].bar(center, frequencies, width=step, align="center")
            axes[k].set_ylabel(str(key))
            axes[k].set_ylim([-0.5,10])

        filename = filename_prefix + "_" + chromosome_name + ".png"
        figure.savefig(filename)
        figure.title()
        pyplot.show()
        pyplot.close()


def test_runlength_counter():
    sequence = "ACCGGGTTTTAAAAACCCCCCGGGGGGGTTTTTTTTAAAAAAAAACCCCCCCCCC"

    character_counts = count_runlength_per_character(sequence)

    print(character_counts)


if __name__ == "__main__":
    # test_runlength_counter()
    main()
