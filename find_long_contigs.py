import sys
from handlers.FastqHandler import FastqHandler
from handlers.FastaHandler import FastaHandler
from matplotlib import pyplot
import argparse
import numpy

'''
Iterate a fastq file and find the read lengths
'''


def main(sequences_path, cutoff):
    fasta = FastaHandler(sequences_path)
    names = fasta.get_contig_names()

    n_reads = 0

    with open("assemble_long_segments.sh", "w") as file:
        for i, name in enumerate(names):
            length = fasta.get_chr_sequence_length(name)

            n_reads += 1

            if length > cutoff:
                print(name, length)
                file.write("../build/shasta-install/bin/AssembleSegment.py " + name + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences",
        type=str,
        required=True,
        help="Fasta/fastq file path of sequences"
    )
    parser.add_argument(
        "--length",
        type=int,
        required=True,
        help="Length threshold"
    )

    args = parser.parse_args()

    main(sequences_path=args.sequences, cutoff=args.length)
