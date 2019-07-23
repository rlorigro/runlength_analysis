import sys
from handlers.FastqReader import FastqReader
from handlers.FastaHandler import FastaHandler
from matplotlib import pyplot
import argparse
import numpy

'''
Iterate a fastq file and find the read lengths
'''


def main(sequences_path, cutoff):
    if sequences_path.endswith(".fastq"):
        reads = FastqReader().iterate_file(path=sequences_path)
    elif sequences_path.endswith(".fasta") or sequences_path.endswith(".fa"):
        reads = FastaHandler(sequences_path).iterate_file()
    else:
        exit("Improper file format: %s" % sequences_path)

    n_reads = 0

    with open("assemble_long_segments.sh", "w") as file:
        for i, item in enumerate(reads):
            n_reads += 1

            if sequences_path.endswith(".fastq"):
                header, sequence, quality = item
            elif sequences_path.endswith(".fasta") or sequences_path.endswith(".fa"):
                header, sequence = item
            else:
                exit("unrecognized file type")

            if len(sequence) > cutoff:
                name = header.split(" ")[0]
                print(name, len(sequence))
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
