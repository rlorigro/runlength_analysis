from handlers.FastqHandler import FastqHandler
import sys
import os


def main():
    fastq_path = "/home/ryan/code/runlength_analysis/data/runnie_subset_test_flipflop_regional_0to10k.fastq"
    output_directory = "output/"
    output_filename_prefix = ".".join(os.path.basename(fastq_path).split(".")[:-1])
    output_filename = output_filename_prefix + ".fasta"

    output_path = os.path.join(output_directory, output_filename)

    print("WRITING:", output_path, "...")

    handler = FastqHandler(fastq_path)

    with open(output_path, "w") as file:
        for r,read in enumerate(handler.iterate_file()):
            sys.stderr.write("\r %d" % r)
            file.write(">"+read.id + "\n")
            file.write(read.sequence + "\n")

    sys.stderr.write("\nDone\n")


if __name__ == "__main__":
    main()