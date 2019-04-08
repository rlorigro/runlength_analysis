from handlers.FileManager import FileManager
from collections import defaultdict
from copy import copy
import random
import csv
import os


def complement_sequence(sequence, reverse=True):
    """
    Complement a sequence of bases in list form (NOT string)
    :param sequence: list of strings (bases)
    :return:
    """
    sequence_complement = list()

    if reverse:
        for base in reversed(sequence):
            sequence_complement.append(complement_base(base))
    else:
        for base in sequence:
            sequence_complement.append(complement_base(base))

    return sequence_complement


def complement_base(base):
    """
    Complement a base
    :param base:
    :return:
    """
    if base == "A":
        complement = "T"

    elif base == "T":
        complement = "A"

    elif base == "C":
        complement = "G"

    elif base == "G":
        complement = "C"

    else:
        exit("ERROR: invalid base encountered in complement fn:", base)

    return complement


def main(output_dir="data/"):
    filename_prefix = "synthetic_runlength_test_" + FileManager.get_datetime_string()
    runlength_reference_path = os.path.join(output_dir, filename_prefix + "_ref.fasta")
    runlength_reads_path = os.path.join(output_dir, filename_prefix + "_reads.fasta")

    reverse_complement = True

    n_repeats = 12
    coverage = 12

    ref_max_runlength = 8

    base_pool = ["A", "T", "G", "C"]

    base_length_offsets = {"A":0, "T":1, "G":2, "C":3}

    ref_sequence = list()

    ref_lengths = list()
    ref_bases = list()

    read_output_lines = list()

    ref_sequence_name = "synthetic_ref_0"
    for i in range(n_repeats):
        ref_runlengths = {b: list(range(1, ref_max_runlength + 1)) for b in base_pool}

        for i in range(ref_max_runlength):
            bases = copy(base_pool)
            random.shuffle(bases)

            if len(ref_bases) > 0:
                while bases[0] == ref_bases[-1]:
                    random.shuffle(bases)

            for base in bases:
                lengths = ref_runlengths[base]
                length = lengths.pop()

                ref_runlengths[base] = lengths

                ref_sequence.extend([base]*length)
                ref_lengths.append(length)
                ref_bases.append(base)

    ref_sequence = "".join(ref_sequence)

    for c in range(coverage):
        read_output_lines.append(">synthetic_read_%d"%c)
        sequence = list()
        for i in range(len(ref_lengths)):
            base = ref_bases[i]
            runlength = ref_lengths[i] + base_length_offsets[base]

            sequence.extend([base]*runlength)

        sequence = "".join(sequence)
        read_output_lines.append(sequence)

        if reverse_complement:
            read_output_lines.append(">synthetic_read_reverse_%d" % c)
            sequence = complement_sequence(sequence=sequence, reverse=True)
            sequence = "".join(sequence)
            read_output_lines.append(sequence)

    print("saving file:", runlength_reference_path)
    with open(runlength_reference_path, "w") as file:
        file.write(">"+ref_sequence_name+"\n")
        file.write(ref_sequence + "\n")

    print("saving file:", runlength_reads_path)
    with open(runlength_reads_path, "w") as file:
        for line in read_output_lines:
            file.write(line + "\n")


if __name__ == "__main__":
    main()
