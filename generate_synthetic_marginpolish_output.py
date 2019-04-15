from handlers.FileManager import FileManager
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


def generate_observations_reverse_complement(read_runlengths, read_bases, n_coverage, observations, scale_coverage=False):
    index = 1

    # Generate observations that correspond to reference
    for r in reversed(range(len(read_runlengths))):
        for base in reversed(read_bases):
            observations_at_current_base = [str(index), complement_base(base)]

            for i in range(n_coverage):
                if random.random() > 0.5:
                    strand = "+"
                    # base_string = base
                else:
                    strand = "-"
                    # base_string = complement_base(base)

                runlength = read_runlengths[r]
                data_string = complement_base(base) + strand + str(runlength)

                weight = 0.3
                weight_string = str(weight)

                observation = ",".join([data_string,weight_string])
                observations_at_current_base.append(observation)

            observations.append(observations_at_current_base)

            index += 1

        if scale_coverage:
            n_coverage -= 1


def generate_observations(read_runlengths, read_bases, n_coverage, observations, scale_coverage=False):
    index = 1

    # Generate observations that correspond to reference
    for r in range(len(read_runlengths)):
        for base in read_bases:
            observations_at_current_base = [str(index), base]

            for i in range(n_coverage):
                if random.random() > 0.5:
                    strand = "+"
                    # base_string = base
                else:
                    strand = "-"
                    # base_string = complement_base(base)

                runlength = read_runlengths[r]
                data_string = base + strand + str(runlength)

                weight = 0.3
                weight_string = str(weight)

                observation = ",".join([data_string,weight_string])
                observations_at_current_base.append(observation)

            observations.append(observations_at_current_base)

            index += 1

        if n_coverage:
            n_coverage -= 1


def generate_sequences(ref_max_runlength, read_max_runlength, n_coverage, scale_coverage=False, reverse_complement=False):
    """
    Make a synthetic reference and a set of reads such that every run length (within a range) is found in the ref and in
    the reads. Additionally, 50% of the reads should be reverse complements.
    :param ref_max_runlength:
    :param read_max_runlength:
    :param n_sequences:
    :return:
    """
    ref_bases = ["A", "T", "G", "C"]
    read_bases = ["A", "T", "G", "C"]

    ref_runlengths = list(range(1,ref_max_runlength+1))
    read_runlengths = list(range(1,read_max_runlength+1))

    ref_sequence = list()
    observations = [["REF_INDEX","REF_BASE","REPEAT_COUNT_OBSxN(READ_BASE:READ_STRAND:REPEAT_COUNT,WEIGHT)"], ["0","N"]]

    # Generate non-RLE reference sequence with runlengths that correspond to the range specified
    for runlength in ref_runlengths:
        for base in ref_bases:
            ref_sequence.extend([base]*runlength)

    if scale_coverage:
        n_coverage = ref_max_runlength

    if reverse_complement:
        generate_observations_reverse_complement(read_runlengths=read_runlengths,
                                                 read_bases=read_bases,
                                                 n_coverage=n_coverage,
                                                 observations=observations,
                                                 scale_coverage=scale_coverage)

    else:
        generate_observations(read_runlengths=read_runlengths,
                              read_bases=read_bases,
                              n_coverage=n_coverage,
                              observations=observations,
                              scale_coverage=scale_coverage)

    ref_sequence = "".join(ref_sequence)

    return ref_sequence, observations


def main():
    """
    Make a synthetic reference and a set of reads and save them to fasta files as reads.fasta and ref.fasta
    :return:
    """
    output_dir = "data/"
    FileManager.ensure_directory_exists(output_dir)

    n_coverage = 2

    ref_max_runlength = 50
    read_max_runlength = 50

    ref_sequence, observations = generate_sequences(ref_max_runlength=ref_max_runlength,
                                                    read_max_runlength=read_max_runlength,
                                                    n_coverage=n_coverage,
                                                    scale_coverage=True)

    datetime_string = FileManager.get_datetime_string()
    filename = "synthetic_coverage_data_marginpolish_" + datetime_string + ".tsv"
    output_path = os.path.join(output_dir, filename)

    file = open(output_path, "w")
    writer = csv.writer(file, delimiter="\t")
    for line in observations:
        writer.writerow(line)
    file.close()

    filename = "synthetic_coverage_data_marginpolish_" + datetime_string + "_ref.fasta"
    output_path = os.path.join(output_dir, filename)

    with open(output_path, "w") as file:
        file.write(">ref_0\n")
        file.write(ref_sequence)


if __name__ == "__main__":
    main()
