from handlers.FileManager import FileManager
from collections import defaultdict
from copy import copy
import random
import csv
import os


def read_weibull_params(weibull_params_path="data/weibull_parameters_by_mode.csv"):
    modal_parameters = defaultdict(list)
    with open(weibull_params_path, "r") as file:
        reader = csv.reader(file)

        for line in reader:
            mode = int(line[0])
            scale = float(line[1])
            shape = float(line[2])

            # print(mode, scale, shape)

            modal_parameters[mode].append([scale,shape])

    return modal_parameters


def main(output_dir="data/"):
    filename_prefix = "synthetic_runnie_test_" + FileManager.get_datetime_string()
    runlength_reference_path = os.path.join(output_dir, filename_prefix + "_ref.fasta")
    runnie_output_path = os.path.join(output_dir, filename_prefix + "_runnie.out")

    modal_parameters = read_weibull_params()

    n_repeats = 12
    coverage = 12

    ref_max_runlength = 8

    base_pool = ["A", "T", "G", "C"]

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
        read_output_lines.append("# synthetic_read_%d"%c)
        for i in range(len(ref_lengths)):
            runlength = ref_lengths[i]
            base = ref_bases[i]

            scale,shape = random.choice(modal_parameters[runlength])

            hex_scale = scale.hex()
            hex_shape = shape.hex()

            line = [base, hex_shape, hex_scale]
            line = list(map(str,line))
            line = "\t".join(line)
            read_output_lines.append(line)

            print(line)

    print(ref_sequence)

    print("saving file:", runlength_reference_path)
    with open(runlength_reference_path, "w") as file:
        file.write(">"+ref_sequence_name+"\n")
        file.write(ref_sequence + "\n")

    print("saving file:", runnie_output_path)
    with open(runnie_output_path, "w") as file:
        for line in read_output_lines:
            file.write(line + "\n")


if __name__ == "__main__":
    main()
