from handlers.RunlengthHandler_v2 import RunlengthHandler
from handlers.FastqHandler import FastqHandler
from handlers.FileManager import FileManager
import random
import sys
import os


def read_names_from_file(file_path):
    names = list()

    with open(file_path, "r") as file:
        for line in file:
            name = line.strip()
            names.append(name)

    return names


def extract_runnie_reads_by_name(runnie_path, output_dir, output_filename_suffix, names):
    output_filename = "runnie_subset_" + output_filename_suffix + ".out"
    output_path = os.path.join(output_dir, output_filename)
    FileManager.ensure_directory_exists(output_dir)

    runnie_handler = RunlengthHandler(runnie_path)

    runnie_handler.extract_reads_by_id(id_set=names, output_path=output_path, print_status=True)

    return output_path


def extract_fastq_reads_by_name(fastq_path, output_dir, output_filename_suffix, names):
    output_filename = "sequence_subset_" + output_filename_suffix + ".fastq"
    output_path = os.path.join(output_dir, output_filename)
    FileManager.ensure_directory_exists(output_dir)

    fastq_handler = FastqHandler(fastq_path)

    fastq_handler.extract_reads_by_id(id_set=names, output_path=output_path, print_status=True)

    return output_path


def partition_names(names, proportion=0.5, shuffle_names=True):
    if shuffle_names:
        random.shuffle(names)

    a = set()
    b = set()

    for n,name in enumerate(names):
        p = n/len(names)

        if p >= proportion:
            b.add(name)
        else:
            a.add(name)

    print("Partition A size: %d" % len(a))
    print("Partition B size: %d" % len(b))

    return a, b


def find_intersection_of_runnie_and_fastq(output_path, fastq_path, runnie_path):
    fastq_handler = FastqHandler(fastq_path)
    runnie_handler = RunlengthHandler(runnie_path)

    guppy_iterator = fastq_handler.iterate_read_names(print_status=True)
    runnie_iterator = runnie_handler.iterate_read_names(print_status=True)

    cutoff = sys.maxsize
    # cutoff = 100000

    guppy_names = set()
    for n,name in enumerate(guppy_iterator):
        name = name.split(" ")[0]
        guppy_names.add(name)

        if n == cutoff - 1:
            break

    print()
    runnie_names = set()
    for n,name in enumerate(runnie_iterator):
        runnie_names.add(name)

        if n == cutoff - 1:
            break

    print()

    names = guppy_names & runnie_names

    print("Total intersection size: %d" % len(names))

    print("SAVING FILE: %s" % output_path)
    with open(output_path, "w") as file:
        for name in names:
            file.write("%s\n" % name)

    return output_path


def read_names_from_fastq(fastq_path, cutoff=sys.maxsize):
    fastq_handler = FastqHandler(fastq_path)

    iterator = fastq_handler.iterate_read_names(print_status=True)

    names = set()
    for n,name in enumerate(iterator):
        name = name.split(" ")[0]
        names.add(name)

        if n == cutoff - 1:
            break

    return names


def main():
    output_dir = "output/" + "read_names_" + FileManager.get_datetime_string()
    output_filename = "read_names.txt"
    output_path = os.path.join(output_dir, output_filename)
    FileManager.ensure_directory_exists(output_dir)

    # STEP 1
    # Find union of read names within runnie and fastq files
    fastq_path = "/home/ryan/data/Nanopore/ecoli/guppy/r94_ec_guppy_rad2.fastq"
    runnie_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_all.out"

    # name_intersection_path = find_intersection_of_runnie_and_fastq(output_path=output_path,
    #                                                                fastq_path=fastq_path,
    #                                                                runnie_path=runnie_path)

    # STEP 2
    # Split sequence names into train/test partition
    name_intersection_path = "/home/ryan/code/runlength_analysis/output/read_names_2019_3_26_11_50_guppy_runnie_intersection/read_names.txt"
    names = read_names_from_file(name_intersection_path)
    names_train, names_test = partition_names(names)

    # STEP 3
    # Extract names and write to files
    runnie_train_subset_path = extract_runnie_reads_by_name(runnie_path=runnie_path,
                                                            output_dir=output_dir,
                                                            output_filename_suffix="train",
                                                            names=names_train)

    fastq_train_subset_path = extract_fastq_reads_by_name(fastq_path=fastq_path,
                                                          output_dir=output_dir,
                                                          output_filename_suffix="train",
                                                          names=names_train)

    runnie_test_subset_path = extract_runnie_reads_by_name(runnie_path=runnie_path,
                                                           output_dir=output_dir,
                                                           output_filename_suffix="test",
                                                           names=names_test)

    fastq_test_subset_path = extract_fastq_reads_by_name(fastq_path=fastq_path,
                                                         output_dir=output_dir,
                                                         output_filename_suffix="test",
                                                         names=names_test)

    # STEP 4
    # Verify
    name_intersection_path = find_intersection_of_runnie_and_fastq(output_path=output_path,
                                                                   fastq_path=fastq_train_subset_path,
                                                                   runnie_path=runnie_train_subset_path)

    name_intersection_path = find_intersection_of_runnie_and_fastq(output_path=output_path,
                                                                   fastq_path=fastq_test_subset_path,
                                                                   runnie_path=runnie_test_subset_path)


def find_runnie_read_line_number():
    runnie_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_all.out"

    runnie_handler = RunlengthHandler(runnie_path)
    runnie_iterator = runnie_handler.iterate_read_names(print_status=True, line_number=True)

    cutoff = sys.maxsize

    for n, (line_number, name) in enumerate(runnie_iterator):
        if name == "76608869-fe22-4b1a-84c6-412166583600":
            print()
            print(line_number, name)
            print()
            break


def find_fastq_read_line_number():
    fastq_path = "/home/ryan/data/Nanopore/ecoli/guppy/r94_ec_guppy_rad2.fastq"

    fastq_handler = FastqHandler(fastq_path, strip_header=True)
    fastq_iterator = fastq_handler.iterate_read_names(print_status=True, line_number=True)

    cutoff = sys.maxsize

    for n, (line_number, name) in enumerate(fastq_iterator):
        if name == "76608869-fe22-4b1a-84c6-412166583600":
            print()
            print(line_number, name)
            print()
            break


def print_file_lines(file_path, start, stop):
    with open(file_path, "r") as file:
        for l,line in enumerate(file):
            if start <= l <= stop:
                print(line)


if __name__ == "__main__":
    # main()

    names = read_names_from_fastq("/home/ryan/code/runlength_analysis/output/guppy_vs_runnie_ecoli_rad2_train_test_sequences/sequence_subset_test_60x_10kb.fastq")
    runnie_path = "/home/ryan/code/runlength_analysis/output/guppy_vs_runnie_ecoli_rad2_train_test_sequences/runnie_subset_test.out"
    output_dir = "/home/ryan/code/runlength_analysis/output/guppy_vs_runnie_ecoli_rad2_train_test_sequences/"
    output_filename_suffix = "test_60x_10kb"
    extract_runnie_reads_by_name(runnie_path, output_dir, output_filename_suffix, names)

    # find_read_line_number()
    # find_fastq_read_line_number()

    # print_file_lines("/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_all.out", 1264104048, 1264104048+100)
    # print_file_lines("/home/ryan/data/Nanopore/ecoli/guppy/r94_ec_guppy_rad2.fastq", 219100, 219100+6)
