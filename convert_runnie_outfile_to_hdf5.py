from handlers.RunlengthHandler_v2 import RunlengthHandler
from handlers.FileManager import FileManager
import numpy
import h5py
import sys
import os


def write_read_to_h5(h5_file, read):
    read_group = h5_file.create_group(read.id)

    sequence = numpy.string_(read.sequence)
    read_group.create_dataset(name="sequence", data=sequence)
    read_group.create_dataset(name="scales", data=read.scales, dtype="f")
    read_group.create_dataset(name="shapes", data=read.shapes, dtype="f")


def main():
    runnie_directory = "/home/ryan/data/Nanopore/Human/runnie/NA12878/part1_outs"

    # output_directory = "output/runnie_hdf5/"
    output_directory = runnie_directory
    output_file_name = "runnie_" + FileManager.get_datetime_string() + ".hdf5"
    output_file_path = os.path.join(output_directory, output_file_name)

    FileManager.ensure_directory_exists(output_directory)

    h5_file = h5py.File(output_file_path, "w")

    runnie_file_paths = FileManager.get_all_file_paths_by_type(runnie_directory, ".out")

    r = -1
    for file_path in runnie_file_paths:
        runnie_file = RunlengthHandler(file_path)
        reads = runnie_file.iterate_file()

        while True:
            try:
                read = next(reads)

                if r % 100 == 0:
                    sys.stderr.write("\rReads converted: %d" % (r + 1))

                write_read_to_h5(h5_file, read)
                r += 1

            except StopIteration:
                break

            except Exception as e:
                print()
                print(e)
                print("ERROR: improper format in file %s" % file_path)
                continue

        sys.stderr.write("\rReads converted: %d" % (r+1))
        sys.stderr.write("\n")

    h5_file.close()


if __name__ == "__main__":
    main()
