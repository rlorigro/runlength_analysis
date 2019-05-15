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
    output_directory = "output/runnie_hdf5/"
    output_file_name = "runnie_" + FileManager.get_datetime_string() + ".hdf5"
    output_file_path = os.path.join(output_directory, output_file_name)

    FileManager.ensure_directory_exists(output_directory)

    runnie_file_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_all.out"
    runnie_file = RunlengthHandler(runnie_file_path)

    h5_file = h5py.File(output_file_path, "w")

    r = -1
    for r,read in enumerate(runnie_file.iterate_file()):
        if r % 100 == 0:
            sys.stderr.write("\rReads converted: %d" % (r+1))

        write_read_to_h5(h5_file, read)

    sys.stderr.write("\rReads converted: %d" % (r+1))
    sys.stderr.write("\n")

    h5_file.close()


if __name__ == "__main__":
    main()
