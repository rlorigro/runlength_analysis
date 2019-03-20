from handlers.RunlengthHandler import RunlengthHandler
from handlers.FileManager import FileManager
import multiprocessing
import sys
import os


def decode(read, output_dir):
    sequence, lengths = RunlengthHandler.convert_runnie_data_to_rle_sequence(read.data)
    expanded_sequence = RunlengthHandler.expand_runlength_sequence(sequence=sequence, lengths=lengths)

    output_filename = read.id + ".fasta"
    output_file_path = os.path.join(output_dir, output_filename)

    with open(output_file_path, "w") as file:
        file.write(">" + read.id + "\n")
        file.write(expanded_sequence + "\n")

    return read.id


def arg_iterator(handler, output_dir):
    for read in handler.iterate_file():
        yield read,output_dir


def arg_unpacker(args):
    read_id = decode(read=args[0], output_dir=args[1])
    return read_id


def main(max_threads=None):
    # runlength_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_runnie_0.out"
    runlength_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_runnie_0_1_10_11_12_13.out"

    output_parent_dir = "output/version_comparison/mode/"
    output_dir = "runlength_matrix_from_assembly_contigs_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    handler = RunlengthHandler(runlength_path)

    if max_threads is None:
        max_threads = max(1, multiprocessing.cpu_count()-2)

    with multiprocessing.Pool(processes=max_threads) as pool:
        for r,read_id in enumerate(pool.imap(arg_unpacker, arg_iterator(handler=handler, output_dir=output_dir))):
            sys.stdout.write("\r%d" % r)
    print()

    print("Concatenating files...")
    output_file_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=output_dir, file_extension=".fasta")

    concatenated_filename = os.path.basename(runlength_path).split(".")[0] + ".fasta"
    concatenated_file_path = os.path.join(output_dir, concatenated_filename)

    print("Saving to file: %s" % concatenated_file_path)

    FileManager.concatenate_files(file_paths=output_file_paths, output_file_path=concatenated_file_path)
    FileManager.delete_files(output_file_paths)


if __name__ == "__main__":
    main()
