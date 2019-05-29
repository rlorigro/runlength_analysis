from measure_runlength_distribution_from_fasta import *
import os


def write_fasta(runlength_reference_path, runlength_ref_sequences, runlength_read_path, runlength_read_sequences):
    print("SAVING run length fasta file:", runlength_reference_path)
    print("SAVING run length fasta file:", runlength_read_path)

    with open(runlength_reference_path, "w") as file:
        for contig_name in runlength_ref_sequences.keys():
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_ref_sequences[contig_name][SEQUENCE] + "\n")

    with open(runlength_read_path, "w") as file:
        for contig_name in runlength_read_sequences.keys():
            file.write(">"+contig_name+" RLE\n")
            file.write(runlength_read_sequences[contig_name][SEQUENCE] + "\n")


def runlength_encode_fasta(ref_fasta_path, read_fasta_path, output_parent_dir="output/"):
    output_dir = "runlength_matrix_from_sequence_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    ref_fasta_filename_prefix = ".".join(os.path.basename(ref_fasta_path).split(".")[:-1])
    runlength_ref_fasta_filename = ref_fasta_filename_prefix + "_rle.fasta"
    runlength_ref_fasta_path = os.path.join(output_dir, runlength_ref_fasta_filename)

    read_fasta_filename_prefix = ".".join(os.path.basename(read_fasta_path).split(".")[:-1])
    runlength_read_fasta_filename = read_fasta_filename_prefix + "_rle.fasta"
    runlength_read_fasta_path = os.path.join(output_dir, runlength_read_fasta_filename)

    runlength_ref_sequences = runlength_encode_fasta_parallel(fasta_sequence_path=ref_fasta_path, min_length=0)
    runlength_read_sequences = runlength_encode_fasta_parallel(fasta_sequence_path=read_fasta_path, min_length=0)

    write_fasta(runlength_reference_path=runlength_ref_fasta_path,
                runlength_ref_sequences=runlength_ref_sequences,
                runlength_read_path=runlength_read_fasta_path,
                runlength_read_sequences=runlength_read_sequences)


def main(ref_fasta_path, read_fasta_path, output_dir):
    if output_dir is None:
        output_dir = "output/"

    runlength_encode_fasta(ref_fasta_path=ref_fasta_path,
                           read_fasta_path=read_fasta_path,
                           output_parent_dir=output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sequences", "-i",
        type=str,
        required=True,
        help="path to fastq file"
    )
    parser.add_argument(
        "--ref",
        type=str,
        required=True,
        help="file path of FASTA reference sequence file"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=False,
        help="path of directory to dump output to"
    )
    args = parser.parse_args()

    main(output_dir=args.output_dir, read_fasta_path=args.sequences, ref_fasta_path=args.ref)
