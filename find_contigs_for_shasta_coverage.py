

def get_contig_names(assembly_path, min_length):
    contig_names = list()

    with open(assembly_path, "r") as file:
        header = None
        found_header = False

        for line in file:
            if line[0] == ">":
                # header
                header = line[1:].strip().split(" ")[0]
                found_header = True

            elif found_header:
                # sequence
                sequence = line.strip()
                found_header = False

                if len(sequence) > min_length:
                    contig_names.append(header)

    return contig_names


def main():
    assembly_path = "/home/ryan/software/shasta/gm24143_chr1/run_2019_4_16_14_20_20_340647/Assembly.fasta"

    names = get_contig_names(assembly_path=assembly_path, min_length=1_000_000)

    output_filename = "assemble_segments.sh"

    with open(output_filename, "w") as file:
        for name in names:
            file.write("../../build/shasta-install/bin/AssembleSegment.py %s\n" % name)


if __name__ == "__main__":
    main()
