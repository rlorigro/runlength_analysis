from discrete_weibull_distribution import *
import numpy
import sys


class Read:
    def __init__(self, read_id, sequence, scales, shapes):
        self.id = read_id
        self.scales = scales
        self.shapes = shapes
        self.length = len(sequence)

        if type(sequence) is list:
            self.sequence = ''.join(sequence)
        else:
            self.sequence = sequence

class RunlengthHandler:
    def __init__(self, path):
        self.path = path
        self.read_id = None
        self.sequence = list()
        self.scales = list()
        self.shapes = list()
        self.buffered_read = None
        self.n = 0

    def extract_reads_by_id(self, id_set, output_path, print_status=False, sequence_cutoff=None):
        if sequence_cutoff is None:
            sequence_cutoff = len(id_set)

        file = open(self.path, "r")
        out_file = open(output_path, "w")

        valid_read = False
        read_lines = list()
        n_sequences = 0
        n_found = 0
        for l,line in enumerate(file):
            # Header
            if line[0] == "#":
                n_sequences += 1

                if valid_read and len(read_lines) > 0:
                    n_found += 1

                    # Reset flag
                    valid_read = False

                    # Write the buffer
                    for read_line in read_lines:
                        out_file.write(read_line)

                    read_lines = list()

                if print_status and l%1000:
                    sys.stderr.write("\r %d %d" % (n_sequences,n_found))

                read_id = " ".join(line.strip().split(" ")[1:])

                if read_id in id_set:
                    valid_read = True

                if n_found == sequence_cutoff:
                    break

            if valid_read:
                read_lines.append(line)

        if valid_read and len(read_lines) > 0:
            # Write the buffer
            for read_line in read_lines:
                out_file.write(read_line)

        sys.stderr.write("\n")
        file.close()
        out_file.close()

    def iterate_read_names(self, print_status=False, line_number=False):
        file = open(self.path, "r")
        n_reads = 0

        for l,line in enumerate(file):
            # Header
            if line[0] == "#":
                read_id = " ".join(line.strip().split(" ")[1:])
                n_reads += 1

                if print_status:
                    sys.stderr.write("\r %d" % n_reads)

                if line_number:
                    yield l, read_id

                else:
                    yield read_id

        sys.stderr.write("\n")

    def iterate_file(self, sequence_cutoff=sys.maxsize, print_status=False):
        if print_status:
            sys.stderr.write("Reading runnie file...\n")

        file = open(self.path, "r")

        for l,line in enumerate(file):
            if l%1000 == 0 and print_status:
                sys.stderr.write("\r line: %d" % l)

            # Header
            if line[0] == "#":
                self.parse_header_line(line)

            # Data
            elif line[0].isalpha():
                self.parse_data_line(line)

            # Unrecognized
            else:
                raise Exception("\nERROR: incorrect format detected at line %d in file %s" % (l, self.path))

            # Generate read data whenever the end of a read is reached
            if self.buffered_read is not None:
                yield self.buffered_read
                self.buffered_read = None

            if self.n == sequence_cutoff:
                break

        file.close()

        if print_status:
            sys.stderr.write("\nCompleted\n")

        self.__init__(self.path)

    def parse_header_line(self, line):
        # Add previous (completed) read to buffer slot
        if self.read_id is not None:
            read = Read(self.read_id, self.sequence, self.scales, self.shapes)
            self.buffered_read = read

            self.n += 1

        # Parse new read
        self.read_id = " ".join(line.strip().split(" ")[1:])
        self.sequence = list()
        self.shapes = list()
        self.scales = list()

    def parse_data_line(self, line):
        line = line.strip().split("\t")

        # print(line)
        base = line[0]
        shape = float.fromhex(line[1])
        scale = float.fromhex(line[2])

        self.sequence.append(base)
        self.shapes.append(shape)
        self.scales.append(scale)

    @staticmethod
    def calculate_numerical_mode(scale, shape, x_range):
        if shape > 0 and scale > 0:
            y = evaluate_discrete_weibull(shape=shape, scale=scale, x=x_range)

            # Get analytical mode for the continuous weibull using parameters
            mode = calculate_mode(scale=scale, shape=shape)

            # Generate window of +1 -1 around analytical mode
            min_index = max(0, round(mode) - 1)
            max_index = min_index + 2

            # Find numerical mode within window
            mode_numerical = min_index + numpy.argmax(y[min_index:max_index])

            # Add one because runnie is 0-based, and 0 is not a valid runlength
            mode_numerical += 1

        else:
            # NULL
            mode_numerical = 0

        return mode_numerical

    @staticmethod
    def convert_runnie_data_to_rle_sequence(data, max_runlength=50):
        sequence = list()
        lengths = list()

        x = numpy.arange(0, max_runlength)

        for item in data:
            mode_numerical = RunlengthHandler.calculate_numerical_mode(scale=item.scale, shape=item.shape, x_range=x)

            sequence.append(item.base)
            lengths.append(mode_numerical)

        sequence = ''.join(sequence)

        return sequence, lengths

    @staticmethod
    def expand_runlength_sequence(sequence, lengths):
        expanded_sequence = list()
        for i in range(len(sequence)):
            character = sequence[i]
            length = lengths[i]

            expanded_sequence.extend([character] * length)

        expanded_sequence = ''.join(expanded_sequence)

        return expanded_sequence


