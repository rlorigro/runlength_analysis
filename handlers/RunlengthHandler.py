from discrete_weibull_distribution import *
import numpy
import sys


class RunlengthData:
    def __init__(self, base, scale, shape):
        self.base = base
        self.scale = scale
        self.shape = shape


class Read:
    def __init__(self, read_id, sequence_data):
        self.id = read_id
        self.data = sequence_data


class RunlengthHandler:
    def __init__(self, path):
        self.path = path
        self.read_id = None
        self.sequence = list()
        self.buffered_read = None
        self.n = 0

    def iterate_file(self, sequence_cutoff=sys.maxsize, print_status=False):
        if print_status:
            sys.stderr.write("\tReading runnie file...\n")

        file = open(self.path, "r")

        for l,line in enumerate(file):
            if l%1000 == 0 and print_status:
                sys.stderr.write("\r\t line: %d" % l)

            # Header
            if line[0] == "#":
                self.parse_header_line(line)

            # Data
            elif line[0].isalpha():
                self.parse_data_line(line)

            # Unrecognized
            else:
                raise Exception("ERROR: incorrect format detected at line %d in file %s" % (l, self.path))

            # Generate read data whenever the end of a read is reached
            if self.buffered_read is not None:
                yield self.buffered_read
                self.buffered_read = None   # empty the cache containing read

            if self.n == sequence_cutoff:
                break

        file.close()

        if print_status:
            sys.stderr.write("\n\tCompleted\n")

    def parse_header_line(self, line):
        # Add previous (completed) read to buffer slot
        if self.read_id is not None:
            read = Read(self.read_id, self.sequence)
            self.buffered_read = read

            self.n += 1

        # Parse new read
        self.read_id = line.strip().split(" ")[-1]
        self.sequence = list()

    def parse_data_line(self, line):
        line = line.strip().split("\t")

        # print(line)
        base = line[0]
        shape = float.fromhex(line[1])
        scale = float.fromhex(line[2])

        data = RunlengthData(base, scale, shape)
        self.sequence.append(data)

    @staticmethod
    def convert_runnie_data_to_rle_sequence(data, max_runlength=50):
        sequence = list()
        lengths = list()

        x = numpy.arange(0, max_runlength, dtype="f4")

        for item in data:
            y = evaluate_discrete_weibull(shape=item.shape, scale=item.scale, x=x)

            # Get analytical mode for the continuous weibull using parameters
            mode = calculate_mode(scale=item.scale, shape=item.shape)

            # Generate window of +1 -1 around analytical mode
            min_index = max(0, round(mode) - 1)
            max_index = min_index + 2

            # Find numerical mode within window
            mode_numerical = min_index + numpy.argmax(y[min_index:max_index])

            # Add one because runnie is 0-based, and 0 is not a valid runlength
            mode_numerical += 1

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


