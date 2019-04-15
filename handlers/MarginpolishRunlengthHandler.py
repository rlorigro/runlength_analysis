import sys
import os


class Read:
    def __init__(self, consensus_sequence, pileup):
        self.sequence = consensus_sequence
        self.pileup = pileup


class ReadData:
    def __init__(self, base, length, reversal, weight):
        self.base = base
        self.length = length
        self.reversal = reversal
        self.weight = weight

    def reversal_to_string(self, reversal):
        """
        For printing
        :param reversal:
        :return:
        """
        if reversal:
            s = "-"
        else:
            s = "+"

        return s

    def __str__(self):
        return str(self.base) + str(self.length) + self.reversal_to_string(self.reversal) + str(self.weight)


class MarginpolishRunlengthHandler:
    def __init__(self, file_path):
        self.file_path = file_path
        self.name = self.get_name_from_file_path(file_path)
        self.consensus_sequence = list()
        self.pileup_data = list()

    def get_name_from_file_path(self, file_path):
        name = os.path.basename(file_path).split(".")[0]

        return name

    def get_pileup_data(self, cutoff=sys.maxsize):
        self.read_file(cutoff=cutoff)

        # print(self.consensus_sequence)

        read = Read(consensus_sequence="".join(self.consensus_sequence),
                    pileup=self.pileup_data)

        return self.name, read

    def parse_strand(self, strand):
        if strand == "-":
            reversal = True
        else:
            reversal = False

        return reversal

    def parse_observation_string(self, observation_string):
        observation_string, weight = observation_string.split(",")
        weight = float(weight)
        base = observation_string[0]
        reversal = self.parse_strand(observation_string[1])
        length = int(observation_string[2:])

        data = ReadData(base=base, length=length, reversal=reversal, weight=weight)

        return data

    def parse_line(self, line):
        data = line.split("\t")

        consensus_index = int(data[0])
        consensus_base = data[1]
        observation_strings = data[2:]

        pileup = list()
        for observation_string in observation_strings:
            data = self.parse_observation_string(observation_string)

            pileup.append(data)

        self.consensus_sequence.append(consensus_base)
        self.pileup_data.append(pileup)

    def read_file(self, cutoff=sys.maxsize):
        with open(self.file_path, "r") as file:
            for l, line in enumerate(file):
                if l <= 1:
                    continue

                self.parse_line(line)

                if l == cutoff:
                    break
