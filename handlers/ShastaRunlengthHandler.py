import sys
import os


class Read:
    def __init__(self, consensus_sequence, consensus_lengths, pileup):
        self.sequence = consensus_sequence
        self.lengths = consensus_lengths
        self.pileup = pileup


class ReadData:
    def __init__(self, base, length, reversal, count):
        self.base = base
        self.length = length
        self.reversal = reversal
        self.count = count

    def reversal_to_string(self, reversal):
        if reversal:
            s = "R"
        else:
            s = "F"

        return s

    def __str__(self):
        return str(self.base) + str(self.length) + self.reversal_to_string(self.reversal) + str(self.count)


class ShastaRunlengthHandler:
    def __init__(self, file_path):
        self.file_path = file_path
        self.name = self.get_name_from_file_path(file_path)
        self.consensus_sequence = list()
        self.consensus_lengths = list()
        self.pileup_data = list()

    def get_name_from_file_path(self, file_path):
        name = os.path.basename(file_path).split(".")[0]

        return name

    def get_pileup_data(self, cutoff=sys.maxsize):
        self.read_file(cutoff=cutoff)

        read = Read(consensus_sequence="".join(self.consensus_sequence),
                    consensus_lengths=self.consensus_lengths,
                    pileup=self.pileup_data)

        return self.name, read

    def parse_strand(self, strand):
        if strand == "-":
            reversal = True
        else:
            reversal = False

        return reversal

    def parse_line(self, line):
        line = line.strip().split(",")
        index = line[0]
        base = line[1]
        length = int(line[2])

        self.consensus_sequence.append(base)
        self.consensus_lengths.append(length)

        read_data = line[3:-1]

        pileup = list()
        for item in read_data:
            item = item.split(" ")

            base = item[0][0]
            length = int(item[0][1:-1])
            reversal = self.parse_strand(item[0][-1])
            count = int(item[1])

            pileup.append(ReadData(base=base, length=length, reversal=reversal, count=count))

        self.pileup_data.append(pileup)

    def read_file(self, cutoff=sys.maxsize):
        with open(self.file_path, "r") as file:
            for l,line in enumerate(file):
                try:
                    self.parse_line(line)
                except Exception as e:
                    print(e)
                    exit("ERROR on line %d in file %s"%(l, self.file_path))

                if l == cutoff:
                    break
