from handlers.RunlengthHandler import RunlengthHandler
from test_runlength_handler_v2 import RunlengthHandler as RunlengthHandler_v2
import numpy
import sys


def main():
    runlength_path = "/home/ryan/data/Nanopore/ecoli/miten/runnie/runnie_phase_1_calls/rad2_pass_runnie_0.out"
    x = numpy.arange(0,30)

    handler_v1 = RunlengthHandler(runlength_path)

    reads_v1 = handler_v1.iterate_file(sequence_cutoff=2)

    print("############# V1 #############")

    for r,read in enumerate(reads_v1):
        data = read.data
        read_id = read.id

        print(r, read_id)
        print([d.base for d in data[:10]])
        print([d.scale for d in data[:10]])
        print([d.shape for d in data[:10]])

    handler_v2 = RunlengthHandler_v2(runlength_path)

    reads_v2 = handler_v2.iterate_file(sequence_cutoff=2)

    print("############# V2 #############")

    for r,read in enumerate(reads_v2):
        length = read.length
        read_id = read.id
        sequence = read.sequence
        scales = read.scales
        shapes = read.shapes

        print(r, read_id)
        print([d for d in sequence[:10]])
        print([d for d in scales[:10]])
        print([d for d in shapes[:10]])


if __name__ == "__main__":
    main()
