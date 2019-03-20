from handlers.RunlengthHandler_v2 import RunlengthHandler
from discrete_weibull_distribution import evaluate_discrete_weibull, calculate_mode
from matplotlib import pyplot
import numpy
import math
import sys


def main():
    runlength_path = "/home/ryan/data/Nanopore/ecoli/miten/runnie/runnie_phase_1_calls/rad2_pass_runnie_0.out"

    handler = RunlengthHandler(runlength_path)

    reads = handler.iterate_file(sequence_cutoff=sys.maxsize)

    x = numpy.arange(0,30)

    for r,read in enumerate(reads):
        length = read.length
        read_id = read.id
        sequence = read.sequence
        scales = read.scales
        shapes = read.shapes

        print(r, read_id)
        print([d for d in sequence[:10]])
        print([d for d in scales[:10]])
        print([d for d in shapes[:10]])

        for i in range(length):
            scale = scales[i]
            shape = shapes[i]

            mode_analytical = calculate_mode(scale=scale, shape=shape)

            x_min = max(0,round(mode_analytical)-1)
            x_max = x_min + 2
            x_subset = numpy.arange(x_min,x_max)

            y_subset = evaluate_discrete_weibull(shape=shape, scale=scale, x=x_subset)

            y = evaluate_discrete_weibull(shape=shape, scale=scale, x=x)

            # pyplot.plot(x,y)
            # pyplot.show()
            # pyplot.close()

            mode_floor = math.floor(mode_analytical)

            mode_numerical = numpy.argmax(y)

            mode_numerical_approximate = x_min + numpy.argmax(y_subset)

            if mode_numerical != mode_numerical_approximate:
                print(mode_numerical, mode_analytical, x_min, x_max)
                print(mode_floor, mode_numerical)

                axes = pyplot.axes()

                y_min, y_max = axes.get_ylim()

                axes.plot(x,y)
                axes.plot([mode_numerical_approximate,mode_numerical_approximate],[y_min,y_max], alpha=0.6)
                # axes.plot([mode_analytical,mode_floor],[y_min,y_max], alpha=0.6, linestyle="--")

                pyplot.show()
                pyplot.close()

        # if mode_floor != mode_numerical:
            #     print(mode_floor, mode_numerical)
            #
            #     axes = pyplot.axes()
            #
            #     y_min, y_max = axes.get_ylim()
            #
            #     axes.plot(x,y)
            #     axes.plot([mode_analytical,mode_analytical],[y_min,y_max], alpha=0.6)
            #     # axes.plot([mode_analytical,mode_floor],[y_min,y_max], alpha=0.6, linestyle="--")
            #
            #     pyplot.show()
            #     pyplot.close()


if __name__ == "__main__":
    main()
