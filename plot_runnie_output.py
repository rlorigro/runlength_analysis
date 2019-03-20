from handlers.RunlengthHandler import RunlengthHandler
from discrete_weibull_distribution import *
import random
import math
import sys


def plot_distribution(axes, x, y, alpha=1):
    color = [0.1,0.4,1]

    axes.bar(x, y, color=color, alpha=alpha)


def main():
    runlength_path = "/home/ryan/data/Nanopore/ecoli/runnie/out/rad2_pass_runnie_0.out"

    handler = RunlengthHandler(runlength_path)

    reads = handler.iterate_file(sequence_cutoff=200, print_status=True)

    n_distributions = 0

    x = numpy.arange(0, 10)

    sys.stderr.write("Binning distributions...\n")
    distribution_bins = [list() for i in range(60)]

    for r,read in enumerate(reads):
        data = read.data
        read_id = read.id

        for i,item in enumerate(data):
            if item.shape < 1:
                sys.stderr.write("WARNING: beta less than 1\n", item.shape)

            y = evaluate_discrete_weibull(shape=item.shape, scale=item.scale, x=x)

            # Get analytical mode for the continuous weibull using parameters
            mode = calculate_mode(scale=item.scale, shape=item.shape)

            # Generate window of +1 -1 around analytical mode
            min_index = max(0, round(mode) - 1)
            max_index = min_index + 2

            # Find numerical mode within window
            mode_numerical = min_index + numpy.argmax(y[min_index:max_index])

            # true_mode = numpy.argmax(y)
            #
            # if true_mode != mode_numerical:
            #     print(true_mode, mode_numerical)

            if mode_numerical < len(distribution_bins):
                distribution_bins[mode_numerical].append([item.scale, item.shape])
                n_distributions += 1

            # print(item.scale, item.shape)
            # print(sum)
            # print(mode)

            # axes = pyplot.axes()
            #
            # plot_distribution(axes, x[:60], y[:60])
            #
            # pyplot.show()
            # pyplot.close()

            # if i == 1000:
            #     break

    n_rows = 8

    figure, axes = pyplot.subplots(nrows=n_rows)

    sys.stderr.write("Plotting...\n")

    sample_size = 1000

    for b,bin in enumerate(distribution_bins[:n_rows]):
        alpha = 1/(sample_size/3)

        bin_sample = list()

        if len(bin) > 0:
            # print(b)
            n_items = min(sample_size, len(bin))

            while len(bin_sample) < n_items:
                bin_sample.append(random.choice(bin)[:n_rows+5])

            for scale, shape in bin_sample:
                print(",".join(list(map(str,[b+1, scale, shape]))))

                y = evaluate_discrete_weibull(scale=scale, shape=shape, x=x[:n_rows+5])

                axes[b].plot(x[:n_rows+5], y, color=[0.1,0.4,1.0], alpha=alpha, linewidth=0.8)

                label = "%d" % (b+1)

                axes[b].set_ylabel(label)

    axes[n_rows-1].set_xlabel("Run length")
    axes[0].set_title("Binned length distributions")

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
