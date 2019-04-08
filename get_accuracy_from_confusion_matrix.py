from matplotlib import pyplot
import numpy
import csv

MAX_RUNLENGTH = 50

BASE_TO_INDEX = {"A":0,
                 "C":1,
                 "G":2,
                 "T":3}

numpy.set_printoptions(suppress=True, linewidth=400)


def calculate_mean(histogram):
    mean_value = 0
    total = sum(histogram)
    for i in range(len(histogram)):
        probability = histogram[i]/total

        mean_value += probability*i

    return mean_value


def calculate_mse(histogram, target):
    """
    For a vector describing a distribution of predictions, find the mean distance from the true value (target)
    :param histogram:
    :param target:
    :return:
    """
    histogram /= sum(histogram)
    mse = 0

    for i in range(len(histogram)):
        probability = histogram[i]

        mse += (target - probability*i)**2

    return mse


def load_base_probability_matrix_from_csv(path):
    """
    Given a CSV of 4 matrices (1 for each base) each with dimensions LxL where L=max_length, load matrices into a matrix
    with shape (4,L,L) where bases are mapped to index 0,1,2,3 via global dict BASE_TO_INDEX
    :param path:
    :return:
    """
    matrices = numpy.zeros([4,MAX_RUNLENGTH+1,MAX_RUNLENGTH+1])
    base_index = None
    row_index = 0

    with open(path, "r") as file:
        is_data = False

        for line in file:
            if not is_data:
                if line[0] == ">":
                    # Header
                    base = line[1]
                    base_index = BASE_TO_INDEX[base]

                    is_data = True
                    row_index = 0
            else:
                if not line[0].isspace():
                    # Data
                    row = list(map(float, line.strip().split(",")))

                    matrices[base_index,row_index,:] = row

                    row_index += 1

                else:
                    # Space
                    is_data = False

    matrices = matrices

    return matrices


def plot_classification_accuracy(accuracies, names):
    """
    :param accuracies: a dict of lists, where lists contain a series of classification accuracies for each class
    :param names:
    :return:
    """
    axes = pyplot.axes()
    for k,key in enumerate(accuracies):
        row_accuracies = accuracies[key]

        # axes.bar(numpy.arange(0,51,1), row_accuracies, width=1, align="center", alpha=0.5)
        axes.plot(numpy.arange(0,51,1), row_accuracies, alpha=0.8, lw=0.5)
        axes.scatter(numpy.arange(0,51,1), row_accuracies, alpha=0.5)

    axes.legend(names, loc="upper right")
    axes.set_title("Classification accuracy per length")
    axes.set_ylabel("Proportion correct")
    axes.set_xlabel("Length")

    axes.set_ylim([0,1])
    axes.set_yticks(numpy.arange(0,1,0.1))
    pyplot.show()
    pyplot.close()


def plot_mean_residual_error(weighted_distances, names):
    """
    :param weighted_distances: a dict of lists, where lists contain a series of MSEs for each class
    :param names:
    :return:
    """
    axes = pyplot.axes()
    for k,key in enumerate(weighted_distances):
        row_weighted_distances = weighted_distances[key]
        print(key, row_weighted_distances[:10])

        # axes.bar(numpy.arange(0,51,1), row_accuracies, width=1, align="center", alpha=0.5)
        axes.plot(numpy.arange(0,51,1), row_weighted_distances, alpha=0.8, lw=0.5)
        axes.scatter(numpy.arange(0,51,1), row_weighted_distances, alpha=0.5)

    axes.legend(names, loc="upper right")
    axes.set_title("Mean residual error")
    axes.set_ylabel("Average distance from true length")
    axes.set_xlabel("Length")

    pyplot.show()
    pyplot.close()


def plot_residual_error(row_centered_confusions, names, y_min, y_max):
    colors = ["blue", "orange"]

    axes = pyplot.axes()

    for k,key in enumerate(row_centered_confusions):
        confusions = row_centered_confusions[key][1:10]
        sampled_confusions = list()

        for confusion in confusions:
            # print(confusion[40:60])
            confusion[numpy.isnan(confusion)] = 0
            confusion = numpy.random.multinomial(1, confusion, size=10000)
            confusion = numpy.argmax(confusion, axis=1)
            sampled_confusions.append(confusion)

        axes.violinplot(sampled_confusions, range(0,len(confusions)),
                        # points=200,
                        widths=0.7,
                        # showmeans=True,
                        showextrema=True,
                        # showmedians=True,
                        bw_method=0.6)
                        # color=colors[k])

    axes.axhline(y_max, linestyle="--", linewidth=0.5, color="k", alpha=0.3)

    # axes.legend(names, loc="upper right")
    axes.set_title("Residual error")
    axes.set_ylabel("Distance from true length")
    axes.set_xlabel("Length")

    axes.set_yticks(list(range(0, abs(y_max)+abs(y_min))))
    axes.set_yticklabels(list(range(y_min, y_max)))
    axes.set_ylim([y_max-10,y_max+10])

    pyplot.show()
    pyplot.close()


def main():
    # csv_paths = ["/home/ryan/code/runlength_analysis/output/confusion_matrix_from_shasta_ecoli_60x_test_NO_BAYES/frequency_matrices_2019_4_2_9_33_38_540396.csv",
    #              "/home/ryan/code/runlength_analysis/output/confusion_matrix_from_shasta_ecoli_60x_test_BAYES/frequency_matrices_2019_4_2_9_31_35_82369.csv",
    #              "/home/ryan/code/runlength_analysis/output/confusion_matrix_from_shasta_ecoli_60x_test_BAYES_NO_PRIOR/frequency_matrices_2019_4_2_13_37_43_678144.csv"]

    # csv_paths = ["/home/ryan/code/runlength_analysis/output/shasta_test/NO_BAYES/frequency_matrices_2019_4_4_16_30_23_984547.csv",
    #              "/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_NEW/frequency_matrices_2019_4_4_16_31_52_691889.csv",
    #              "/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_STRANDED/frequency_matrices_2019_4_5_13_0_31_314707.csv"]

    csv_paths = ["/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_STRANDED/frequency_matrices_2019_4_5_13_0_31_314707.csv",
                 # "/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_STRANDED_NO_PRIOR/frequency_matrices_2019_4_5_13_24_29_183724.csv",
                 "/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_STRANDED_FIX/frequency_matrices_2019_4_5_13_49_19_144874.csv",
                 "/home/ryan/code/runlength_analysis/output/shasta_test/BAYES_STRANDED_FIX_NO_PRIOR/frequency_matrices_2019_4_5_14_1_10_425455.csv"]

    # csv_paths = ["/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_v1/frequency_matrices_2019_3_8_15_10_19_914103.csv",
    #              "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_v2/frequency_matrices_2019_3_8_15_0_34_372610.csv",
    #              "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_vMode/frequency_matrices_2019_3_12_13_1_22_387365.csv"]

    # names = ["mean (v1)", "mean (v2)", "mode"]
    # names = ["No Bayes", "Bayes"]
    names = ["No_Bayes", "Bayes_new", "Bayes_default"]

    accuracies = dict()
    weighted_distances = dict()
    row_confusions = dict()
    for path in csv_paths:
        matrix = load_base_probability_matrix_from_csv(path)
        matrix = numpy.sum(matrix, axis=0)

        diagonal = numpy.eye(matrix.shape[0], matrix.shape[1])

        diagonal_mask = (diagonal == 1)
        off_diagonal_mask = (diagonal == 0)

        mask_template = numpy.zeros(matrix.shape, dtype=numpy.bool)

        row_accuracies = list()             # True/non-true predictions for each class
        row_weighted_distances = list()     # MSE of predictions for each class
        row_centered_confusions = list()    # zero-centered confusion (target=0) for each true value

        # Accumulate true positives and non true positives for all true values
        for i in range(matrix.shape[0]):
            row_mask = numpy.copy(mask_template)
            row_mask[i,:] = numpy.True_

            row_diagonal_confusion = matrix[diagonal_mask & row_mask]
            row_off_diagonal_confusion = matrix[off_diagonal_mask & row_mask]

            row_diagonal_sum = numpy.sum(row_diagonal_confusion)
            row_off_diagonal_sum = numpy.sum(row_off_diagonal_confusion)
            row_sum = numpy.sum(matrix[row_mask])

            row_accuracies.append(float(row_diagonal_sum/row_sum))

            weighted_average = calculate_mean(matrix[row_mask])
            distance = weighted_average - i
            row_weighted_distances.append(distance)

            # Initialize a vector that has the capacity to contain 100% negative and 100% positive bias
            row_centered_confusion = numpy.zeros(matrix.shape[0]*2)

            # This index will effectively center the row-wise confusion on zero
            start_index = matrix.shape[0] - i
            stop_index = start_index + matrix.shape[0]

            row_centered_confusion[start_index:stop_index] = matrix[row_mask]/row_sum
            row_centered_confusions.append(row_centered_confusion)

        diagonal_confusion = matrix[diagonal_mask]
        off_diagonal_confusion = matrix[off_diagonal_mask]

        diagonal_sum = numpy.sum(diagonal_confusion)
        off_diagonal_sum = numpy.sum(off_diagonal_confusion)

        matrix_sum = numpy.sum(matrix)
        accuracies[path] = row_accuracies
        weighted_distances[path] = row_weighted_distances
        row_confusions[path] = row_centered_confusions

        print()
        print(path)
        print(diagonal_sum/matrix_sum, off_diagonal_sum/matrix_sum)

    # plot_mean_residual_error(weighted_distances, names)
    plot_residual_error(row_confusions, names, y_min=-matrix.shape[0], y_max=matrix.shape[0])


if __name__ == "__main__":
    main()
