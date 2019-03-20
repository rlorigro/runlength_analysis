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
    histogram /= sum(histogram)
    mse = 0

    for i in range(len(histogram)):
        probability = histogram[i]

        mse += (target - probability*i)**2

    return mse


def load_base_probability_matrix_from_csv(path):
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


def main():
    csv_paths = ["/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_v1/frequency_matrices_2019_3_8_15_10_19_914103.csv",
                 "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_v2/frequency_matrices_2019_3_8_15_0_34_372610.csv",
                 "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_vMode/frequency_matrices_2019_3_12_13_1_22_387365.csv"]

    # csv_paths = ["/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_vMode_0_ONLY/frequency_matrices_2019_3_12_12_25_46_220680.csv",
    #              "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_v2_0_ONLY/frequency_matrices_2019_3_12_12_24_6_369476.csv"]

    names = ["mean (v1)", "mean (v2)", "mode"]

    accuracies = dict()
    weighted_distances = dict()
    for path in csv_paths:
        matrix = load_base_probability_matrix_from_csv(path)
        matrix = numpy.sum(matrix, axis=0)

        diagonal = numpy.eye(matrix.shape[0], matrix.shape[1])

        diagonal_mask = (diagonal == 1)
        off_diagonal_mask = (diagonal == 0)

        mask_template = numpy.zeros(matrix.shape, dtype=numpy.bool)

        row_accuracies = list()
        row_weighted_distances = list()

        for i in range(matrix.shape[0]):
            row_mask = numpy.copy(mask_template)
            row_mask[i,:] = numpy.True_

            row_diagonal_confusion = matrix[diagonal_mask & row_mask]
            row_off_diagonal_confusion = matrix[off_diagonal_mask & row_mask]

            row_diagonal_sum = numpy.sum(row_diagonal_confusion)
            row_off_diagonal_sum = numpy.sum(row_off_diagonal_confusion)
            row_sum = numpy.sum(matrix[row_mask])

            print(matrix[row_mask][:10])

            row_accuracies.append(float(row_diagonal_sum/row_sum))

            weighted_average = calculate_mean(matrix[row_mask])
            distance = weighted_average - i
            row_weighted_distances.append(distance)

            print(calculate_mse(matrix[row_mask], target=i))
            print(i, row_diagonal_sum/row_sum, row_off_diagonal_sum/row_sum)

        diagonal_confusion = matrix[diagonal_mask]
        off_diagonal_confusion = matrix[off_diagonal_mask]

        diagonal_sum = numpy.sum(diagonal_confusion)
        off_diagonal_sum = numpy.sum(off_diagonal_confusion)

        matrix_sum = numpy.sum(matrix)
        accuracies[path] = row_accuracies
        weighted_distances[path] = row_weighted_distances

        # print()
        # print(path)
        # print(diagonal_confusion.shape, off_diagonal_confusion.shape)
        # print(diagonal_sum, off_diagonal_sum)
        # print(diagonal_sum/matrix_sum, off_diagonal_sum/matrix_sum)

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

    axes = pyplot.axes()
    for k,key in enumerate(weighted_distances):
        row_weighted_distances = weighted_distances[key]

        # axes.bar(numpy.arange(0,51,1), row_accuracies, width=1, align="center", alpha=0.5)
        axes.plot(numpy.arange(0,51,1), row_weighted_distances, alpha=0.8, lw=0.5)
        axes.scatter(numpy.arange(0,51,1), row_weighted_distances, alpha=0.5)

    axes.legend(names, loc="upper right")
    axes.set_title("Mean residual error")
    axes.set_ylabel("Average distance from true length")
    axes.set_xlabel("Length")

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
