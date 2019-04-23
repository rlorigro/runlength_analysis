from handlers.FileManager import FileManager
from matplotlib import pyplot, colors
import matplotlib
import platform
import numpy
import sys
import os

if os.environ.get("DISPLAY", "") == "":
   print("no display found. Using non-interactive Agg backend")
   matplotlib.use("Agg")

if platform.system() == "Darwin":
   matplotlib.use("macosx")

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


def load_base_length_matrix_from_csv(path, max_runlength):
    matrices = numpy.zeros([4, max_runlength + 1, max_runlength + 1])
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

                    matrices[base_index, row_index, :] = row

                    row_index += 1

                else:
                    # Space
                    is_data = False

    matrices = matrices

    return matrices


def print_frequency_matrices(frequency_matrices, cutoff=None):
    if cutoff is None:
        cutoff = frequency_matrices.shape[-1]

    if frequency_matrices.ndim == 4:
        for i in range(frequency_matrices.shape[1]):
            print(numpy.array2string(frequency_matrices[0,i,:cutoff,:cutoff], max_line_width=sys.maxsize, separator='\t', formatter={'int': lambda x: "%d"%x}))
    elif frequency_matrices.ndim == 3:
        for i in range(frequency_matrices.shape[0]):
            print(numpy.array2string(frequency_matrices[i,:cutoff,:cutoff], max_line_width=sys.maxsize, separator='\t', formatter={'int': lambda x: "%d"%x}))


def plot_directional_complementary_diff(frequency_matrices, max_runlength, cutoff=20):
    frequency_matrices = diff_complementary_matrices(frequency_matrices, max_runlength=max_runlength)
    # print_frequency_matrices(frequency_matrices=frequency_matrices, cutoff=10)

    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    vmax = numpy.max(frequency_matrices)
    vmin = -vmax

    colormap = pyplot.cm.RdBu
    norm = colors.SymLogNorm(linthresh=0.03, vmin=vmin, vmax=vmax)

    c_axes = axes[0][0].imshow(frequency_matrices[0, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][0].set_title(INDEX_TO_BASE[0])
    figure.colorbar(c_axes, ax=axes[0][0])

    c_axes = axes[1][0].imshow(frequency_matrices[3, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][0].set_title(INDEX_TO_BASE[3])
    figure.colorbar(c_axes, ax=axes[1][0])

    c_axes = axes[0][1].imshow(frequency_matrices[2, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][1].set_title(INDEX_TO_BASE[2])
    figure.colorbar(c_axes, ax=axes[0][1])

    c_axes = axes[1][1].imshow(frequency_matrices[1, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][1].set_title(INDEX_TO_BASE[1])
    figure.colorbar(c_axes, ax=axes[1][1])

    pyplot.show()
    pyplot.close()


def plot_directional_diff(frequency_matrices, max_runlength, cutoff=20):
    frequency_matrices = diff_directional_matrices(frequency_matrices, max_runlength=max_runlength)
    # print_frequency_matrices(frequency_matrices=frequency_matrices, cutoff=10)

    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    vmax = numpy.max(frequency_matrices)
    vmin = -vmax

    colormap = pyplot.cm.RdBu
    norm = colors.SymLogNorm(linthresh=0.03, vmin=vmin, vmax=vmax)

    c_axes = axes[0][0].imshow(frequency_matrices[0, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][0].set_title(INDEX_TO_BASE[0])
    figure.colorbar(c_axes, ax=axes[0][0])

    c_axes = axes[1][0].imshow(frequency_matrices[3, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][0].set_title(INDEX_TO_BASE[3])
    figure.colorbar(c_axes, ax=axes[1][0])

    c_axes = axes[0][1].imshow(frequency_matrices[2, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[0][1].set_title(INDEX_TO_BASE[2])
    figure.colorbar(c_axes, ax=axes[0][1])

    c_axes = axes[1][1].imshow(frequency_matrices[1, :cutoff, :cutoff], cmap=colormap, norm=norm)
    axes[1][1].set_title(INDEX_TO_BASE[1])
    figure.colorbar(c_axes, ax=axes[1][1])

    pyplot.show()
    pyplot.close()


def diff_complementary_matrices(directional_matrices, max_runlength):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by subtracting the counts of any reverse bases from their base complement's matrix
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, max_runlength+1, max_runlength+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        if r == 0:
            reversal = False
        else:
            reversal = True

        for b in range(n_bases):
            if reversal:
                destination_index = complement_base_index(b)
                nondirectional_matrices[destination_index,:,:] -= directional_matrices[r,b,:,:]
            else:
                destination_index = b
                nondirectional_matrices[destination_index,:,:] += directional_matrices[r,b,:,:]

    return nondirectional_matrices


def diff_directional_matrices(directional_matrices, max_runlength):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by subtracting the counts of any reverse bases from their base complement's matrix
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, max_runlength+1, max_runlength+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        if r == 0:
            reversal = False
        else:
            reversal = True

        for b in range(n_bases):
            if reversal:
                destination_index = b
                nondirectional_matrices[destination_index,:,:] -= directional_matrices[r,b,:,:]
            else:
                destination_index = b
                nondirectional_matrices[destination_index,:,:] += directional_matrices[r,b,:,:]

    return nondirectional_matrices


def plot_base_matrices(matrix, test_spot=False, cutoff=20, normalize_matrices=False):
    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    if matrix.ndim == 4:
        matrix_A = matrix[0, 0, :cutoff, :cutoff]

        if test_spot:
            matrix_A[cutoff-1,0] = numpy.max(numpy.log10(matrix_A)) + 1

        if normalize_matrices:
            matrix_A = normalize(matrix_A, pseudocount=0)
            axes[0][0].imshow(matrix_A)
        else:
            axes[0][0].imshow(numpy.log10(matrix_A))

        axes[0][0].set_title(INDEX_TO_BASE[0])

        matrix_T = matrix[0, 3, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_T  = normalize(matrix_T, pseudocount=0)
            axes[1][0].imshow(matrix_T)
        else:
            axes[1][0].imshow(numpy.log10(matrix_T))

        axes[1][0].set_title(INDEX_TO_BASE[3])

        matrix_G = matrix[0, 2, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_G  = normalize(matrix_G, pseudocount=0)
            axes[0][1].imshow(matrix_G)
        else:
            axes[0][1].imshow(numpy.log10(matrix_G))

        axes[0][1].set_title(INDEX_TO_BASE[2])

        matrix_C = matrix[0, 1, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_C  = normalize(matrix_C, pseudocount=0)
            axes[1][1].imshow(matrix_C)
        else:
            axes[1][1].imshow(numpy.log10(matrix_C))

        axes[1][1].set_title(INDEX_TO_BASE[1])

    elif matrix.ndim == 3:
        matrix_A = matrix[0, :cutoff, :cutoff]

        if test_spot:
            matrix_A[cutoff, 0] = numpy.max(matrix_A)

        if normalize_matrices:
            matrix_A = normalize(matrix_A, pseudocount=0)
            axes[0][0].imshow(matrix_A)
        else:
            axes[0][0].imshow(numpy.log10(matrix_A))

        axes[0][0].set_title(INDEX_TO_BASE[0])

        matrix_T = matrix[3, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_T  = normalize(matrix_T, pseudocount=0)
            axes[1][0].imshow(matrix_T)
        else:
            axes[1][0].imshow(numpy.log10(matrix_T))
        axes[1][0].set_title(INDEX_TO_BASE[3])

        matrix_G = matrix[2, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_G  = normalize(matrix_G, pseudocount=0)
            axes[0][1].imshow(matrix_G)
        else:
            axes[0][1].imshow(numpy.log10(matrix_G))
        axes[0][1].set_title(INDEX_TO_BASE[2])

        matrix_C = matrix[1, :cutoff, :cutoff]

        if normalize_matrices:
            matrix_C  = normalize(matrix_C, pseudocount=0)
            axes[1][1].imshow(matrix_C)
        else:
            axes[1][1].imshow(numpy.log10(matrix_C))
        axes[1][1].set_title(INDEX_TO_BASE[1])

    axes[1][1].set_xlabel("Observed length")
    axes[1][0].set_xlabel("Observed length")
    axes[1][0].set_ylabel("True length")
    axes[0][0].set_ylabel("True length")

    axes[1][1].set_yticks(numpy.arange(0,cutoff, 2))
    axes[1][0].set_yticks(numpy.arange(0,cutoff, 2))
    axes[0][0].set_yticks(numpy.arange(0,cutoff, 2))
    axes[0][1].set_yticks(numpy.arange(0,cutoff, 2))

    axes[1][1].set_xticks(numpy.arange(0, cutoff, 2))
    axes[1][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][1].set_xticks(numpy.arange(0, cutoff, 2))

    pyplot.show()
    pyplot.close()


def sum_complementary_matrices(directional_matrices, max_runlength):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by adding the counts of any reverse bases into their base complement's matrix
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, max_runlength+1, max_runlength+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        if r == 0:
            reversal = False
        else:
            reversal = True

        for b in range(n_bases):
            destination_index = b
            if reversal:
                destination_index = complement_base_index(b)

            nondirectional_matrices[destination_index,:,:] += directional_matrices[r,b,:,:].squeeze()

    return nondirectional_matrices


def sum_directional_matrices(directional_matrices, max_runlength):
    """
    Take a frequency matrix with shape [n_directions, n_bases, n_y, n_x] and convert it to a non-directional matrix,
    by adding the counts of any reverse bases into their forward base's matrix (WITHOUT COMPLEMENTING!)
    :param directional_matrices:
    :return:
    """
    nondirectional_matrices = numpy.zeros([4, max_runlength+1, max_runlength+1], dtype=numpy.float64)

    n_directions, n_bases, n_y, n_x = directional_matrices.shape

    for r in range(n_directions):
        for b in range(n_bases):
            destination_index = b

            nondirectional_matrices[destination_index,:,:] += directional_matrices[r,b,:,:].squeeze()

    return nondirectional_matrices


def complement_base_index(base_index):
    return 3 - base_index


def normalize(frequency_matrix, pseudocount, diagonal_bias=0):
    frequency_matrix = frequency_matrix.astype(numpy.float32)

    frequency_matrix += pseudocount

    diagonal_pseudocount = diagonal_bias
    diagonal_mask = numpy.eye(frequency_matrix.shape[0], frequency_matrix.shape[1], dtype=numpy.bool)
    frequency_matrix[diagonal_mask] += diagonal_pseudocount

    sum_y = numpy.sum(frequency_matrix, axis=1)

    probability_matrix = frequency_matrix / sum_y[:, numpy.newaxis]
    probability_matrix = numpy.log10(probability_matrix)

    return probability_matrix


def save_directional_frequency_matrices_as_delimited_text(output_dir, frequency_matrices, chromosome_name=None, delimiter=",", log_normalize=False, plot=False, pseudocount=0, diagonal_bias=0, default_type=int):
    if chromosome_name is not None:
        name_suffix = chromosome_name + "_"
    else:
        name_suffix = ""

    if log_normalize:
        filename = "probability_matrices_directional_" + name_suffix + FileManager.get_datetime_string() + ".csv"
    else:
        filename = "frequency_matrices_directional_" + name_suffix + FileManager.get_datetime_string() + ".csv"

    reversal_suffixes = ["F", "R"]
    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    print("SAVING: %s" % output_path)

    for reversal in [0,1]:
        for base_index in range(4):
            base = INDEX_TO_BASE[base_index]
            suffix = reversal_suffixes[reversal]

            matrix = numpy.squeeze(frequency_matrices[reversal,base_index,:,:])

            type = default_type
            if log_normalize:
                matrix = normalize(matrix, pseudocount=pseudocount, diagonal_bias=diagonal_bias)
                type = float

            if plot:
                pyplot.imshow(matrix)
                pyplot.show()
                pyplot.close()

            matrix_name = "_".join([base, suffix])
            header = ">" + matrix_name + "\n"

            # print(type)
            file.write(header)
            for r in range(matrix.shape[0]):
                row = [str(type(x)) for x in matrix[r]]

                # if r < 4 and not log_normalize:
                    # print(row)

                row = delimiter.join(row) + "\n"

                file.write(row)

            file.write("\n")

    file.close()


def save_nondirectional_frequency_matrices_as_delimited_text(output_dir, frequency_matrices, chromosome_name=None, delimiter=",", log_normalize=False, pseudocount=1e-12, diagonal_bias=0, plot=False, default_type=int):
    if chromosome_name is not None:
        name_suffix = chromosome_name + "_"
    else:
        name_suffix = ""

    if log_normalize:
        filename = "probability_matrices_" + name_suffix + FileManager.get_datetime_string() + ".csv"
    else:
        filename = "frequency_matrices_" + name_suffix + FileManager.get_datetime_string() + ".csv"

    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    for base_index in range(4):
        base = INDEX_TO_BASE[base_index]

        matrix = numpy.squeeze(frequency_matrices[base_index,:,:])

        type = default_type
        if log_normalize:
            matrix = normalize(matrix, pseudocount=pseudocount, diagonal_bias=diagonal_bias)
            type = float

        if plot:
            pyplot.imshow(matrix)
            pyplot.show()
            pyplot.close()

        matrix_name = base
        header = ">" + matrix_name + "\n"

        file.write(header)

        for r in range(matrix.shape[0]):
            row = [str(type(x)) for x in matrix[r]]
            row = delimiter.join(row) + "\n"

            file.write(row)

        file.write("\n")

    file.close()
