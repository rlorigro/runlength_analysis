from handlers.MarginpolishRunlengthHandler import MarginpolishRunlengthHandler
from handlers.FastaHandler import FastaHandler
from handlers.BamHandler import BamHandler
from multiprocessing import Manager, Pool, cpu_count
from modules.align import align_minimap
from modules.matrix import *
import argparse
import numpy
import copy
import sys
import os
from multithread import *

MAX_RUNLENGTH = 50

# Indexes for "observation" data tuples
OBS_BASE = 0
OBS_REVERSAL = 1
OBS_LENGTH = 2
OBS_WEIGHT = 3

# Dimension labels for matrix of runlengths
MATRIX_REVERSAL = 0
MATRIX_BASE = 1
MATRIX_TRUE_LENGTH = 2
MATRIX_OBSERVED_LENGTH = 3

# Indexing RLE tuples
SEQUENCE = 0
LENGTHS = 1

DEBUG_INDEX = 0


def main(pickle_dir, output_parent_dir):

    output_dir = "runlength_matrix_from_slurped_pickels_" + FileManager.get_datetime_string()
    output_dir = os.path.join(output_parent_dir, output_dir)
    FileManager.ensure_directory_exists(output_dir)

    matrix_file_paths = FileManager.get_all_file_paths_by_type(parent_directory_path=pickle_dir,
                                                               file_extension=".pkl")

    all_matrices = list()
    for path in matrix_file_paths:
        all_matrices.append(numpy.load(path))

    matrix = numpy.stack(all_matrices, axis=0)
    matrix = numpy.sum(matrix, axis=0)

    save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                          frequency_matrices=matrix,
                                                          chromosome_name="genomic",
                                                          log_normalize=False,
                                                          diagonal_bias=0,
                                                          plot=False,
                                                          default_type=float)

    save_directional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                          frequency_matrices=matrix,
                                                          chromosome_name="genomic",
                                                          log_normalize=True,
                                                          pseudocount=15,
                                                          diagonal_bias=0,
                                                          plot=False,
                                                          default_type=float)

    save_directional_frequency_matrices_as_marginPolish_params(output_dir=output_dir,
                                                               frequency_matrices=matrix,
                                                               chromosome_name="genomic",
                                                               pseudocount=15,
                                                               diagonal_bias=0)

    nondirectional_matrix = sum_complementary_matrices(matrix, max_runlength=MAX_RUNLENGTH)

    save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                             frequency_matrices=nondirectional_matrix,
                                                             chromosome_name="genomic",
                                                             log_normalize=False,
                                                             plot=False)

    save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                             frequency_matrices=nondirectional_matrix,
                                                             chromosome_name="genomic",
                                                             log_normalize=True,
                                                             plot=False)

    # zero_mask = (matrix == 0)
    # nonzero_mask = numpy.invert(zero_mask)
    # matrix[zero_mask] += numpy.min(matrix[nonzero_mask])

    plot_directional_complementary_diff(matrix, max_runlength=MAX_RUNLENGTH)
    plot_base_matrices(matrix, test_spot=False, normalize_matrices=False)
    plot_base_matrices(matrix, test_spot=False, normalize_matrices=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pickle_dir", "-i",
        type=str,
        required=True,
        help="path to directory containing mah pickles"
    )
    parser.add_argument(
        "--output_dir", '-o',
        type=str,
        required=False,
        default="output",
        help="base directory for output"
    )
    args = parser.parse_args()

    main(pickle_dir=args.pickle_dir, output_parent_dir=args.output_dir)
