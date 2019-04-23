from modules.matrix import *
import numpy


def main():
    matrix_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_shasta_human_chr1_GM24143/frequency_matrices_genomic_2019_4_23_21_38_57_361118.csv"

    matrix = load_base_length_matrix_from_csv(path=matrix_path, max_runlength=50)

    plot_base_matrices(matrix=matrix, cutoff=50)


if __name__ == "__main__":
    main()