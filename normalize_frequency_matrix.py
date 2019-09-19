from modules.matrix import *
from handlers.FileManager import FileManager
import numpy


def main():
    matrix_path = "/home/ryan/code/runlength_analysis/output/runlength_matrix_shasta_human_chr1_GM24143/frequency_matrices_genomic_2019_4_23_21_38_57_361118.csv"
    output_dir = "output/runlength_matrix_shasta_human_chr1_GM24143/pseudocounts/"
    output_filename_prefix = "probability_matrices_GM24143_chr1_shasta_pseudocount_"

    FileManager.ensure_directory_exists(output_dir)

    matrix = load_base_length_matrix_from_csv(path=matrix_path, max_runlength=50)

    pseudocounts = [1, 4, 8, 16]

    for pseudocount in pseudocounts:
        filename = output_filename_prefix + str(pseudocount) + ".csv"
        save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                                 frequency_matrices=matrix,
                                                                 chromosome_name="genomic",
                                                                 log_normalize=True,
                                                                 filename=filename,
                                                                 pseudocount=pseudocount,
                                                                 plot=False)






if __name__ == "__main__":
    main()
