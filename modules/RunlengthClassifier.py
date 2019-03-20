import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
from measure_runlength_modal_distribution_from_runnie_sequence import MAX_RUNLENGTH
from handlers.FileManager import FileManager
from matplotlib import pyplot
import numpy


numpy.set_printoptions(precision=9, linewidth=sys.maxsize, suppress=True, threshold=sys.maxsize)

MATRIX_PATH = "/home/ryan/code/runnie_parser/output/runlength_matrix_from_sequence_runnie_vMode2/probability_matrices_2019_3_12_13_52_19_93991.csv"    # E coli stranded guppy wg first50k reads

# Indexing RLE tuples
SEQUENCE = 0
LENGTHS = 1

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


class RunlengthClassifier:
    """
    Calculate the probability of true runlength given the following:
        1) a vector of observed runlengths for a set of aligned reads
        2) a matrix of dimensions [true_runlength, observed_runlength] containing raw true/observed frequencies
    """

    def __init__(self, path, log_scale=True):
        self.log_scale = log_scale

        self.pseudocount = 15

        self.probability_matrices = self.load_base_probability_matrix_from_csv(path)

        # self.y_maxes = [matrix.shape[0] for matrix in self.probability_matrices[0]]
        # self.x_maxes = [matrix.shape[1] for matrix in self.probability_matrices[0]]

        # self.export_matrices_to_fasta_one_liner()

    def plot_matrix(self, probability_matrix, title=""):
        axes = pyplot.axes()
        pyplot.imshow(probability_matrix)
        axes.set_xlabel("Observed Runlength")
        axes.set_ylabel("True Runlength")
        pyplot.title(title)
        pyplot.show()
        pyplot.close()

    def load_base_probability_matrix_from_csv(self, path):
        matrices = numpy.zeros([4, MAX_RUNLENGTH + 1, MAX_RUNLENGTH + 1])
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

    def log_sum_exp(self, x):
        """
        Non-log addition in log-space vector of values... doesn't work for signed addition? currently unused
        :param x:
        :return:
        """
        b = numpy.max(x[(x<sys.maxsize)])   # ignore inf values

        s = b + numpy.log(numpy.sum(numpy.exp(x-b)))

        return s

    def normalize_likelihoods(self, log_likelihood_y, max_index):
        """
        Given a vector of log likelihood values for each Y, and the index of the maximum p(y), normalize each value wrt
        the maximum: p(Y_i|x)/p(Y_max|x)
        :param log_likelihood_y:
        :param max_index:
        :return:
        """
        max_value = log_likelihood_y[max_index,:]

        normalized_likelihoods = log_likelihood_y - max_value

        return normalized_likelihoods

    def factor_repeats(self, x, reversal):
        """
        Given a vector of repeats factor them into counts for each unique repeat length
        e.g. [1,1,1,1,2,2,2] -> [1,2], [4,3]
        :param x:
        :return:
        """
        # print("x\t\t", x)
        # print("reversal\t", reversal)
        reversal = reversal.squeeze()
        x = x.squeeze()

        x_forward = x[numpy.invert(reversal)]
        x_reverse = x[reversal]

        # print(x_forward)
        # print(x_reverse)

        unique_forward, inverse_forward = numpy.unique(x_forward, return_inverse=True)
        unique_reverse, inverse_reverse = numpy.unique(x_reverse, return_inverse=True)

        # print(unique_forward)
        # print(inverse_forward)
        # print(unique_reverse)
        # print(inverse_reverse)

        bincount_forward = numpy.bincount(inverse_forward)
        bincount_reverse = numpy.bincount(inverse_reverse)

        unique = numpy.concatenate([unique_forward, unique_reverse])
        bincount = numpy.concatenate([bincount_forward, bincount_reverse])

        n_forward = unique_forward.shape[0]
        n_reverse = unique_reverse.shape[0]

        reversal = numpy.concatenate([numpy.zeros(n_forward, dtype=bool), numpy.ones(n_reverse, dtype=bool)])

        return unique, bincount, reversal

    def predict(self, character_index, x, reversal, skip_zeros=True, randomize_reversal=False):
        """
        for a vector of observations x, find the product of likelihoods p(x_i|y_j) for x_i in x, for all possible Y
        values, and return the maximum value
        :param x: array/vector of repeat values
        :param reversal: array/vector of binary reversal statuses (1=reversed)
        :param character_index: the character (base) of the consensus of this column of aligned reads (should be on
        per-base basis?)
        :return:
        """
        if randomize_reversal:
            l = reversal.size()
            print(l)
            print(reversal)
            reversal = numpy.random.binomial(1, 0.5, l).astype(numpy.bool)
            print(reversal)

        # factor the repeats to avoid iterating probability lookups multiple times
        x, counts, reversal = self.factor_repeats(x, reversal)

        log_likelihood_y = numpy.zeros([MAX_RUNLENGTH+1, 1])

        for y_j in range(0, MAX_RUNLENGTH + 1):
            # initialize log likelihood for this (jth) y value
            log_sum = 0

            for i in range(x.shape[0]):
                # iterate x values, summing log likelihood for every p(x_i|y_j)
                x_i = x[i]
                c_i = counts[i]
                r_i = reversal[i]

                # convert to indices
                x_i = int(x_i)
                y_j = int(y_j)

                if skip_zeros:
                    if x_i == 0:
                        continue

                if x_i > MAX_RUNLENGTH:
                    x_i = MAX_RUNLENGTH

                if r_i:
                    character_index = 3 - character_index

                # retrieve conditional probability for this x|y... index assumes zeros (row+col) are included matrix
                prob_x_i_given_y_j = self.probability_matrices[character_index, y_j, x_i]

                # print(r_i, x_i, y_j, prob_x_i_given_y_j, 10**prob_x_i_given_y_j)

                # exponentiate by the number of independently observed repeats of this value
                log_sum += c_i*float(prob_x_i_given_y_j)

            if numpy.isnan(log_sum):    # if the frequency matrix has empty rows, consider p(y_j) to be 0
                log_sum = -numpy.inf

            # store result of log sum of likelihoods
            log_likelihood_y[y_j,0] = log_sum

        j_max = numpy.argmax(log_likelihood_y)

        normalized_posterior = self.normalize_likelihoods(log_likelihood_y=log_likelihood_y, max_index=j_max)

        # print(10**normalized_posterior)
        # print(j_max)

        return normalized_posterior, j_max

    def print_normalized_likelihoods(self, normalized_likelihoods):
        for i in range(normalized_likelihoods.shape[0]):
            print("%d:\t%.3f" % (i,float(normalized_likelihoods[i,:])))


def test():
    runlength_classifier = RunlengthClassifier(MATRIX_PATH)

    while True:
        sys.stdout.write("Enter character: \n")
        input_string = sys.stdin.readline().strip().upper()
        character_index = BASE_TO_INDEX[input_string]

        sys.stdout.write("Enter space-separated repeat observations: \n")
        input_string = sys.stdin.readline()
        x = input_string.strip().split(" ")
        x = numpy.array(list(map(int, x)))

        sys.stdout.write("Enter space-separated reversal statuses: \n")
        input_string = sys.stdin.readline()
        reversal = input_string.strip().split(" ")
        reversal = numpy.array(list(map(int, reversal)), dtype=numpy.bool)

        print(x, character_index, reversal)

        normalized_y_log_likelihoods, y_max = runlength_classifier.predict(x=x,
                                                                           character_index=character_index,
                                                                           reversal=reversal)

        runlength_classifier.print_normalized_likelihoods(10**normalized_y_log_likelihoods)

        print("\nMost likely runlength: ", y_max)


if __name__ == "__main__":
    test()
