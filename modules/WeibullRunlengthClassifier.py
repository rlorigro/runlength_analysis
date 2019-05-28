import sys
import os
sys.path.append(os.path.dirname(sys.path[0]))
from measure_runlength_modal_distribution_from_runnie import MAX_RUNLENGTH
from discrete_weibull_distribution import evaluate_discrete_weibull
from modules.matrix import normalize
from matplotlib import pyplot
import numpy
import math


numpy.set_printoptions(precision=9, linewidth=sys.maxsize, suppress=True, threshold=sys.maxsize)

MATRIX_PATH = "/home/ryan/code/runnie_parser/output/runlength_matrix_from_runnie_sequence_0_1_2_3_TRAIN_ecoli/probability_matrices_2019_3_19_15_33_44_765752.csv"    # E coli stranded guppy wg first50k reads

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


def plot_matrices(matrix, cutoff=12):
    figure, axes = pyplot.subplots(nrows=2, ncols=2)
    figure.set_size_inches(8, 8)

    matrix_A = matrix[0, :cutoff, :cutoff]
    axes[0][0].imshow(matrix_A)
    axes[0][0].set_title(INDEX_TO_BASE[0])

    matrix_T = matrix[3, :cutoff, :cutoff]

    axes[1][0].imshow(matrix_T)
    axes[1][0].set_title(INDEX_TO_BASE[3])

    matrix_G = matrix[2, :cutoff, :cutoff]

    axes[0][1].imshow(matrix_G)
    axes[0][1].set_title(INDEX_TO_BASE[2])

    matrix_C = matrix[1, :cutoff, :cutoff]

    axes[1][1].imshow(matrix_C)
    axes[1][1].set_title(INDEX_TO_BASE[1])

    axes[1][1].set_xlabel("Observed length")
    axes[1][0].set_xlabel("Observed length")
    axes[1][0].set_ylabel("True length")
    axes[0][0].set_ylabel("True length")

    axes[1][1].set_yticks(numpy.arange(0, cutoff, 2))
    axes[1][0].set_yticks(numpy.arange(0, cutoff, 2))
    axes[0][0].set_yticks(numpy.arange(0, cutoff, 2))
    axes[0][1].set_yticks(numpy.arange(0, cutoff, 2))

    axes[1][1].set_xticks(numpy.arange(0, cutoff, 2))
    axes[1][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][0].set_xticks(numpy.arange(0, cutoff, 2))
    axes[0][1].set_xticks(numpy.arange(0, cutoff, 2))


class RunlengthClassifier:
    """
    Calculate the probability of true runlength given the following:
        1) a vector of observed runlengths for a set of aligned reads
        2) a matrix of dimensions [true_runlength, observed_runlength] containing raw true/observed frequencies
    """

    def __init__(self, path, normalize_matrix=True, pseudocount=1):
        self.pseudocount = pseudocount

        self.probability_matrices = self.load_base_probability_matrix_from_csv(path)

        if normalize_matrix:
            print(self.probability_matrices.shape)
            for i in range(4):
                self.probability_matrices[i,:,:] = normalize(self.probability_matrices[i,:,:], pseudocount=self.pseudocount)

        plot_matrices(self.probability_matrices, cutoff=20)

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
        matrices = numpy.zeros([4, MAX_RUNLENGTH+1, MAX_RUNLENGTH+1], dtype=numpy.float)
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

        # pyplot.imshow(matrices[0,:,:])
        # pyplot.show()
        # pyplot.close()

        # trim zeros :(
        matrices = matrices[:, 1:, 1:]

        # pyplot.imshow(matrices[0,:,:])
        # pyplot.show()
        # pyplot.close()

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

    def plot_inference(self, observations, cutoff=None):
        if cutoff is None:
            n_observations = len(observations)
        else:
            n_observations = min(len(observations), cutoff)

        figure, axes = pyplot.subplots(nrows=n_observations)

        for i in range(n_observations):
            scale, shape, distribution = observations[i]
            print(scale, shape, distribution[:8])
            axes[i].imshow(numpy.atleast_2d(distribution))

        pyplot.show()
        pyplot.close()

        return

    def predict(self, character_index, x_scales, x_shapes, reversal, randomize_reversal=False, prior=False):
        """
        for a vector of observations x, find the product of likelihoods p(x_i|y_j) for x_i in x, for all possible Y
        values, and return the maximum value
        :param x: array/vector of repeat values
        :param reversal: array/vector of binary reversal statuses (1=reversed)
        :param character_index: the character (base) of the consensus of this column of aligned reads (should be on
        per-base basis?)
        :return:
        """
        x_range = numpy.arange(0,MAX_RUNLENGTH)

        if randomize_reversal:
            l = reversal.size()
            reversal = numpy.random.binomial(1, 0.5, l).astype(numpy.bool)

        log_likelihood_y = numpy.zeros([MAX_RUNLENGTH, 1])

        # Precompute full distributions for all params in input
        x_weibull_distributions = list()
        for i in range(x_scales.shape[0]):
            if x_scales[i] < 0 and x_shapes[i] < 0:
                x_i_weibull = None
            else:
                # each observation contains a distribution over x
                x_i_weibull = evaluate_discrete_weibull(scale=x_scales[i], shape=x_shapes[i], x=x_range)

            x_weibull_distributions.append(x_i_weibull)

        # Use weibull distributions for inference
        for y_j in range(0, MAX_RUNLENGTH):
            # initialize log likelihood for this (jth) y value, use prior if specified
            if prior:
                log_sum = math.log(3.5, 10)*y_j if y_j > 0 else -1e9
            else:
                log_sum = 0

            # observations = list()
            for i in range(x_scales.shape[0]):
                if x_scales[i] < 0 and x_shapes[i] < 0:
                    continue

                # each observation contains a distribution over x
                x_i_weibull = x_weibull_distributions[i]

                # observations.append([x_scales[i], x_shapes[i], x_i_weibull])

                x_i_weibull = numpy.log10(x_i_weibull)

                # iterate x values, summing log likelihood for every p(x_i|y_j)
                r_i = reversal[i]

                # convert to indices
                y_j = int(y_j)

                # complement character if the pileup contains a reversed strand
                if r_i:
                    character_index = 3 - character_index

                # retrieve conditional probability for this x|y
                prob_x_i_given_y_j = self.probability_matrices[character_index, y_j, :] + x_i_weibull

                # print(r_i, x_i, y_j, prob_x_i_given_y_j, 10**prob_x_i_given_y_j)

                # update the log likelihood for this Y value, using the product of the observed weibull and the model
                # sum over all possible (x,y) products for the distribution of x
                log_sum += float(self.log_sum_exp(prob_x_i_given_y_j))

            if numpy.isnan(log_sum):    # if the frequency matrix has empty rows, consider p(y_j) to be 0
                log_sum = -numpy.inf

            # store result of log sum of likelihoods
            log_likelihood_y[y_j,0] = log_sum

        # self.plot_inference(observations, cutoff=10)

        j_max = numpy.argmax(log_likelihood_y) + 1  # weibull runlength is 0 based, should be 1 based

        normalized_posterior = self.normalize_likelihoods(log_likelihood_y=log_likelihood_y, max_index=j_max)

        # print(10**normalized_posterior)
        # print(j_max)

        return normalized_posterior, j_max

    def print_normalized_likelihoods(self, normalized_likelihoods):
        for i in range(normalized_likelihoods.shape[0]):
            print("%d:\t%.3f" % (i,float(normalized_likelihoods[i,:])))


if __name__ == "__main__":
    pass
