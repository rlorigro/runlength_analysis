import numpy


def main():
    coverage = 10               # How many bases to sample per trial
    n_samples = 1000000         # This should be much greater than the number of significant digits in target quality
    p_error = 0.07              # Quality
    p_success = 1 - p_error     # For binomial parameter

    # a vector of size n_samples, where each value is the rate of success (proportion of accurate bases)
    trials = numpy.random.binomial(coverage, p_success, n_samples)/coverage

    # Boolean mask. Binarizes the rate of positives for each trial based on whether they are > the cutoff
    success_mask = (trials > 0.5)

    # Count the proportion of consensuses that would have given the correct base
    p_success_consensus = numpy.sum(success_mask)/n_samples
    print(p_success_consensus)


if __name__ == "__main__":
    main()
