from matplotlib import pyplot
import numpy


def weibull_cdf(x, shape, scale):
    """  Tail cumulative probability of Weibull distribution

        :param x: points at which to evaluate
        :param shape: Shape parameter
        :param scale: Scale parameter
    """
    return numpy.exp(-numpy.power(x / scale, shape))


def evaluate_discrete_weibull(scale, shape, x):
    """
    Evaluate the discrete form of the distribution at points specified by x vector
    :param scale:
    :param shape:
    :param x:
    :return:
    """
    distribution = weibull_cdf(x, shape, scale) - weibull_cdf(x + 1, shape, scale)

    return distribution


def calculate_mode(scale, shape):
    """
    From wiki, lambda=scale, k=shape
    :param scale:
    :param shape:
    :param x:
    :return:
    """

    if shape > 1:
        mode = scale * ((shape-1) / shape)**(1 / shape)
        mode = max(0, mode)

    else:
        mode = 0

    return mode


# def plot_distribution(axes, x, y):
#     color = [0,0.05,0.95]
#
#     axes.bar(x, y, color=color, alpha=0.2)


def main():
    pass


if __name__ == "__main__":
    main()
