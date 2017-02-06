import logging
import numpy as np
import scipy.stats


class StatisticsUtils(object):
    """
    Class providing statistic utility functions
    """

    @staticmethod
    def mad(array, median):
        """
        Calculate median absolute deviation (MAD)
        :param array: data array
        :param median: the median value of array
        :return: MAD
        """
        return StatisticsUtils.median(np.absolute([x - median for x in array]))

    @staticmethod
    def interquantile(array):
        """
        Calculate interquantile of the data
        :param array: data array
        :return: interquantile of the data (25%~75%)
        """
        numpy_array = np.array(array)
        return np.percentile(numpy_array, 75, interpolation='higher') - np.percentile(numpy_array, 25, interpolation='lower')

    @staticmethod
    def median(array):
        """
        Returns the median value of a list.
        :param array: data array
        :return: Median value
        """
        sorted_list = sorted(array)
        middle = len(array) // 2
        logging.debug("Median: len {} middle {} sorted_len {}".format(len(array), middle, len(sorted_list)))
        if len(array) % 2:
            return sorted_list[middle]
        else:
            return (sorted_list[middle] + sorted_list[middle - 1]) / 2

    @staticmethod
    def cov(array):
        """
        Calculate coefficient of variation of the data
        :param array: data array
        :return: coefficient of variation of array
        """
        return scipy.stats.variation(array)

    @staticmethod
    def std(array):
        """
        Calculate standard deviation of the data
        :param array: data array
        :return: standard deviation of array
        """
        return np.std(np.array(array))
