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
        interquantile = np.percentile(numpy_array, 75, interpolation='higher') - np.percentile(
            numpy_array, 25, interpolation='lower')
        # Note: numpy func 'to_list' convert a numpy array into a normal python data type
        return interquantile.to_list()

    @staticmethod
    def median(array):
        """
        Returns the median value of a list.
        :param array: data array
        :return: Median value
        """
        return np.median(array).to_list()

    @staticmethod
    def cov(array):
        """
        Calculate coefficient of variation of the data
        :param array: data array
        :return: coefficient of variation of array
        """
        return scipy.stats.variation(array).to_list()

    @staticmethod
    def std(array):
        """
        Calculate standard deviation of the data
        :param array: data array
        :return: standard deviation of array
        """
        return np.std(np.array(array)).to_list()
