"""
The Iso module is the class for performing calculation using plant wax stable
isotope data imported as a 2D array-like object
"""

from enum import Enum
import pandas as pd
import numpy as np
import warnings
from ..utils import validate_data
from ..utils import data_type_enum


class Iso:


    def __init__(self, input_data):

        self.data = input_data


    def iso_range():

        iso_range = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            iso_range[row] = np.max(self.data[row,:]) - np.min(self.data[row,:])

        return iso_range


    def iso_avg(chain_data):

        if np.shape(self.data) != np.shape(chain_data):
            raise ValueError("Input isotope and chain-length distribution data have the same number of rows and columns")

        iso_avg = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                if self.data[row,col] == np.nan:
                    chain_data[row,col] = 0

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                iso_avg[row] += self.data[row,col] * chain_data[row,col]

            iso_avg[row] = iso_avg[row]/np.sum(chain_data[row,:])

        return iso_avg


    def epsilon(epsilon_numer, epsilon_denom):

        epsilon = np.zeros(len(epsilon_numer[:,0]))

        for row in range(0, len(epsilon_numer[:,0])):
            epsilon[row] = (((1000+epsilon_numer[row])/(1000+epsilon_denom[row]))-1)*1000

        return epsilon
