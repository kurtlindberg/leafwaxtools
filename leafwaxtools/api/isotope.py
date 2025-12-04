"""
The Iso module is the class for performing calculation using plant wax stable
isotope data imported as a 2D array-like object
"""

# import pandas as pd
import numpy as np
import scipy.stats
# from ..utils import validate_data


class Isotope:


    def __init__(self, input_data):

        self.data = input_data


    def iso_range(self):

        iso_range = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            iso_range[row] = np.max(self.data[row,:]) - np.min(self.data[row,:])

        return iso_range


    def iso_avg(self, chain_data):

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

    
    def epsilon(self, epsilon_numerator=None, epsilon_denominator=None):

        if epsilon_numerator is None:
            epsilon_numerator = self.data

        if epsilon_denominator is None:
            epsilon_denominator = self.data

        epsilon = (((1000+epsilon_numerator)/(1000+epsilon_denominator))-1)*1000

        return epsilon
    

    def corr_rvals(self, minimum_obs=2):

        r_vals = np.zeros((len(self.data[0,:]), len(self.data[0,:])))

        for row in range(0, len(r_vals[:,0])):
            for col in range(0, len(r_vals[0,:])):

                x_corr = np.array(self.data[:,row])
                y_corr = np.array(self.data[:,col])

                if (len(x_corr) >= minimum_obs) and (len(y_corr) >= minimum_obs):
                    r_vals[row,col] = scipy.stats.pearsonr(x_corr, y_corr)[0]
                else:
                    r_vals[row,col] = np.nan

        return r_vals


    def corr_pvals(self, minimum_obs=2):

        p_vals = np.zeros((len(self.data[0,:]), len(self.data[0,:])))

        for row in range(0, len(p_vals[:,0])):
            for col in range(0, len(p_vals[0,:])):

                x_corr = np.array(self.data[:,row])
                y_corr = np.array(self.data[:,col])

                if (len(x_corr) >= minimum_obs) and (len(y_corr) >= minimum_obs):
                    p_vals[row,col] = scipy.stats.pearsonr(x_corr, y_corr)[1]
                else:
                    p_vals[row,col] = np.nan

        return p_vals
