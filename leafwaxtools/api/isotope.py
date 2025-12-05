"""
The Iso module is the class for performing calculation using plant wax stable
isotope data imported as a 2D array-like object
"""

# import pandas as pd
import numpy as np
import scipy.stats
# from ..utils import validate_data


class Isotope:
    """
    Represents leaf wax compound-specific stable isotope data imported as a 
    2-D, array-like object (i.e., list, array) with rows representing unique 
    samples and columns representing unique data types (carbon chain-length 
    number).
    
    Parameters
    ----------
    input_data : 2-D array-like
        User leaf wax isotope data.
        
    Attributes
    ----------
    input_data : 2-D array-like
        User leaf wax isotope data.
        
    
    Examples
    --------
    
    .. jupyter-execute::
        
        from leafwaxtools import Isotope
        
    """

    def __init__(self, input_data):

        self.data = input_data
        
        if self.data.ndim != 2:
            raise TypeError("'input_data' must be 2-dimensional")


    def value_range(self):
        """
        Calculates the maximum isotope value range (max. - min.) between all 
        chain-lengths (columns) of each sample (rows).

        Returns
        -------
        value_range : numpy.ndarray
            1-D Numpy array of maximum isotope value ranges for each sample 
            (row).

        """

        value_range = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            value_range[row] = np.max(self.data[row,:]) - np.min(self.data[row,:])

        return value_range


    def conc_avg(self, chain_data):
        """
        Calculates the chain-length concentration-weighted average isotope 
        value of each sample (rows) using .

        Parameters
        ----------
        chain_data : 2-D array-like
            User leaf wax chain-length concentration/abundance data. Must be 
            the same shape (same number of rows and columns) as the user 
            isotope input data. Chain-length data with a value of NaN is 
            treated as having a concentration/abundance of 0.

        Raises
        ------
        ValueError
            Raises an error when Isotope.data and 'chain_data' are not the 
            same shape (same number of rows and columns).

        Returns
        -------
        conc_avg : numpy.ndarray
            1-D Numpy array of chain-length concentration-weighted average 
            isotope values for each sample (row).

        """
        
        if np.shape(self.data) != np.shape(chain_data):
            raise ValueError("Input isotope and chain-length distribution data must have the same number of rows and columns")

        conc_avg = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                if self.data[row,col] == np.nan:
                    chain_data[row,col] = 0

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                conc_avg[row] += self.data[row,col] * chain_data[row,col]

            conc_avg[row] = conc_avg[row]/np.sum(chain_data[row,:])

        return conc_avg

    
    def epsilon(self, epsilon_numerator=None, epsilon_denominator=None):
        """
        Calculates the isotopic fractionation factor (epsilon) between stable 
        isotope values/arrays (permil units) in the numerator and denominator 
        of the following equation (e.g. Diefendorf & Freimuth, 2017):
            
        epsilon = (((1000 + numerator) / (1000 + denominator)) - 1) * 1000
        
        This equation is based on epsilon calculations for leaf wax stable 
        hydrogen and carbon isotope vlaues.
        
        References:
            
        Diefendorf, A. F., & Freimuth, E. J. (2017). Extracting the most from 
        terrestrial plant-derived n-alkyl lipids and their carbon isotopes 
        from the sedimentary record: A review. Organic Geochemistry, 103, 
        1-21. https://doi.org/10.1016/j.orggeochem.2016.10.016

        Parameters
        ----------
        epsilon_numerator : 1-D or 2-D array-like, optional
            Numerator stable isotope value/array. Uses Isotope.data by default 
            if no argument is passed. The default is None.
        epsilon_denominator : 1-D or 2-D array-like, optional
            Denominator stable isotope value/array. Uses Isotope.data by 
            default if no argument is passed. The default is None.

        Returns
        -------
        epsilon : 1-D or 2-D array-like
            1-D or 2-D array-like of epsilon values for each sample (row) and 
            chain-length (column; if applicable).

        """
        
        if epsilon_numerator is None:
            epsilon_numerator = self.data

        if epsilon_denominator is None:
            epsilon_denominator = self.data

        epsilon = (((1000+epsilon_numerator)/(1000+epsilon_denominator))-1)*1000

        return epsilon
    

    def corr_rvals(self, minimum_obs=2):
        """
        Calculates the Pearson correlation r-values between each leaf wax 
        chain-length (columns). To be extended with other correlation methods 
        (Spearman, Kendall Tau) in a future version.

        Parameters
        ----------
        minimum_obs : int, optional
            Minimum number of observations (samples/rows) required to return a
            Pearson r-value. The default is 2.

        Returns
        -------
        r_vals : numpy.ndarray
            2-D Numpy array of Pearson correlation r-values between each leaf 
            wax chain-length (column) with all values in the major diagonal 
            equal to 1.

        """

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
        """
        Calculates the Pearson correlation p-values between each leaf wax 
        chain-length (columns). To be extended with other correlation methods 
        (Spearman, Kendall Tau) in a future version.

        Parameters
        ----------
        minimum_obs : int, optional
           Minimum number of observations (samples/rows) required to return a
           Pearson r-value. The default is 2.

        Returns
        -------
        p_vals : numpy.ndarray
            2-D Numpy array of Pearson correlation p-values between each leaf 
            wax chain-length (column).

        """

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
