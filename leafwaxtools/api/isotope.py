"""
The Iso module is the class for performing calculation using plant wax stable
isotope data imported as a 2D array-like object
"""

import pandas as pd
import numpy as np
from ..utils import validate_init
from ..utils import correlation


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
        
        validate_init.validate_data_dimensions(self.data)


    def value_range(self):
        """
        Calculates the maximum isotope value range (max - min) between all 
        chain-lengths (columns) of each sample (rows).

        Returns
        -------
        value_range : numpy.ndarray
            1-D Numpy array of maximum isotope value ranges for each sample 
            (row).

        """

        value_range = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            value_range[row] = np.nanmax(self.data[row,:]) - np.nanmin(self.data[row,:])
        
        return value_range


    def concentration_avg(self, chain_data):
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
        concentration_avg : numpy.ndarray
            1-D Numpy array of chain-length concentration-weighted average 
            isotope values for each sample (row).

        """
        
        if np.shape(self.data) != np.shape(chain_data):
            raise ValueError("Input isotope and chain-length distribution data must have the same number of rows and columns")

        concentration_avg = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                if self.data[row,col] == np.nan:
                    chain_data[row,col] = 0

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                concentration_avg[row] += self.data[row,col] * chain_data[row,col]

            concentration_avg[row] = concentration_avg[row]/np.sum(chain_data[row,:])

        return concentration_avg

    
    def epsilon(self, epsilon_numerator=None, epsilon_denominator=None):
        """
        Calculates the isotopic fractionation factor (epsilon) between stable 
        isotope values/arrays (permil units) in the numerator and denominator 
        of the following equation (e.g., Diefendorf & Freimuth, 2017):
            
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
    
    
    def wax_to_source(self, epsilon, epsilon_numerator=None):
        """
        Calculates the isotopic value of a source material to a leaf wax using 
        an isotopic fractionation factor (epsilon). A common application is the
        calculation of source water (e.g., precipitation, lake water) stable
        hydrogen isotope (d2H) values using leaf wax d2H values and an
        associated epsilon value between the source water and leaf wax (e.g., 
        Feakins, 2013; Holtzman et al., 2025).
        
        References:
            
        Feakins, S. J. (2013). Pollen-corrected leaf wax D/H reconstructions of
        northeast African hydrological changes during the late Miocene. 
        Palaeogeography, Palaeoclimatology, Palaeoecology, 374, 62-71.
        https://doi.org/10.1016/j.palaeo.2013.01.004
        
        Holtzman, H., Thomas, E. K., Erb, M., Marshall, L., Castañeda, I. S., 
        Kaufman, D., ... & Melles, M. (2025). Early Holocene atmospheric 
        circulation changes over northern Europe based on isotopic and 
        biomarker evidence from Kola Peninsula. Paleoceanography and 
        Paleoclimatology, 40(3), e2024PA005076. 
        https://doi.org/10.1029/2024PA005076

        Parameters
        ----------
        epsilon : 1-D or 2-D array-like
            Isotopic fractionation factor (epsilon) value/array.
        epsilon_numerator : 1-D or 2-D array-like, optional
            Numerator stable isotope value/array. Uses Isotope.data by default 
            if no argument is passed. The default is None.

        Returns
        -------
        source_isotope : 1-D or 2-D array-like
            1-D or 2-D array-like of source isotope values for each sample 
            (row) and chain-length (column; if applicable).

        """        
        
        if epsilon_numerator is None:
            epsilon_numerator = self.data
            
        source_isotope = ((1000+epsilon_numerator)/((epsilon/1000)+1))-1000
        
        return source_isotope


    def correlation_rvals(self, minimum_obs=2):
        """
        Calculates the Pearson correlation r-values between each leaf wax 
        chain-length (columns). To be extended with other correlation methods 
        (Spearman, Kendall Tau) in a future version. This functionality is 
        identical between the Chain and Isotope API classes.

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

        r_vals = correlation.corr_r(data=self.data, min_obs=minimum_obs)
        
        return r_vals


    def correlation_pvals(self, minimum_obs=2):
        """
        Calculates the Pearson correlation p-values between each leaf wax 
        chain-length (columns). To be extended with other correlation methods 
        (Spearman, Kendall Tau) in a future version. This functionality is 
        identical between the Chain and Isotope API classes.

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

        p_vals = correlation.corr_p(data=self.data, min_obs=minimum_obs)
                
        return p_vals
