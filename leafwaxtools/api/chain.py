"""
The Chain module is the class for performing calculations using wax carbon 
chain-length concentration/abundance data imported as a 2D array-like object.
"""


import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from composition_stats import clr, closure, multiplicative_replacement
import scipy.stats
# import warnings
# from ..utils import validate_data


class Chain:
    """
    Represents leaf wax carbon chain-length concentration/abundance data
    imported as a 2-D, array-like object (i.e., list, array) with rows 
    representing unique samples and columns representing unique data types 
    (carbon chain-length number).
    
    Parameters
    ----------
    input_data : 2-D array-like
        User leaf wax chain-length concentration/abundance data.
        
    Attributes
    ----------
    input_data : 2-D array-like
        User leaf wax chain-length concentration/abundance data.
        
    
    Examples
    --------
    
    .. jupyter-execute::
        
        from leafwaxtools import Chain
        
    """
   
    def __init__(self, input_data):

        self.data = input_data
        
        if self.data.ndim != 2:
            raise TypeError("'input_data' must be 2-dimensional")


    def total_conc(self, zero_total=0, calculate_log=False):
        """
        Calculates the total concentration of each sample (rows).

        Parameters
        ----------
        zero_total : int, optional
            Return value if the sum of all columns in a row = 0. The default 
            is 0.
        calculate_log : bool, optional
            Returns log (base e) of the sum of each row instead of just the 
            sum. The default is False.

        Raises
        ------
        ValueError
            Raises an error when 'calculate_log' is neither True nor False.

        Returns
        -------
        total_conc : numpy.ndarray
            1-D Numpy array of total leaf wax concentrations for each sample 
            (row).

        """
        
        total_conc = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            total_conc[row] = np.nansum(self.data[row,:])

            if total_conc[row] == 0:
                total_conc[row] = zero_total

        match calculate_log:
            case True:
                total_conc = np.log(total_conc)
            case False:
                total_conc = total_conc
            case _:
                raise ValueError("'calculate_log' must either be True or False (default)")

        return total_conc


    def relative_abd(self, calculate_percent=False):
        """
        Calculates the relative abundance (fraction out of 1 or percentage) of 
        each leaf wax carbon chain-length (columns) for each sample (rows).

        Parameters
        ----------
        calculate_percent : bool, optional
            Calculate each chain-length relative abundance as a percentage 
            instead of a fraction of 1. The default is False.

        Raises
        ------
        ValueError
            Raises an error when 'calculate_percent' is neither True nor False.

        Returns
        -------
        rel_abd : numpy.ndarray
            2-D Numpy array of leaf wax chain-length relative abundances 
            (columns) for each sample (row).

        """

        rel_abd = np.zeros(np.shape(self.data))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):

                rel_abd[row,col] = self.data[row,col]/np.sum(self.data[row,:])
                
        match calculate_percent:
            case True:
                for row in range(0, len(self.data[:,0])):
                    for col in range(0, len(self.data[0,:])):
                        rel_abd[row,col] = rel_abd[row,col]*100
                
            case False:
                rel_abd = rel_abd
            
            case _:
                raise ValueError("'calculate_percent' must either be True or False (default)")

        return rel_abd


    def acl(self, chain_lengths):
        """
        Calculates the Average Chain-Length (ACL; Bray & Evans, 1961; Bush & 
        McInerney, 2013) of each sample (rows).
        
        References:
            
        Bray, E. E., & Evans, E. D. (1961). Distribution of n-paraffins as a 
        clue to recognition of source beds. Geochimica et Cosmochimica Acta, 
        22(1), 2-15. https://doi.org/10.1016/0016-7037(61)90069-2
        
        Bush, R. T., & McInerney, F. A. (2013). Leaf wax n-alkane 
        distributions in and across modern plants: Implications for 
        paleoecology and chemotaxonomy. Geochimica et Cosmochimica Acta, 117, 
        161-179. https://doi.org/10.1016/j.gca.2013.04.016

        Parameters
        ----------
        chain_lengths : list
            List of integers or floats representing the carbon chain-length 
            number of each column.

        Raises
        ------
        TypeError
            Raises an error if 'chain_lengths' is not a list.
        ValueError
            Raises an error if 'chain_lengths' is an empty list.

        Returns
        -------
        acl : numpy.ndarray
            1-D Numpy array of ACL values for each sample (row).

        """

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 24, 26, 28]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        acl_numer = np.zeros(len(self.data[:,0]))
        acl = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):

                acl_numer[row] += self.data[row,col] * chain_lengths[col]

            acl[row] = acl_numer[row]/np.sum(self.data[row,:])

        return acl


    def cpi(self, chain_lengths, even_over_odd=True):
        """
        Calculates the Carbon Preference Index (CPI; Marzi et al., 1993) of 
        each sample (rows).
        
        References:
            
        Marzi, R., Torkelson, B. E., & Olson, R. K. (1993). A revised carbon 
        preference index. Organic Geochemistry, 20(8), 1303-1306.
        https://doi.org/10.1016/0146-6380(93)90016-5

        Parameters
        ----------
        chain_lengths : list
            List of integers or floats representing the carbon chain-length 
            number of each column.
        even_over_odd : bool, optional
            Calculates the CPI of even-chain over odd-chain leaf waxes (use 
            case for n-alkanoic acids). Change to False to calculate the CPI 
            of odd-chain over even-chain waxes (use case for n-alkanes). The 
            default is True.

        Raises
        ------
        TypeError
            Raises an error if 'chain_lengths' is not a list.
        ValueError
            Raises an error if 'chain_lengths' is an empty list or if 
            'even_over_odd' is neither True nor False.

        Returns
        -------
        cpi : numpy.ndarray
            1-D Numpy array of CPI values for each sample (row).

        """

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 23, 24, 25, 26]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        '''
        EKT: use warnings to flag if even over odd order is wrong
        '''

        chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
        chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]

        data = pd.DataFrame(data=self.data, columns=(map(str, chain_lengths)))
        data_even = np.array(data.filter(items=(map(str, chain_lengths_even))))
        data_odd = np.array(data.filter(items=(map(str, chain_lengths_odd))))
        cpi = np.zeros(len(self.data[:,0]))

        match even_over_odd:
            case True:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_even[row,0:-1]) + np.nansum(data_even[row,1:])) / (2 * np.nansum(data_odd[row,:]))

            case False:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_odd[row,0:-1]) + np.nansum(data_odd[row,1:])) / (2 * np.nansum(data_even[row,:]))

            case _:
                raise ValueError("'even_over_odd' must be True (default) or False")

        return cpi


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



    def pca(self, chain_lengths, use_clr=True):
        """
        Performs a Principal Component Analysis (PCA) on the leaf wax 
        chain-length data with the centered log-ratio transform (clr; 
        Aitchison, 1982) applied to the input compositional data (Gloor et al. 
        2017).
                                                                  
        References:
            
        Aitchison, J. (1982). The statistical analysis of compositional data. 
        Journal of the Royal Statistical Society: Series B (Methodological), 
        44(2), 139-160. https://doi.org/10.1111/j.2517-6161.1982.tb01195.x
        
        Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. 
        (2017). Microbiome datasets are compositional: and this is not 
        optional. Frontiers in microbiology, 8, 2224.
        https://doi.org/10.3389/fmicb.2017.02224

        Parameters
        ----------
        chain_lengths : list
            List of integers or floats representing the carbon chain-length 
            number of each column..
        use_clr : bool, optional
            Calculates the clr of the leaf wax chain-length abundance data, 
            replacing 0 values with 1/N where N is the number of chain-lengths 
            (columns). The default is True.

        Raises
        ------
        TypeError
            Raises an error if 'chain_lengths' is not a list.
        ValueError
            Raises an error if 'chain_lengths' is an empty list or if 
            'use_clr' is neither True nor False.

        Returns
        -------
        pca_dict : dict
            Dictionary of PCA results including the PC scores for each factor 
            loading (chain-lengths/columns) and sample (rows). For more 
            details on all of the returns in pca_dict['pca'], please see the 
            documentation for the scikit-learn PCA function (sklearn.PCA()) 
            https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html

        """

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 23, 24, 25, 26]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        '''
        # deal with missing data
        for row in range(0, len(self.data[:,0]):
            for col i range(0, len(self.data[0,:]):

                if self.data[row,col] == np.nan:
                    self.data[row,col] = 0

            if np.sum(self.data[row,:]) == 0:
        '''


        match use_clr:
            case True:
                wax_relabd = closure(multiplicative_replacement(self.data))
                wax_clr = clr(wax_relabd)
                wax_data = pd.DataFrame(data=wax_clr, columns=chain_lengths)

            case False:
                wax_data = pd.DataFrame(data=self.relative_abd(), columns=chain_lengths)
                
            case _:
                raise ValueError("'use_clr' must be True or False (default)")

        wax_scaler = StandardScaler()
        wax_scaler.fit(wax_data)
        wax_data_scaled = wax_scaler.transform(wax_data)

        wax_pca = PCA(n_components=len(chain_lengths))
        wax_pca.fit_transform(wax_data_scaled)

        # wax_PC_scores = pd.DataFrame(
        #     wax_pca.fit_transform(wax_data_scaled),
        #     columns=chain_lengths
        # )
        # wax_loadings = pd.DataFrame(
        #     wax_pca.components_.T,
        #     columns=chain_lengths,
        #     index=wax_data.columns
        # )

        wax_ldings = wax_pca.components_
        wax_features = wax_data.columns
        wax_pc_values = np.arange(wax_pca.n_components_) + 1

        pca_dict = {
            "pca": wax_pca,
            "pc_values": wax_pc_values,
            "features": wax_features,
            "loadings": wax_ldings
        }

        for i in range(0, len(chain_lengths)):

            wax_pc = wax_pca.fit_transform(wax_data_scaled)[:,i]
            wax_scale_pc = 1.0 / (wax_pc.max() - wax_pc.min())
            wax_pc_score = wax_pc * wax_scale_pc

            # pca_dict.update({f"wax_pc{i+1}": wax_pc})
            # pca_dict.update({f"wax_scale_pc{i+1}": wax_scale_pc})
            pca_dict.update({f"wax_pc{i+1}_score": wax_pc_score})


        return pca_dict
