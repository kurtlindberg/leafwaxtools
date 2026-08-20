"""
The Chain module is the class for performing calculations using wax carbon 
chain-length concentration/abundance data imported as a 2D array-like object.
"""


import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from skbio.stats.composition import clr, multi_replace
import warnings
from ..utils import validate_init
from ..utils import correlation


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

        validate_init.validate_data_dimensions(input_data)
        
        input_data_df = pd.DataFrame(data=input_data)
        input_data_dfnan = input_data_df.apply(pd.to_numeric, errors='coerce')
        input_data_arr = np.array(input_data_dfnan)

        self.data = input_data_arr


    def total_conc(self, calculate_log=False):
        """
        Calculates the total concentration of each sample (rows). This function
        utilizes numpy.nansum(), which ignores all NaN values in each sample
        and returns 0 for a row of all NaNs.

        Parameters
        ----------
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

        if calculate_log is True:
            total_conc = np.log(total_conc)

        elif calculate_log is False:
            total_conc = total_conc

        else:
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
                rel_abd[row,col] = self.data[row,col]/np.nansum(self.data[row,:])
                    
        if calculate_percent is True:
            for row in range(0, len(self.data[:,0])):
                for col in range(0, len(self.data[0,:])):
                    rel_abd[row,col] = rel_abd[row,col]*100
                
        elif calculate_percent is False:
            rel_abd = rel_abd

        else:
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
        chain_lengths : array-like
            Array-like of integers or floats representing the carbon 
            chain-length number of each column.

        Raises
        ------
        ValueError
            Raises an error if 'chain_lengths' is not the same length as the
            number of chain-lengths (columns).

        Returns
        -------
        acl : numpy.ndarray
            1-D Numpy array of ACL values for each sample (row).

        """

        if len(chain_lengths) != len(self.data[0,:]):
            raise ValueError(
                "'chain_lengths' must be the same length as the number of data columns"    
            )

        acl_numer = np.zeros(len(self.data[:,0]))
        acl = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):
                acl_numer[row] += self.data[row,col] * chain_lengths[col]

            acl[row] = acl_numer[row]/np.nansum(self.data[row,:])

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
        chain_lengths : array-like
            Array-like of integers or floats representing the carbon 
            chain-length number of each column.
        even_over_odd : bool, optional
            Calculates the CPI of even-chain over odd-chain leaf waxes (use 
            case for n-alkanoic acids). Change to False to calculate the CPI 
            of odd-chain over even-chain waxes (use case for n-alkanes). The 
            default is True.

        Raises
        ------
        ValueError
            Raises an error if 'chain_lengths' is not the same length as the
            number of chain-lengths (columns) or if 'even_over_odd' is neither 
            True nor False.

        Returns
        -------
        cpi : numpy.ndarray
            1-D Numpy array of CPI values for each sample (row).

        """

        if len(chain_lengths) != len(self.data[0,:]):
            raise ValueError(
                "'chain_lengths' must be the same length as the number of data columns"    
            )

        if chain_lengths[0] % 2 != 0 and even_over_odd is True:
            warnings.warn(
                f"The first chain-length '{chain_lengths[0]}' is an odd number, but this cpi function will be dividing even number chain-lengths over odd number ones"
            )
        
        if chain_lengths[0] % 2 == 0 and even_over_odd is False:
            warnings.warn(
                f"The first chain-length '{chain_lengths[0]}' is an even number, but this cpi function will be dividing odd number chain-lengths over even number ones"
            )

        chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
        chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]

        data = pd.DataFrame(data=self.data, columns=(map(str, chain_lengths)))
        data_even = np.array(data.filter(items=(map(str, chain_lengths_even))))
        data_odd = np.array(data.filter(items=(map(str, chain_lengths_odd))))
        cpi = np.zeros(len(self.data[:,0]))

        if even_over_odd is True:
            for row in range(0, len(self.data[:,0])):
                cpi[row] = (np.nansum(data_even[row,0:-1]) + np.nansum(data_even[row,1:])) / (2 * np.nansum(data_odd[row,:]))

        elif even_over_odd is False:
             for row in range(0, len(self.data[:,0])):
                 cpi[row] = (np.nansum(data_odd[row,0:-1]) + np.nansum(data_odd[row,1:])) / (2 * np.nansum(data_even[row,:]))
                 
        else:
             raise ValueError("'even_over_odd' must be True (default) or False")

        return cpi


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


    def pca(self, chain_lengths, scaling_method='z-score'):

        if len(chain_lengths) != len(self.data[0,:]):
            raise ValueError(
                "'chain_lengths' must be the same length as the number of data columns"    
            )
        
        for row in range(0, len(self.data[:,0])):
            if np.sum(self.data[row,:]) == 0:
                raise ValueError(f"Sample in row {row} does not contain any leaf wax chain-length data (concentration or abundances == NaN or 0). Please remove these samples from the input data array before performing PCA.")

        # Apply data scaling before PCA
        match scaling_method:
            case 'z-score':
                data_scaler = StandardScaler()
                data_scaler.fit(self.data)
                data_scaled = data_scaler.transform(self.data)
                data_df_scaled = pd.DataFrame(data=data_scaled, columns=chain_lengths)

            case 'clr':
                data_multi_replace = multi_replace(self.data)
                data_clr = clr(data_multi_replace)
                data_df_scaled = pd.DataFrame(data=data_clr, columns=chain_lengths)

            # case 'ilr':
            #     data_multi_replace = multi_replace(self.data)
            #     data_ilr = ilr(data_multi_replaced)
            #     Figure out how to back-transform ilr matrix for PCA per Filzmoser et al. (2009)
            #     data_df_scaled = pd.DataFrame(data=data_ilr_transform, columns=chain_lengths)

            case None:
                warnings.warn("It is recommended that the user apply a scaling method to their data for Principal Component Analysis.")
                data_df_scaled = pd.DataFrame(data=self.data, columns=chain_lengths)

            case _:
                raise ValueError("'scaling_method' must be set to 'z-score' (default), 'clr', or None")

        # PCA procedure regardless of data scaling method
        data_pca = PCA(n_components=len(chain_lengths))
        data_pca.fit_transform(data_df_scaled)
            
        data_pca_loadings = data_pca.components_
        data_pca_features = data_df_scaled.columns
        data_pca_values = np.arange(data_pca.n_components_) + 1

        pca_dict = {
            "pca": data_pca,
            "pc_values": data_pca_values,
            "features": data_pca_features,
            "loadings": data_pca_loadings
        }

        for i in range(0, len(chain_lengths)):

            data_pc = data_pca.fit_transform(data_df_scaled)[:,i]
            data_scale_pc = 1.0 / (data_pc.max() - data_pc.min())
            data_pc_scores = data_pc * data_scale_pc

            # pca_dict.update({f"wax_pc{i+1}": wax_pc})
            # pca_dict.update({f"wax_scale_pc{i+1}": wax_scale_pc})
            pca_dict.update({f"pc{i+1}_scores": data_pc_scores})

        return pca_dict
