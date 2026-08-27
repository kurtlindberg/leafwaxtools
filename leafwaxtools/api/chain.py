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
from ..utils import preprocessing, correlation


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
    data : 2-D array-like
        User leaf wax chain-length concentration/abundance data.
        
    See also
    --------
    
    leafwaxtools.utils.preprocessing.validate_data_dimensions: Checks to make sure input user data is 2-dimensional.
    
    leafwaxtools.utils.preprocessing.coerce_nan: Converts all non-numeric values to NaNs.
    
    Examples
    --------
    
    .. jupyter-execute::
        
        from leafwaxtools import Chain
        
    """
   
    def __init__(self, input_data):

        preprocessing.validate_data_dimensions(input_data)
        
        self.data = preprocessing.coerce_nan(input_data)


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
        chain_lengths : 1-D array-like
            1-D array-like of integers or floats representing the 
            carbon chain-length number of each column.

        Returns
        -------
        acl : numpy.ndarray
            1-D Numpy array of ACL values for each sample (row).
            
        See also
        --------
        
        leafwaxtools.utils.preprocessing.validate_chain_lengths: Checks to make sure 'chain_lengths' is the same length as the number of user data columns

        """

        preprocessing.validate_chain_lengths(self.data, chain_lengths)

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
        chain_lengths : 1-D array-like
            1-D array-like of integers or floats representing the 
            carbon chain-length number of each column.
        even_over_odd : bool, optional
            Calculates the CPI of even-chain over odd-chain leaf waxes (use 
            case for n-alkanoic acids). Change to False to calculate the CPI 
            of odd-chain over even-chain waxes (use case for n-alkanes). The 
            default is True.

        Raises
        ------
        ValueError
            Raises an error if 'even_over_odd' is neither True nor False.

        Returns
        -------
        cpi : numpy.ndarray
            1-D Numpy array of CPI values for each sample (row).
            
        See also
        --------
        
        leafwaxtools.utils.preprocessing.validate_chain_lengths: Checks to make sure 'chain_lengths' is the same length as the number of user data columns

        """

        preprocessing.validate_chain_lengths(self.data, chain_lengths)

        if chain_lengths[0] % 2 != 0 and even_over_odd is True:
            warnings.warn(
                f"even_over_odd: The first chain-length '{chain_lengths[0]}' is an odd number, but this cpi function will be dividing even number chain-lengths over odd number ones"
            )
        
        if chain_lengths[0] % 2 == 0 and even_over_odd is False:
            warnings.warn(
                f"even_over_odd: The first chain-length '{chain_lengths[0]}' is an even number, but this cpi function will be dividing odd number chain-lengths over even number ones"
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
            
        See also
        --------
        
        leafwaxtools.utils.correlation.corr_r: Calculates the correlation r-values for both the Chain and Isotope API classes.

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
            
        See also
        --------
        
        leafwaxtools.utils.correlation.corr_p: Calculates the correlation p-values for both the Chain and Isotope API classes.

        """
        
        p_vals = correlation.corr_p(data=self.data, min_obs=minimum_obs)
        
        return p_vals


    def pca(self, chain_lengths, scaling_method='z-score', drop_nans=False):
        """
        Performs a Principal Component Analysis (PCA) on the leaf wax 
        chain-length data.
                                                                  
        References:
        Aitchison, J. (1982). The statistical analysis of compositional data. 
        Journal of the Royal Statistical Society: Series B (Methodological), 
        44(2), 139-160. https://doi.org/10.1111/j.2517-6161.1982.tb01195.x
        
        Aton, M., McDonald, D., Cañardo Alastuey, J., Azom, R., Batra, P., 
        Bezshapkin, V., ... & Zhu, Q. (2026). Scikit-bio: a fundamental Python 
        library for biological omic data analysis. Nature Methods, 23(2), 
        274-276. https://doi.org/10.1038/s41592-025-02981-z
        
        Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V., & Egozcue, J. J. 
        (2017). Microbiome datasets are compositional: and this is not 
        optional. Frontiers in microbiology, 8, 2224.
        https://doi.org/10.3389/fmicb.2017.02224
        
        Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., 
        Grisel, O., ... & Duchesnay, É. (2011). Scikit-learn: Machine learning 
        in Python. the Journal of machine Learning research, 12, 2825-2830.

        Parameters
        ----------
        chain_lengths : 1-D array-like
            1-D array-like of integers or floats representing the 
            carbon chain-length number of each column.
        scaling_method : str or None, optional
            Scaling method applied to user data prior to PCA. Available 
            methods include 'z-score' using 
            sklearn.preprocessing.StandardScaler (Pedregosa et al., 2011), 
            'clr' (centered log-ratio transformation; Aitchison, 1982) using 
            skbio.stats.composition.clr (Aton et al., 2026), and None. The 
            default is 'z-score'.
        drop_nans : bool, optional
            Drops all rows of user input data (self.data) containing NaNs. If 
            NaNs are not dropped from the data, sklearn.decomposition.PCA() 
            raise a ValueError. The default is False.

        Raises
        ------
        ValueError
            Raises an error if 'drop_nans' is neither True nor False.

        Returns
        -------
        pca_dict : dict
            A dictionary containing the following keys: "pca" (the full set of
            parameters and returns from the sklearn.decomposition.PCA class 
            after fitting it to the user data), "features" (array of names of 
            each loading provided by 'chain_lengths'), "loadings" (Pandas 
            DataFrame of the vectors/principal component scores of each 
            loading feature (rows) organized by decreasing explained variance 
            (columns)),"scores" (Pandas DataFrame of scores in each principal 
            component (columns) for every input data sample (rows)), 
            "scores_scaled" (Pandas DataFrame of "scores" scaled by the 
            minimum and maximum score of each principal component; useful for 
            creating PCA biplots with loadings and scores set to the same 
            scale).
            
        See also
        --------
        
        leafwaxtools.utils.preprocessing.validate_chain_lengths: Checks to make sure 'chain_lengths' is the same length as the number of user data columns
        
        leafwaxtools.utils.preprocessing.drop_nan: Removes rows (samples) containing NaN values.

        """
        
        preprocessing.validate_chain_lengths(self.data, chain_lengths)
            
        if drop_nans is True:
            data = preprocessing.drop_nan(self.data)
        
        elif drop_nans is False:
            data = self.data
            
        else:
            raise ValueError("'drop_nans' must either be True or False (default)")

        match scaling_method:
            case 'z-score':
                data_scaler = StandardScaler()
                data_scaled = data_scaler.fit_transform(data)
                data_df_scaled = pd.DataFrame(data=data_scaled, columns=chain_lengths)

            case 'clr':
                data_clr = clr(multi_replace(data))
                data_df_scaled = pd.DataFrame(data=data_clr, columns=chain_lengths)

            case None:
                warnings.warn("scaling_method: It is recommended that the user apply a scaling method to their data for Principal Component Analysis.")
                data_df_scaled = pd.DataFrame(data=data, columns=chain_lengths)

            case _:
                raise ValueError("'scaling_method' must be set to 'z-score' (default), 'clr', or None")

        pca_model = PCA(n_components=len(chain_lengths))
        pca_scores = pca_model.fit_transform(data_df_scaled)
        
        pc_columns = [None] * pca_model.n_components_
        for i in range(pca_model.n_components_):
            pc_columns[i] = "PC" + str(i+1)

        pca_dict = {
            "pca": pca_model,
            "features": chain_lengths,
            "loadings": pd.DataFrame(data=pca_model.components_.T, index=chain_lengths, columns=pc_columns),
            "scores": pd.DataFrame(data=pca_scores, columns=pc_columns)
        }

        pca_scores_scaled = np.zeros(shape=np.shape(pca_scores))
        for col in range(pca_model.n_components_):
            pca_scores_col = pca_scores[:,col]
            pca_scores_col_scale = 1.0 / (pca_scores_col.max() - pca_scores_col.min())
            pca_scores_scaled[:,col] = pca_scores_col * pca_scores_col_scale            
        
        pca_dict.update({"scores_scaled": pd.DataFrame(data=pca_scores_scaled, columns=pc_columns)})
            
        return pca_dict
