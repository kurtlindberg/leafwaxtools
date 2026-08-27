import scipy.stats
import numpy as np
import pandas as pd


def corr_r(data, min_obs):
    """
    Calculates the Pearson correlation r-values between each leaf wax 
    chain-length (columns). To be extended with other correlation methods 
    (Spearman, Kendall Tau) in a future version. This function utilises the 
    correlation features from SciPy (Virtanen et al., 2020). Additionally, this 
    function drops all NaN values from the two data columns used in each 
    correlation, similar to the implementation of pandas.DataFrame.corr().

    References:
    
    Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., 
    Cournapeau, D., ... & Van Mulbregt, P. (2020). SciPy 1.0: fundamental 
    algorithms for scientific computing in Python. Nature methods, 17(3), 
    261-272.
    https://doi.org/10.1038/s41592-019-0686-2

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) variable from the API class' __init__.
    min_obs : int
        Minimum number of observations (samples/rows) required to return a
        Pearson r-value..

    Returns
    -------
    corr_rvals : numpy.ndarray
        2-D Numpy array of Pearson correlation r-values between each leaf 
        wax chain-length (column) with all values in the major diagonal 
        equal to 1.

    """
    
    data_df = pd.DataFrame(data=data)
    corr_rvals = np.zeros((len(data[0,:]), len(data[0,:])))

    for row in data_df.columns:
        for col in data_df.columns:
            if col == row:
                data_df_corr = data_df[[row]].dropna()
            else:
                data_df_corr = data_df[[row, col]].dropna()
            
            if (len(data_df_corr[row]) >= min_obs) and (len(data_df_corr[col]) >= min_obs):
                corr_rvals[row,col] = scipy.stats.pearsonr(data_df_corr[row], data_df_corr[col])[0]
            else:
                corr_rvals[row,col] = np.nan
    
    return corr_rvals

    
def corr_p(data, min_obs):
    """
    Calculates the Pearson correlation p-values between each leaf wax 
    chain-length (columns). To be extended with other correlation methods 
    (Spearman, Kendall Tau) in a future version. This function utilises the 
    correlation features from SciPy (Virtanen et al., 2020). Additionally, this 
    function drops all NaN values from the two data columns used in each 
    correlation, similar to the implementation of pandas.DataFrame.corr().

    References:
    
    Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., 
    Cournapeau, D., ... & Van Mulbregt, P. (2020). SciPy 1.0: fundamental 
    algorithms for scientific computing in Python. Nature methods, 17(3), 
    261-272.
    https://doi.org/10.1038/s41592-019-0686-2

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) variable from the API class' __init__.
    min_obs : int
        Minimum number of observations (samples/rows) required to return a
        Pearson p-value.

    Returns
    -------
    corr_pvals : numpy.ndarray
        2-D Numpy array of Pearson correlation p-values between each leaf 
        wax chain-length (column).

    """
    
    data_df = pd.DataFrame(data=data)
    corr_pvals = np.zeros((len(data[0,:]), len(data[0,:])))

    for row in data_df.columns:
        for col in data_df.columns:
            if col == row:
                data_df_corr = data_df[[row]].dropna()
            else:
                data_df_corr = data_df[[row, col]].dropna()
            
            if (len(data_df_corr[row]) >= min_obs) and (len(data_df_corr[col]) >= min_obs):
                corr_pvals[row,col] = scipy.stats.pearsonr(data_df_corr[row], data_df_corr[col])[1]
            else:
                corr_pvals[row,col] = np.nan
    
    return corr_pvals
