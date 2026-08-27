import pandas as pd
import numpy as np
import warnings


def validate_data_dimensions(data):
    """
    Raises a TypeError if 'data' is not 2-dimensional.

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) attribute from the Chain and Isotope 
        API classes.

    Raises
    ------
    TypeError
        Raises and error when Numpy array-transformed 'data' is not 
        2-dimensional.

    Returns
    -------
    None.

    """
    
    data_arr = np.array(data)
    
    if data_arr.ndim != 2:
        raise TypeError("'input_data' must be 2-dimensional")
        

def validate_chain_lengths(data, chain):
    """
    Raises a ValueError if 'chain' is not the same length as the number of 'data' columns.

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) attribute from the Chain and Isotope 
        API classes.
    chain : 1-D array-like
        1-D array-like of integers or floats representing the 
        carbon chain-length number of each column.

    Raises
    ------
    ValueError
        Raises an error when 'chain' is not the same length as number of 
        columns (2nd dimension) as 'data'.

    Returns
    -------
    None.

    """
    
    if len(chain) != len(data[0,:]):
        raise ValueError(
            "'chain_lengths' must be the same length as the number of data columns"    
        )
    

def coerce_nan(data):
    """
    Converts all non-numeric values in 'data' into NaN-like types 
    (numpy.float64(nan)). 

    Parameters
    ----------
    data : 2-D array-like
        Uses the 'input_data' during the instantiation (__init__) of the Chain 
        and Isotope API classes.

    Returns
    -------
    data_arr : numpy.ndarray
        2-D Numpy array of 'data' where all non-numeric values have been 
        converted to the same NaN-like type (numpy.float64(nan)).

    """
    
    data_df = pd.DataFrame(data=data)
    data_dfnan = data_df.apply(pd.to_numeric, errors='coerce')
    data_arr = np.array(data_dfnan)
    
    return data_arr


def drop_nan(data):
    """
    Removes all rows from 'data' that contain any NaN values. Prints the 
    number of rows dropped, if greater than 0, to the console.

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) attribute from the Chain and Isotope 
        API classes.

    Returns
    -------
    data_arr : numpy.ndarray
        2-D Numpy array of 'data' with all rows containing any NaN values 
        removed.

    """
    
    data_og = np.array(data)
    data_df = pd.DataFrame(data=data)
    data_df.dropna(inplace=True)
    
    data_arr = np.array(data_df)
    
    rows_dropped = len(data_og[:,0]) - len(data_arr[:,0])
    if rows_dropped > 0:
        warnings.warn(f"drop_nans: {rows_dropped} row(s)/sample(s) removed due to having NaN values")
    
    return data_arr
