import numpy as np


def validate_data_dimensions(data):
    """
    Raises a TypeError if the data arg is not 2-dimensional.

    Parameters
    ----------
    data : 2-D array-like
        Uses the self.data (input_data) variable in the __init__ for the Chain
        and Isotope API classes.

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
