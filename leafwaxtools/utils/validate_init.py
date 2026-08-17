import numpy as np


def validate_data_dimensions(data):
    
    data_arr = np.array(data)
    
    if data_arr.ndim != 2:
        raise TypeError("'input_data' must be 2-dimensional")