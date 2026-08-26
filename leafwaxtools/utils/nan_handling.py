import pandas as pd
import numpy as np
import warnings


def drop_nan(data):
    
    data_og = np.array(data)
    data_df = pd.DataFrame(data=data)
    data_df.dropna(inplace=True)
    
    data_arr = np.array(data_df)
    
    rows_dropped = len(data_og[:,0]) - len(data_arr[:,0])
    if rows_dropped > 0:
        warnings.warn(f"drop_nans: {rows_dropped} row(s)/sample(s) removed due to having NaN values")
    
    return data_arr


def coerce_nan(data):
    
    data_df = pd.DataFrame(data=data)
    data_dfnan = data_df.apply(pd.to_numeric, errors='coerce')
    data_arr = np.array(data_dfnan)
    
    return data_arr
