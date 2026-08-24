import pandas as pd
import numpy as np
import warnings


def drop_nans(data):
    
    data_df = pd.DataFrame(data=data)
    data_df.dropna(inplace=True)
    
    data_arr = np.array(data_df)
    
    rows_dropped = len(data[:,0]) - len(data_arr[:,0])
    warnings.warn(f"drop_nans: {rows_dropped} row(s)/sample(s) removed due to having NaN values")
    
    return data_arr