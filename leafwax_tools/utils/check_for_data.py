"""
This module checks the input_data DataFrame in class WaxData for appropriate
column names denoting leaf wax compound classes (FAMEs/n-alknaoic acids or
n-alkanes) and data types (chain-length concentration, d2H, d13C)
"""

import warnings

def check_for_data(input_data):
    
    f_data = input_data.filter(regex="f")
    if len(f_data.columns) == 0:
        warnings.warn("Warning: column names do not label FAMEs (n-alkanoic acids) compound class with an 'f'; i.e. c20_fconc")
        
    a_data = input_data.filter(regex="a")
    if len(a_data.columns) == 0:
        warnings.warn("Warning: column names do not label n-alkanes compound class with an 'a'\n; i.e. c20_aconc")

    # add checks for column names