"""
This module checks the input_data DataFrame in class WaxData for appropriate
column names denoting leaf wax compound classes (FAMEs/n-alknaoic acids or
n-alkanes) and data types (chain-length concentration, d2H, d13C)
"""


def validate_data(input_data, data_type):
    
    if ((data_type != "f") and (data_type != "a")):
        raise ValueError("""data_type must either be "f" (FAMEs/n-alkanoic acids) or "a" (n-alkanes)""")
    
    data_filter = input_data.filter(regex=data_type)
    if len(data_filter.columns) == 0:
        raise ValueError(f"Column names do not label compound class with an '{data_type}' i.e. c24_{data_type}conc, {data_type}c24")
    
    # add checks for column names
