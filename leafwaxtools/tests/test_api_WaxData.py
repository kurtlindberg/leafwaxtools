"""
Tests for Cat Class
"""

''' Tests for pyCatSim.api.cat.Cat

Naming rules:
1. class: Test{filename}{Class}{method} with appropriate camel case
2. function: test_{method}_t{test_id}
Notes on how to test:
0. Make sure [pytest](https://docs.pytest.org) has been installed: `pip install pytest`
1. execute `pytest {directory_path}` in terminal to perform all tests in all testing files inside the specified directory
2. execute `pytest {file_path}` in terminal to perform all tests in the specified file
3. execute `pytest {file_path}::{TestClass}::{test_method}` in terminal to perform a specific test class/method inside the specified file
4. after `pip install pytest-xdist`, one may execute "pytest -n 4" to test in parallel with number of workers specified by `-n`
5. for more details, see https://docs.pytest.org/en/stable/usage.html
'''

import os
from pathlib import Path
import pytest
import pandas as pd
import numpy as np
from leafwaxtools import WaxData

# Path to test data
DATA_DIR = Path(__file__).parents[1].joinpath("data").resolve()
data_path = os.path.join(DATA_DIR, 'Hollister_et_al_2022_leafwax_data.xlsx')

class TestwaxdataWaxDataInit:
    ''' Test for WaxData instantiation '''
    
    def test_init_t0(self):
        test_df = pd.read_excel(data_path)
        test_data = WaxData(test_df)
        
        assert type(test_data.data) == pd.core.frame.DataFrame
        #assert WaxData.data == test_data
        
    @pytest.mark.xfail
    def test_init_t1(self):
       wax_df = pd.read_excel(data_path)
       wax_arr = np.array(wax_df)
       wax_data = WaxData(wax_arr)
       
       assert type(wax_data.data) == pd.core.frame.DataFrame
        
        
