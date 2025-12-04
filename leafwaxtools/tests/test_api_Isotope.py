"""
Tests for Isotope Class
"""

''' Tests for leafwaxtools.api.isotope.Isotope

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
from leafwaxtools import Isotope

# Path to test data
DATA_DIR = Path(__file__).parents[1].joinpath("data").resolve()
data_path = os.path.join(DATA_DIR, "gorbey2021qpt.csv")


class TestIsotopeIsotopeInit:
    ''' Test for Isotope instantiation '''

    def test_init_t0(self):
        
        qpt_df = pd.read_csv(data_path)
        qpt_isotope_df = qpt_df[
            [
                'c22d2h',
                'c24d2h',
                'c26d2h',
                'c28d2h'
            ]
        ]
        qpt_isotope_arr = np.array(qpt_isotope_df)
        qpt_isotope_obj = Isotope(qpt_isotope_arr)
        
        assert qpt_isotope_obj.data.all() == qpt_isotope_arr.all()
        
    
    def test_init_t1(self):
        
        qpt_df = pd.read_csv(data_path)
        qpt_isotope_df = qpt_df[
            [
                'c22d2h',
                'c24d2h',
                'c26d2h',
                'c28d2h'
            ]
        ]
        qpt_isotope_arr = np.array(qpt_isotope_df)
        qpt_isotope_obj = Isotope(qpt_isotope_arr)
        
        assert qpt_isotope_obj.data.ndim == 2
    
        
    @pytest.mark.xfail    
    def test_init_t2(self):
        
        qpt_df = pd.read_csv(data_path)
        qpt_isotope_df = qpt_df[
            [
                'c22d2h',
            ]
        ]
        qpt_isotope_arr = np.array(qpt_isotope_df)
        qpt_isotope_obj = Isotope(qpt_isotope_arr)


# class TestIsotopeIsotopeIso_range:
#     Test for Isotope.iso_range()

#     # def test_iso_range_t0(self):


# class TestIsotopeIsotopeIso_avg:
#     Test for Isootpe.iso_avg()

#     # def test_iso_avg_t0(self):


# class TestIsotopeIsotopeEpsilon:
#     Test for Isotope.epsilon()

#     # def test_epsilon_t0(self):
