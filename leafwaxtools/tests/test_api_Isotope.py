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
arctic_data_path = os.path.join(DATA_DIR, "Lindberg_Arctic_terrestrial_plantwaxes.csv")
qpt_data_path = os.path.join(DATA_DIR, "LakeQaupatPlantWaxData.csv")

arctic_data_df = pd.read_csv(arctic_data_path)
arctic_acid_iso_df = arctic_data_df[
    [
        'c20_fd2h',
        'c22_fd2h',
        'c24_fd2h',
        'c26_fd2h',
        'c28_fd2h',
        'c30_fd2h',
        'c32_fd2h',
    ]
]

arctic_acid_chain_df = arctic_data_df[
    [
        'c20_fconc',
        'c22_fconc',
        'c24_fconc',
        'c26_fconc',
        'c28_fconc',
        'c30_fconc',
        'c32_fconc',
    ]
]

qpt_data_df = pd.read_csv(qpt_data_path)
qpt_acid_iso_df = qpt_data_df[
    [
         'c20d2h',
         'c22d2h',
         'c24d2h',
         'c26d2h',
         'c28d2h',
         'c30d2h',
    ]
]

qpt_acid_chain_df = qpt_data_df[
    [
        'c20concentration',
        'c22concentration',
        'c24concentration',
        'c26concentration',
        'c28concentration',
        'c30concentration',
    ]
]


class TestisotopeIsotopeInit:
    ''' Test Isotope instantiation '''
    
    def test_init_t0(self):
        
        arctic_acid_iso_obj = Isotope(arctic_acid_iso_df)
        
        assert arctic_acid_iso_obj.data.ndim == 2


    def test_init_t1(self):
        
        with pytest.raises(TypeError):
            arctic_c22_ser = pd.Series(data=arctic_acid_iso_df.c22_fd2h)
            arctic_c22_arr = np.array(arctic_c22_ser)
            arctic_c22_obj = Isotope(arctic_c22_arr)


class TestisotopeIsotopeValue_range:
    ''' Test Isotope.value_range() '''

    def test_value_range_t0(self):
        
        arctic_acid_iso_obj = Isotope(arctic_acid_iso_df)
        arctic_acid_iso_range = arctic_acid_iso_obj.value_range()
        
        assert np.isnan(arctic_acid_iso_range[0]) == True
        assert arctic_acid_iso_range[6] == 15
        assert arctic_acid_iso_range[-2] == 14
        assert np.isnan(arctic_acid_iso_range[-1]) == True


class TestisotopeIsotopeConcentration_avg:
    ''' Test Isotope.concentration_avg() '''

    def test_concentration_avg_t0(self):
        
        arctic_acid_iso_obj = Isotope(arctic_acid_iso_df)
        arctic_acid_conc_avg = arctic_acid_iso_obj.concentration_avg(chain_data=arctic_acid_chain_df)
        
        assert np.isnan(arctic_acid_conc_avg[0]) == True
        assert np.round(arctic_acid_conc_avg[6], decimals=3) == -216.111
        assert arctic_acid_conc_avg[-2] == -211
        assert np.isnan(arctic_acid_conc_avg[-1]) == True
        
        

    def test_concentration_avg_t1(self):
        
        with pytest.raises(ValueError):
            arctic_acid_iso_obj = Isotope(arctic_acid_iso_df)
            arctic_acid_somechain_df = arctic_data_df[
                [
                    'c22_fconc',
                    'c24_fconc',
                    'c26_fconc',
                    'c28_fconc'
                ]
            ]
            arctic_acid_conc_avg = arctic_acid_iso_obj.concentration_avg(chain_data=arctic_acid_somechain_df)
        

# class TestisotopeIsotopeEpsilon:
#     Test Isotope.epsilon()

#     # def test_epsilon_t0(self):

    
# class TestisotopeIsotopeWax_to_source:
#     Test Isotope.wax_to_source()

#     # def test_wax_to_source_t0(self):
    

class TestisotopeIsotopeCorrelation_rvals:
    ''' Test Isotope.corr_rvals() '''

    def test_correlation_rvals_t0(self):
        
        qpt_acid_iso_obj = Isotope(qpt_acid_iso_df)
        qpt_acid_rvals = qpt_acid_iso_obj.correlation_rvals(minimum_obs=2)
        
        assert np.round(qpt_acid_rvals[1,0], decimals=5) == 0.46539
            

    def test_correlation_rvals_t1(self):
        
        qpt_acid_iso_obj = Isotope(qpt_acid_iso_df)
        qpt_acid_rvals = qpt_acid_iso_obj.correlation_rvals(minimum_obs=2)
        
        for col in range(len(qpt_acid_rvals[0,:])):
            assert np.round(qpt_acid_rvals[col,col], decimals=5) == 1


    def test_correlation_rvals_t2(self):
        
        qpt_acid_iso_obj = Isotope(qpt_acid_iso_df)
        qpt_acid_rvals = qpt_acid_iso_obj.correlation_rvals(minimum_obs=2)
        
        for row in range(len(qpt_acid_rvals[:,0])):
            for col in range(len(qpt_acid_rvals[0,:])):
                assert qpt_acid_rvals[row,col] == qpt_acid_rvals[col,row]
    
    
class TestisotopeIsotopeCorrelation_pvals:
    ''' Test Chain.corr_pvals() '''
    
    def test_correlation_pvals_t0(self):
        
        qpt_acid_iso_obj = Isotope(qpt_acid_iso_df)
        qpt_acid_pvals = qpt_acid_iso_obj.correlation_pvals(minimum_obs=2)
        
        for col in range(len(qpt_acid_pvals[0,:])):
            assert np.round(qpt_acid_pvals[col,col], decimals=5) == 0
            
            
    def test_correlation_pvals_t1(self):
        
        qpt_acid_iso_obj = Isotope(qpt_acid_iso_df)
        qpt_acid_pvals = qpt_acid_iso_obj.correlation_pvals(minimum_obs=2)
        
        for row in range(len(qpt_acid_pvals[:,0])):
            for col in range(len(qpt_acid_pvals[0,:])):
                assert qpt_acid_pvals[row,col] == qpt_acid_pvals[col,row]
