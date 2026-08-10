"""
Tests for Chain Class
"""

''' Tests for leafwaxtools.api.chain.Chain

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
from leafwaxtools import Chain

# Path to test data
DATA_DIR = Path(__file__).parents[1].joinpath("data").resolve()
data_path = os.path.join(DATA_DIR, "Lindberg_Arctic_terrestrial_plantwaxes.csv")

arctic_df = pd.read_csv(data_path)
arctic_acid_chain_df = arctic_df[
    [
        'c20_fconc',
        'c21_fconc',
        'c22_fconc',
        'c23_fconc',
        'c24_fconc',
        'c25_fconc',
        'c26_fconc',
        'c27_fconc',
        'c28_fconc',
        'c29_fconc',
        'c30_fconc',
        'c31_fconc',
        'c32_fconc',
    ]
]

# Create array for chain_lengths arg in Chain.acl and Chain.cpi
arctic_chain_lengths = np.zeros(shape=len(arctic_acid_chain_df.columns))
for col in range(len(arctic_acid_chain_df.columns)):
    col_name = arctic_acid_chain_df.columns[col]
    arctic_chain_lengths[col] = int(col_name[1:3])


class TestChainChainInit:
    ''' Test for Chain instantiation '''

    def test_init_t0(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        assert arctic_acid_chain_obj.data.all() == arctic_acid_chain_arr.all()
    
    
    def test_init_t1(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        assert arctic_acid_chain_obj.data.ndim == 2


    @pytest.mark.xfail
    def test_init_t2(self):
        
        arctic_c22_ser = pd.Series(data=arctic_acid_chain_df.c22_fconc)
        arctic_c22_arr = np.array(arctic_c22_ser)
        arctic_c22_obj = Chain(arctic_c22_arr)
        
        assert arctic_c22_obj.data.ndim == 2
        

class TestChainChainTotal_conc:
    ''' Test Chain.total_conc() '''

    def test_total_conc_t0(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log=False)
        
        assert np.round(arctic_acid_total_conc[0], decimals=3) == 41.6
        assert np.round(arctic_acid_total_conc[14], decimals=3) == 0
        assert np.round(arctic_acid_total_conc[-1], decimals=3) == 315.4
        
    
    def test_total_conc_t1(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log=True)
        
        assert np.round(arctic_acid_total_conc[0], decimals=3) == np.round(np.log(41.6), decimals=3)
        assert np.round(arctic_acid_total_conc[14], decimals=3) == np.round(np.log(0), decimals=3)
        assert np.round(arctic_acid_total_conc[-1], decimals=3) == np.round(np.log(315.4), decimals=3)
        

    @pytest.mark.xfail
    def test_total_conc_t2(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log="False")


class TestChainChainRelative_abd:
    ''' Test Chain.relative_abd() '''

    def test_relative_abd_t0(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent=False)
        
        assert np.round(np.sum(arctic_acid_relative_abd[0,:]), decimals=5) == 1
        assert np.isnan(np.sum(arctic_acid_relative_abd[14,:])) == True
        assert np.round(np.sum(arctic_acid_relative_abd[-1,:]), decimals=5) == 1
        
        
    def test_relative_abd_t1(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent=True)
        
        assert np.round(np.sum(arctic_acid_relative_abd[0,:]), decimals=5) == 100
        assert np.isnan(np.sum(arctic_acid_relative_abd[14,:])) == True
        assert np.round(np.sum(arctic_acid_relative_abd[-1,:]), decimals=5) == 100
    
    
    @pytest.mark.xfail
    def test_relative_abd_t2(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent="False")


class TestChainChainAcl:
    ''' Test Chain.acl() '''
    
    def test_acl_t0(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_acl = arctic_acid_chain_obj.acl(chain_lengths=arctic_chain_lengths)
        
        assert np.round(arctic_acid_acl[0], decimals=2) == 24.55
        assert np.isnan(np.sum(arctic_acid_acl[14])) == True
        assert np.round(arctic_acid_acl[-1], decimals=2) == 24.27
        
    
    @pytest.mark.xfail
    def test_acl_t1(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_acl = arctic_acid_chain_obj.acl(chain_lengths=np.arange(arctic_chain_lengths[0], arctic_chain_lengths[-1]))


class TestChainChainCpi:
    ''' Test Chain.cpi() '''

    def test_cpi_t0(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd=True)
        
        assert np.round(arctic_acid_cpi[0], decimals=2) == 2.09
        assert arctic_acid_cpi[6] == np.inf
        assert np.isnan(np.sum(arctic_acid_cpi[14])) == True
        assert np.round(arctic_acid_cpi[-1], decimals=2) == 13.72
        
        
    def test_cpi_t1(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd=False)
        
        assert np.round(arctic_acid_cpi[0], decimals=2) == 0.36
        assert arctic_acid_cpi[6] == 0
        assert np.isnan(np.sum(arctic_acid_cpi[14])) == True
        assert np.round(arctic_acid_cpi[-1], decimals=2) == 0.07
        
    
    @pytest.mark.xfail
    def test_cpi_t2(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=np.arange(arctic_chain_lengths[0], arctic_chain_lengths[-1]), even_over_odd=True)


    @pytest.mark.xfail
    def test_cpi_t3(self):
        
        arctic_acid_chain_arr = np.array(arctic_acid_chain_df)
        arctic_acid_chain_obj = Chain(arctic_acid_chain_arr)
        
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd="True")
        

# class TestChainChainCorr_rvals:
#     Test Chain.corr_rvals()

#     # def test_corr_rvals_t0(self):


# class TestChainChainCorr_pvals:
#     Test Chain.corr_pvals()

#     # def test_corr_pvals_t0(self):


# class TestChainChainPca:
#     Test Chain.pca()

#     # def test_pca_t0(self):
