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
arctic_data_path = os.path.join(DATA_DIR, "Lindberg_Arctic_terrestrial_plantwaxes.csv")
qpt_data_path = os.path.join(DATA_DIR, "LakeQaupatPlantWaxData.csv")

arctic_data_df = pd.read_csv(arctic_data_path)
arctic_acid_chain_df = arctic_data_df[
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

qpt_data_df = pd.read_csv(qpt_data_path)
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

# Create array for chain_lengths arg in Chain.acl and Chain.cpi
arctic_chain_lengths = np.zeros(shape=len(arctic_acid_chain_df.columns))
for col in range(len(arctic_acid_chain_df.columns)):
    col_name = arctic_acid_chain_df.columns[col]
    arctic_chain_lengths[col] = int(col_name[1:3])

qpt_chain_lengths = np.zeros(shape=len(qpt_acid_chain_df.columns))
for col in range(len(qpt_acid_chain_df.columns)):
    col_name = qpt_acid_chain_df.columns[col]
    qpt_chain_lengths[col] = int(col_name[1:3])
    

class TestChainChainInit:
    ''' Test for Chain instantiation '''
    
    def test_init_t0(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        
        assert arctic_acid_chain_obj.data.ndim == 2


    def test_init_t1(self):
        
        with pytest.raises(TypeError):
            arctic_c22_ser = pd.Series(data=arctic_acid_chain_df.c22_fconc)
            arctic_c22_arr = np.array(arctic_c22_ser)
            arctic_c22_obj = Chain(arctic_c22_arr)
        

class TestchainChainTotal_conc:
    ''' Test Chain.total_conc() '''

    def test_total_conc_t0(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log=False)
        
        assert np.round(arctic_acid_total_conc[0], decimals=3) == 41.6
        assert np.round(arctic_acid_total_conc[14], decimals=3) == 0
        assert np.round(arctic_acid_total_conc[-1], decimals=3) == 315.4
        
    
    def test_total_conc_t1(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log=True)
        
        assert np.round(arctic_acid_total_conc[0], decimals=3) == np.round(np.log(41.6), decimals=3)
        assert np.round(arctic_acid_total_conc[14], decimals=3) == np.round(np.log(0), decimals=3)
        assert np.round(arctic_acid_total_conc[-1], decimals=3) == np.round(np.log(315.4), decimals=3)
        

    def test_total_conc_t2(self):
        
        with pytest.raises(ValueError):
            arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
            arctic_acid_total_conc = arctic_acid_chain_obj.total_conc(calculate_log="False")


class TestchainChainRelative_abd:
    ''' Test Chain.relative_abd() '''

    def test_relative_abd_t0(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent=False)
        
        assert np.round(np.sum(arctic_acid_relative_abd[0,:]), decimals=5) == 1
        assert np.isnan(np.sum(arctic_acid_relative_abd[14,:])) == True
        assert np.round(np.sum(arctic_acid_relative_abd[-1,:]), decimals=5) == 1
        
        
    def test_relative_abd_t1(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent=True)
        
        assert np.round(np.sum(arctic_acid_relative_abd[0,:]), decimals=5) == 100
        assert np.isnan(np.sum(arctic_acid_relative_abd[14,:])) == True
        assert np.round(np.sum(arctic_acid_relative_abd[-1,:]), decimals=5) == 100
    
    
    def test_relative_abd_t2(self):
        
        with pytest.raises(ValueError):
            arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
            arctic_acid_relative_abd = arctic_acid_chain_obj.relative_abd(calculate_percent="False")


class TestchainChainAcl:
    ''' Test Chain.acl() '''
    
    def test_acl_t0(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_acl = arctic_acid_chain_obj.acl(chain_lengths=arctic_chain_lengths)
        
        assert np.round(arctic_acid_acl[0], decimals=2) == 24.55
        assert np.isnan(np.sum(arctic_acid_acl[14])) == True
        assert np.round(arctic_acid_acl[-1], decimals=2) == 24.27
        
    
    def test_acl_t1(self):
        
        with pytest.raises(ValueError):
            arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
            arctic_acid_acl = arctic_acid_chain_obj.acl(chain_lengths=np.arange(arctic_chain_lengths[0], arctic_chain_lengths[-1]))


class TestchainChainCpi:
    ''' Test Chain.cpi() '''

    def test_cpi_t0(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd=True)
        
        assert np.round(arctic_acid_cpi[0], decimals=2) == 2.09
        assert arctic_acid_cpi[6] == np.inf
        assert np.isnan(np.sum(arctic_acid_cpi[14])) == True
        assert np.round(arctic_acid_cpi[-1], decimals=2) == 13.72
        
    
    @pytest.mark.filterwarnings("ignore:The first")
    def test_cpi_t1(self):
        
        arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
        arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd=False)
        
        assert np.round(arctic_acid_cpi[0], decimals=2) == 0.36
        assert arctic_acid_cpi[6] == 0
        assert np.isnan(np.sum(arctic_acid_cpi[14])) == True
        assert np.round(arctic_acid_cpi[-1], decimals=2) == 0.07
        
    
    def test_cpi_t2(self):
        
        with pytest.raises(ValueError):
            arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
            arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=np.arange(arctic_chain_lengths[0], arctic_chain_lengths[-1]), even_over_odd=True)


    def test_cpi_t3(self):
    
        with pytest.raises(ValueError):        
            arctic_acid_chain_obj = Chain(arctic_acid_chain_df)
            arctic_acid_cpi = arctic_acid_chain_obj.cpi(chain_lengths=arctic_chain_lengths, even_over_odd="True")
        

class TestchainChainCorrelation_rvals:
    ''' Test Chain.corr_rvals() '''

    def test_correlation_rvals_t0(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_rvals = qpt_acid_chain_obj.correlation_rvals(minimum_obs=2)
        
        assert np.round(qpt_acid_rvals[1,0], decimals=5) == 0.89322
            

    def test_correlation_rvals_t1(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_rvals = qpt_acid_chain_obj.correlation_rvals(minimum_obs=2)
        
        for col in range(len(qpt_acid_rvals[0,:])):
            assert np.round(qpt_acid_rvals[col,col], decimals=5) == 1


    def test_correlation_rvals_t2(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_rvals = qpt_acid_chain_obj.correlation_rvals(minimum_obs=2)
        
        for row in range(len(qpt_acid_rvals[:,0])):
            for col in range(len(qpt_acid_rvals[0,:])):
                assert qpt_acid_rvals[row,col] == qpt_acid_rvals[col,row]
                

class TestchainChainCorrelation_pvals:
    ''' Test Chain.corr_pvals() '''
    
    def test_correlation_pvals_t0(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_pvals = qpt_acid_chain_obj.correlation_pvals(minimum_obs=2)
        
        for col in range(len(qpt_acid_pvals[0,:])):
            assert np.round(qpt_acid_pvals[col,col], decimals=5) == 0
            
            
    def test_correlation_pvals_t1(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_pvals = qpt_acid_chain_obj.correlation_pvals(minimum_obs=2)
        
        for row in range(len(qpt_acid_pvals[:,0])):
            for col in range(len(qpt_acid_pvals[0,:])):
                assert qpt_acid_pvals[row,col] == qpt_acid_pvals[col,row]


class TestchainChainPca:
    ''' Test Chain.pca() '''

    def test_pca_t0(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=qpt_chain_lengths, scaling_method='z-score')


    def test_pca_t1(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=qpt_chain_lengths, scaling_method='clr')


    @pytest.mark.filterwarnings("ignore:It is")
    def test_pca_t2(self):
        
        qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
        qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=qpt_chain_lengths, scaling_method=None)


    def test_pca_t3(self):
        
        with pytest.raises(ValueError):
            qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
            qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=np.arange(qpt_chain_lengths[0], qpt_chain_lengths[-1]))


    def test_pca_t4(self):
        
        with pytest.raises(ValueError):
            qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
            qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=qpt_chain_lengths, scaling_method="True")
            
            
    def test_pca_t5(self):
        
        with pytest.raises(ValueError):
            qpt_acid_chain_df.iloc[0,0] = np.nan
            qpt_acid_chain_obj = Chain(qpt_acid_chain_df)
            qpt_acid_pca = qpt_acid_chain_obj.pca(chain_lengths=qpt_chain_lengths, scaling_method="z-score")
            
