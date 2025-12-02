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
from leafwaxtools import Chain, Isotope

# Path to test data
DATA_DIR = Path(__file__).parents[1].joinpath("data").resolve()
data_path = os.path.join(DATA_DIR, "Lindberg_Arctic_terrestrial_plantwaxes.csv")


class TestChainChainInit:
    ''' Test for Chain instantiation '''

    # @pytest.mark.xfail
    def test_init_t0(self):
        test_df = pd.read_excel(data_path)
        test_data = Chain(np.array(test_df))

        assert type(test_data.data) == pd.core.frame.DataFrame
        #assert WaxData.data == test_data


class TestChainChainTotal_conc:
    ''' Test Chain.total_conc() '''

    # def test_total_conc_t0(self):


class TestChainChainRelative_abd:
    ''' Test Chain.relative_abd() '''

    # def test_relative_abd_t0(self):


class TestChainChainAcl:
    ''' Test Chain.acl() '''

    # def test_acl_t0(self):


class TestChainChainCpi:
    ''' Test Chain.cpi() '''

    # def test_cpi_t0(self):


class TestChainChainCorr_rvals:
    ''' Test Chain.corr_rvals() '''

    # def test_corr_rvals_t0(self):


class TestChainChainCorr_pvals:
    ''' Test Chain.corr_pvals() '''

    # def test_corr_pvals_t0(self):


class TestChainChainPca:
    ''' Test Chain.pca() '''

    # def test_pca_t0(self):
