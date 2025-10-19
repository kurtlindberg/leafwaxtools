"""
The Chain module is the class for performing calculations using plant wax
chain-length concentration/abundance data imported as a 2D array-like object
"""

from enum import Enum
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from composition_stats import clr, closure, multiplicative_replacement
import warnings
from ..utils import validate_data
from ..utils.data_type_enum import DataType


class Chain:


    def __init__(self, input_data):

        self.data = input_data


    def total_conc(self, should_calc_log=False):

        total_conc = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            total_conc[row] = np.nansum(self.data[row,:])

            if total_conc[row] == 0:
                total_conc[row] = np.nan

        match should_calc_log:
            case True:
                total_conc = np.log(total_conc)
            case False:
                total_conc = total_conc

        return total_conc


    def relative_abd(self):

        rel_abd = np.zeros(np.shape(self.data))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):

                rel_abd[row,col] = self.data[row,col]/np.sum(self.data[row,:])

        return rel_abd


    def acl(self, start_chain, end_chain, data_type=DataType.NALKANOIC_ACID, use_all_chains=True):

        chain_lengths = list(range(start_chain, end_chain+1))
        chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
        chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]

        data = pd.DataFrame(data=self.data, columns=(map(str, chain_lengths)))

        match use_all_chains:
            case True:
                chain_lengths = chain_lengths

            case False:
                match data_type:
                    case DataType.NALKANOIC_ACID:
                        data = np.array(data.filter(items=(map(str, chain_lengths_even))))
                        chain_lengths = chain_lengths_even

                    case DataType.NALKANE:
                        data = np.array(data.filter(items=(map(str, chain_lengths_odd))))
                        chain_lengths = chain_lengths_odd

        acl_numer = np.zeros(len(data[:,0]))
        acl = np.zeros(len(data[:,0]))

        for row in range(0, len(data[:,0])):
            for col in range(0, len(data[0,:])):

                acl_numer[row] += data[row,col] * chain_lengths[col]

            acl[row] = acl_numer[row]/np.sum(data[row,:])

        return acl


    def cpi(self, start_chain, end_chain, data_type=DataType.NALKANOIC_ACID):

        chain_lengths = list(range(start_chain, end_chain+1))
        chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
        chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]

        data = pd.DataFrame(data=self.data, columns=(map(str, chain_lengths)))
        data_even = np.array(data.filter(items=(map(str, chain_lengths_even))))
        data_odd = np.array(data.filter(items=(map(str, chain_lengths_odd))))
        cpi = np.zeros(len(self.data[:,0]))

        match data_type:
            case DataType.NALKANOIC_ACID:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_even[row,0:-1]) + np.nansum(data_even[row,1:])) / (2 * np.nansum(data_odd[row,:]))

                return cpi

            case DataType.NALKANE:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_odd[row,0:-1]) + np.nansum(data_odd[row,1:])) / (2 * np.nansum(data_even[row,:]))

                return cpi


    '''
    def pca(self, start_chain, end_chain, use_all_chains=True, use_clr=False):

        chain_lengths = list(range(start_chain, end_chain+1))

        match use_all_chains:
            case True:
                chain_lengths = chain_lengths

            case False:
                match data_type:
                    case DataType.NALKANOIC_ACID:
                        chain_lengths = [num for num in chain_lengths if num % 2 == 0]

                    case DataType.NALKANE:
                        chain_lengths = [num for num in chain_lengths if num % 2 == 1]

        # add code that drops rows with all missing data

        match use_clr:
            case True:
                wax_relabd = closure(multiplicative_replacement(self.data))
                wax_clr = clr(wax_relabd)
                wax_data = pd.DataFrame(data=wax_clr, columns=[])

            case False:
                wax_data = pd.DataFrame(data=relative_abd(self.data), columns=[])

        wax_scaler = StandardScaler()
        wax_scaler.fit(wax_data)
        wax_data_scaled = wax_scaler.transform(wax_data)

        wax_pca = PCA(n_components=len(chain_lengths))
        wax_PC_scores = pd.DataFrame(
            wax_pca.fit_transform(wax_data_scaled),
            columns=[]
        )
        wax_loadings = pd.DataFrame(
            wax_pca.components_.T,
            columns=[],
            index=wax_data.columns
        )

        for i in range(0, len(chain_lengths)):
            # figure out how to create new, sequential variables names in loop
            wax_pc{i+1} = wax_pca.fit_transform(wax_data_scaled)[:,i]
            wax_scale_pc{i+1} = 1.0/(wax_pc{i+1}.max() - wax_pc{i+1}.min())
            wax_pc{i+1}_scores = pd.DataFrame(data=(wax_pc{i+1} * wax_scale_pc{i+1}), columns=[f'pc{i+1}'])

        wax_ldings = wax_pca.components_
        wax_features = wax_data.columns
        wax_pc_values = np.arange(wax_pca.n_components_) + 1

        pca_dict = {}

        return pca_dict
        '''
