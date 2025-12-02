"""
The Chain module is the class for performing calculations using plant wax
chain-length concentration/abundance data imported as a 2D array-like object
"""


import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from composition_stats import clr, closure, multiplicative_replacement
import scipy.stats
# import warnings
# from ..utils import validate_data


class Chain:


    def __init__(self, input_data):

        self.data = input_data


    def total_conc(self, calculate_log=False, zero_total=0):

        total_conc = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            total_conc[row] = np.nansum(self.data[row,:])

            if total_conc[row] == 0:
                total_conc[row] = zero_total

        match calculate_log:
            case True:
                total_conc = np.log(total_conc)
            case False:
                total_conc = total_conc
            case _:
                raise ValueError("'calculate_log' must either be True or False (default)")

        return total_conc


    def relative_abd(self, calculate_percent=False):

        rel_abd = np.zeros(np.shape(self.data))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):

                rel_abd[row,col] = self.data[row,col]/np.sum(self.data[row,:])
                
        match calculate_percent:
            case True:
                for row in range(0, len(self.data[:,0])):
                    for col in range(0, len(self.data[0,:])):
                        rel_abd[row,col] = rel_abd[row,col]*100
                
            case False:
                rel_abd = rel_abd
            
            case _:
                raise ValueError("'calculate_percent' must either be True or False (default)")

        return rel_abd


    def acl(self, chain_lengths):

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 24, 26, 28]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        acl_numer = np.zeros(len(self.data[:,0]))
        acl = np.zeros(len(self.data[:,0]))

        for row in range(0, len(self.data[:,0])):
            for col in range(0, len(self.data[0,:])):

                acl_numer[row] += self.data[row,col] * chain_lengths[col]

            acl[row] = acl_numer[row]/np.sum(self.data[row,:])

        return acl


    def cpi(self, chain_lengths, even_over_odd=True):

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 23, 24, 25, 26]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        '''
        EKT: use warnings to flag if even over odd order is wrong
        '''

        chain_lengths_even = [num for num in chain_lengths if num % 2 == 0]
        chain_lengths_odd = [num for num in chain_lengths if num % 2 == 1]

        data = pd.DataFrame(data=self.data, columns=(map(str, chain_lengths)))
        data_even = np.array(data.filter(items=(map(str, chain_lengths_even))))
        data_odd = np.array(data.filter(items=(map(str, chain_lengths_odd))))
        cpi = np.zeros(len(self.data[:,0]))

        match even_over_odd:
            case True:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_even[row,0:-1]) + np.nansum(data_even[row,1:])) / (2 * np.nansum(data_odd[row,:]))

            case False:
                for row in range(0, len(self.data[:,0])):
                    cpi[row] = (np.nansum(data_odd[row,0:-1]) + np.nansum(data_odd[row,1:])) / (2 * np.nansum(data_even[row,:]))

            case _:
                raise ValueError("'even_over_odd' must be True (default) or False")

        return cpi


    def corr_rvals(self, minimum_obs=2):

        r_vals = np.zeros((len(self.data[0,:]), len(self.data[0,:])))

        for row in range(0, len(r_vals[:,0])):
            for col in range(0, len(r_vals[0,:])):

                x_corr = np.array(self.data[:,row])
                y_corr = np.array(self.data[:,col])

                if (len(x_corr) >= minimum_obs) and (len(y_corr) >= minimum_obs):
                    r_vals[row,col] = scipy.stats.pearsonr(x_corr, y_corr)[0]
                else:
                    r_vals[row,col] = np.nan

        return r_vals


    def corr_pvals(self, minimum_obs=2):

        p_vals = np.zeros((len(self.data[0,:]), len(self.data[0,:])))

        for row in range(0, len(p_vals[:,0])):
            for col in range(0, len(p_vals[0,:])):

                x_corr = np.array(self.data[:,row])
                y_corr = np.array(self.data[:,col])

                if (len(x_corr) >= minimum_obs) and (len(y_corr) >= minimum_obs):
                    p_vals[row,col] = scipy.stats.pearsonr(x_corr, y_corr)[1]
                else:
                    p_vals[row,col] = np.nan

        return p_vals



    def pca(self, chain_lengths, use_clr=True):

        if type(chain_lengths) is not type(list()):
            raise TypeError(
                "'chain_lengths' must be a list() type containing integers or floats; Example: [22, 23, 24, 25, 26]"
            )

        if len(chain_lengths) < 1:
            raise ValueError(
                "'chain_lengths' is currently an empty list. Please make sure 'chain_lengths' contains at least 1 integer or float."
            )

        '''
        # deal with missing data
        for row in range(0, len(self.data[:,0]):
            for col i range(0, len(self.data[0,:]):

                if self.data[row,col] == np.nan:
                    self.data[row,col] = 0

            if np.sum(self.data[row,:]) == 0:
        '''


        match use_clr:
            case True:
                wax_relabd = closure(multiplicative_replacement(self.data))
                wax_clr = clr(wax_relabd)
                wax_data = pd.DataFrame(data=wax_clr, columns=chain_lengths)

            case False:
                wax_data = pd.DataFrame(data=self.relative_abd(), columns=chain_lengths)
                
            case _:
                raise ValueError("'use_clr' must be True or False (default)")

        wax_scaler = StandardScaler()
        wax_scaler.fit(wax_data)
        wax_data_scaled = wax_scaler.transform(wax_data)

        wax_pca = PCA(n_components=len(chain_lengths))
        wax_pca.fit_transform(wax_data_scaled)

        # wax_PC_scores = pd.DataFrame(
        #     wax_pca.fit_transform(wax_data_scaled),
        #     columns=chain_lengths
        # )
        # wax_loadings = pd.DataFrame(
        #     wax_pca.components_.T,
        #     columns=chain_lengths,
        #     index=wax_data.columns
        # )

        wax_ldings = wax_pca.components_
        wax_features = wax_data.columns
        wax_pc_values = np.arange(wax_pca.n_components_) + 1

        pca_dict = {
            "pca": wax_pca,
            "pc_values": wax_pc_values,
            "features": wax_features,
            "loadings": wax_ldings
        }

        for i in range(0, len(chain_lengths)):

            wax_pc = wax_pca.fit_transform(wax_data_scaled)[:,i]
            wax_scale_pc = 1.0 / (wax_pc.max() - wax_pc.min())
            wax_pc_score = wax_pc * wax_scale_pc

            pca_dict.update({f"wax_pc{i+1}": wax_pc})
            pca_dict.update({f"wax_scale_pc{i+1}": wax_scale_pc})
            pca_dict.update({f"wax_pc{i+1}_score": wax_pc_score})


        return pca_dict
