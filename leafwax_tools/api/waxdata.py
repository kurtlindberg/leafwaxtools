"""
The WaxData module is the main class for manipulating leaf wax data
imported as a pandas.core.frame.DataFrame
"""

import pandas as pd
import numpy as np
import warnings

class WaxData:

    
    def __init__(self, input_data):

        self.data = input_data

        if type(input_data) != pd.core.frame.DataFrame:
            raise TypeError("Expecting type: pandas.core.frame.DataFrame")
            
        f_data = self.data.filter(regex="f")
        if len(f_data.columns) == 0:
            warnings.warn("Warning: column names do not label FAMEs (n-alkanoic acids) compound class with an 'f'; i.e. c20_fconc")
            
        a_data = self.data.filter(regex="a")
        if len(a_data.columns) == 0:
            warnings.warn("Warning: column names do not label n-alkanes compound class with an 'a'\n; i.e. c20_aconc")

        # add checks for column names


    def tot_conc(self, data_type, conc_name="conc", log=False, ret_name="tot_conc"):
        
        conc_data_type = data_type + conc_name
        conc_df = self.data.filter(regex=conc_data_type)
        conc_arr = np.array(conc_df)
        tot_conc_arr = np.zeros(len(conc_arr[:,0]))
        
        for row in range(0, len(conc_arr[:,0])):
            
            tot_conc_arr[row] = np.nansum(conc_arr[row,:])
                
            if tot_conc_arr[row] == 0:
                tot_conc_arr[row] = np.nan
                    
        if log is True:
            tot_conc_arr = np.log(tot_conc_arr)
        else:
            tot_conc_arr = tot_conc_arr
        
        tot_conc = pd.Series(data=tot_conc_arr, name=ret_name)
        
        return tot_conc


    def rel_abd(self, data_type, conc_name="conc", start=20, end=30, all_chain=True):

        conc_data_type = data_type + conc_name

        wax_conc_all = self.data.filter(regex=conc_data_type).fillna(0)
        wax_conc = pd.DataFrame()
        chain_lengths = list(range(start, end+1))

        if all_chain is True:
            chain_lengths = chain_lengths

        else:
            if data_type == "f":
                chain_lengths = [num for num in chain_lengths if num % 2 == 0]
            elif data_type == "a":
                chain_lengths = [num for num in chain_lengths if num % 2 == 1]
            else:
                raise ValueError("data_type can only be 'f' or 'a'")

        ## Filter for carbon chain-length concentration data within start-end range
        for n in chain_lengths:
            wax_chain = pd.DataFrame(
                data=np.array(wax_conc_all.filter(items=["c"+str(n)+"_"+conc_data_type])),
                columns=["c"+str(n) + '_' + data_type + "abd"]
            )
            wax_conc = pd.concat([wax_conc, wax_chain], axis=1)

        wax_conc_arr = np.array(wax_conc)
        rel_abd_arr = np.zeros(np.shape(wax_conc_arr))

        for row in range(0, len(wax_conc_arr[:,0])):
            for col in range(0, len(wax_conc_arr[0,:])):

                rel_abd_arr[row,col] = wax_conc_arr[row,col]/np.sum(wax_conc_arr[row,:])

        rel_abd = pd.DataFrame(data=rel_abd_arr, columns = wax_conc.columns)

        return rel_abd

    '''
    def acl(self, data_type, conc_name="conc", start=20, end=30, ret_name="acl"):

        conc_data_type = data_type + conc_name
        wax_conc_all = self.data.filter(regex=conc_data_type).fillna(0)
        # wax_conc_all = wax_conc_all.fillna(0)
        chain_lengths = list(range(start, end+1))
        wax_conc = pd.DataFrame()

        ## Filter for carbon chain-length concentration data within start-end range
        for n in chain_lengths:
            chain=pd.DataFrame(
                data=np.array(wax_conc_all.filter(regex=str(n))),
                columns=[str(n)]
            )
            wax_conc = pd.concat([wax_conc, chain], axis=1)

        wax_conc_arr = np.array(wax_conc)
        acl_numer = np.zeros(len(wax_conc_arr[:,0]))
        acl_arr = np.zeros(len(wax_conc[:,0]))

        ## Calculate ACL for each row (sample)
        for row in range(0, len(wax_conc_arr[:,0])):
            for col in range(0, len(wax_conc_arr[0,:])):

                acl_numer[row] += wax_conc_arr[row,col] * chain_lengths[col]

                acl_arr[row] = acl_numer[row]/np.sum(wax_conc_arr[row,:])

        acl = pd.Series(data=acl_arr, name=ret_name)

        return acl
    '''

    # def cpi():

    # def tot_conc():

    # def paq():

    # def iso_avg():

    # def iso_diff():
