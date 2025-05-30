"""
The WaxData module is the main class for manipulating leaf wax data
imported as a pandas.core.frame.DataFrame
"""


import pandas as pd
import numpy as np
import warnings
from ..utils import validate_data


class WaxData:
    
    """
    Represents leaf wax data imported as a pandas DataFrame.
    
    Parameters
    ----------
    input_data : pandas.core.frame.DataFrame
        Pandas DataFrame of leaf wax data
        
    Attributes
    ----------
    input_data : pandas.core.frame.DataFrame
        Pandas DataFrame of leaf wax data
        
    
    Examples
    --------
    
    .. jupyter-execute::
        
        import leafwax_tools as leafwax
        wax_data = leafwax.WaxData(input_data=wax_df)
        
    """
    
    
    def __init__(self, input_data):

        self.data = input_data

        if type(input_data) != pd.core.frame.DataFrame:
            raise TypeError("Expecting type: pandas.core.frame.DataFrame")


    def total_conc(self, data_type, conc_name="conc", log=False, ret_name="tot_conc"):
        """
        Calculates the total leaf wax concentration of each sample (DataFrame row)

        Parameters
        ----------
        data_type : str
            Leaf wax compound class to search for; FAMEs/n-alkanoic acids ("f"), n-alkanes ("a")
        conc_name : str, optional
            String within self.data column names denoting leaf wax chain-length concentrations. The default is "conc".
        log : bool, optional
            Whether or not to calulcate the log of total leaf wax concentration. The default is False.
        ret_name : str, optional
            Name of the returned Pandas Series 'total_conc'. The default is "tot_conc".

        Returns
        -------
        total_conc : pandas.core.series.Series
            Pandas Series of leaf wax total concentrations per sample.

        """
        
        validate_data(self.data, data_type)
        
        conc_data_type = data_type + conc_name
        conc_df = self.data.filter(regex=conc_data_type)
        conc_arr = np.array(conc_df)
        total_conc_arr = np.zeros(len(conc_arr[:,0]))
        
        # add check for if self.data.conc_ugg_plant exists
        # warn user if not exists
        
        for row in range(0, len(conc_arr[:,0])):
            
            total_conc_arr[row] = np.nansum(conc_arr[row,:])
                
            if total_conc_arr[row] == 0:
                total_conc_arr[row] = np.nan
                    
        if log is True:
            total_conc_arr = np.log(total_conc_arr)
        else:
            total_conc_arr = total_conc_arr
        
        total_conc = pd.Series(data=total_conc_arr, name=ret_name)
        
        return total_conc


    '''
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

    # def paq():

    # def iso_avg():

    # def iso_diff():
