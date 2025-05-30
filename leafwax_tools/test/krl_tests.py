from leafwax_tools import WaxData
import pandas as pd
# import numpy as np

test_df = pd.read_excel("Hollister_et_al_2022_leafwax_data.xlsx")
test_data = WaxData(test_df)

total_conc = test_data.tot_conc(data_type="f", log=False)
total_logconc = test_data.tot_conc(data_type="f", log=True)
# print(total_conc)

# test_data2 = WaxData(3)

bad_df = pd.read_excel("bad_leafwax_data.xlsx")
bad_data = WaxData(bad_df)
# bad_conc = bad_data.tot_conc(data_type="f")