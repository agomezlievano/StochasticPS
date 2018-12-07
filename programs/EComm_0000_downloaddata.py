
# EComm_0000_downloaddata.py

# #########################################################################
# Download some of the data
# #########################################################################


# =========================================================================
# LIBRARIES
# -------------------------------------------------------------------------
import pandas as pd
import numpy as np
import numpy.random as npr
from numpy.linalg import matrix_power as matpower

import numpy.linalg as linalg
from scipy.linalg import eig
from operator import itemgetter

from sklearn.preprocessing import StandardScaler, RobustScaler
from statsmodels.regression.linear_model import OLS

from numpy.random import exponential, negative_binomial, randint, choice, binomial
from random import shuffle

import os
import sys
from pathlib import Path

# =========================================================================
# PATHS
# -------------------------------------------------------------------------

# Allow imports from parent directory 
# https://stackoverflow.com/questions/34478398/import-local-function-from-a-module-housed-in-another-directory-with-relative-im/35273613#35273613
#module_path = os.path.abspath(os.path.join(os.pardir))
module_path = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\programs\\'
if module_path not in sys.path:
    sys.path.append(module_path)
sys.path


import EComm_0001_complexities 

# Paths
path_fig = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\figures\\'
path_data = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\data\\'
path_outputdata = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\outputdata\\'
path_inputdata = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\inputdata\\'

ccpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CCPY/"
ccpy_filename = "S2_final_{yr}.dta"
ccpy_filenameCSV = "S2_final_{yr}.csv"
cpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CPY/S2_final_cpy_all.dta"
ctyregions_filepath = "https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv"
gdp_filepath = "C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\data\\WorldBank_GDPperCapita_1962_to_2015.json"



# =========================================================================
# COUNTRY WITH REGIONS
# -------------------------------------------------------------------------
ctyregs_df = pd.read_csv(ctyregions_filepath)
filename = "ctyregions.csv"
ctyregs_df.to_csv(path_inputdata + filename)


# =========================================================================
# EXPORTER-IMPORTER DATA
# -------------------------------------------------------------------------

years = np.arange(1962, 2017)


for year in years:
    print("Starting the year {}".format(year))
    my_file = Path(path_inputdata + ccpy_filenameCSV.format(yr=year))
    if(my_file.is_file()):
        continue

    # 1) Load the data cpi for the year
    longdf_ccpy = pd.read_stata(ccpy_filepath + ccpy_filename.format(yr=year))

    # saving
    longdf_ccpy.to_csv(path_inputdata + ccpy_filenameCSV.format(yr=year))

    print("end of cycle")
    print("")

 

