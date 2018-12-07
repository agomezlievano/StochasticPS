import os
import sys


# ===============================================================================================
# Allow imports from parent directory 
# https://stackoverflow.com/questions/34478398/import-local-function-from-a-module-housed-in-another-directory-with-relative-im/35273613#35273613
#module_path = os.path.abspath(os.path.join(os.pardir))
module_path = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\programs\\'
if module_path not in sys.path:
    sys.path.append(module_path)
sys.path

# ===============================================================================================
# libraries
import time
import datetime

import pandas as pd
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

import seaborn as sns

import itertools
import collections
import warnings
import IPython.display
import scipy.stats
import networkx as nx
from operator import itemgetter

plt.style.use('seaborn-white')
plt.rc('font', family='serif', serif='Helvetica')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=8)
plt.rc('ytick', labelsize=8)
plt.rc('axes', labelsize=16, linewidth=0.5)


from sklearn.preprocessing import StandardScaler, RobustScaler, FunctionTransformer, PolynomialFeatures
from sklearn.decomposition import PCA, NMF
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error

from numpy.random import exponential, negative_binomial, randint, choice, binomial
from random import shuffle

import statsmodels.formula.api as smf
from statsmodels.regression.linear_model import OLS
from statsmodels.iolib.summary2 import summary_col

from scipy.spatial import distance

LETTERS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]


import EComm_0001_complexities 

# ===============================================================================================
# Paths
path_fig = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\figures\\'
path_data = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\data\\'
path_outputdata = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\outputdata\\'
path_inputdata = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\inputdata\\'

# format of figures
figformat = "pdf"
save2file = False

#ccpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CCPY/S2_final_{yr}.dta"
ccpy_filepath = path_inputdata + "S2_final_{yr}.csv"

#cpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CPY/S2_final_cpy_all.dta"
cpy_filepath = path_inputdata + "S2_final_cpy_all.dta"

#ctyregions_filepath = "https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv"
ctyregions_filepath = path_data + "all.csv"

gdp_filepath = path_data + "WorldBank_GDPperCapita_1962_to_2015.json"

# ===============================================================================================
# Some setup
colimporter = 'importer'
colprod = 'commoditycode'
newcol = 'imp_prod'

rowvarstring = 'exporter'
colvarstring = 'commoditycode'
valvarstring = 'export_value'

rcaval = 'mcp'

CountriesToRemove = ["GRL", "SPM", "SOM", "BLZ", "AFG", "NCL", "ATA", "SMR",
"PCN", "SWZ", "RWA", "KIR", "GNB", "MLT", "IMN", "MSR", "GRD", "SUR", "BDI", "NFK",
"BHR", "VUT", "MHL", "CPV", "TON", "HTI", "MYT", "MDV", "GUF", "CXR", "MCO", "LSO",
"ESH", "BRN", "AIA", "IRQ", "KNA", "ANT", "AND", "VCT", "DJI", "SGS", "GNQ", "ATG",
"MWI", "TKL", "ERI", "LIE", "MNG", "VGB", "ABW", "NPL", "COK", "GMB", "SHN", "TLS",
"BFA", "MAF", "BVT", "NIU", "MNP", "YEM", "FRO", "CYP", "FJI", "SLE", "WLF", "WSM",
"BRB", "SLB", "ASM", "HMD", "DMA", "GUM", "STP", "GIB", "CUW", "BMU", "TCD", "SYC",
"PSE", "MNE", "NRU", "TCA", "LUX", "BHS", "FSM", "CAF", "FLK", "ISL", "BTN", "PLW",
"BEN", "VAT", "TUV", "COM", "VIR", "SXM", "GUY", "MUS", "LCA", "NER", "MTQ", "CYM",
"CCK", "UMI", "PYF"]


# ===============================================================================================
###### LOADING THE DATA #######
ctyregs_df = pd.read_csv(ctyregions_filepath)

ctyregs_df = ctyregs_df[((~ctyregs_df['sub-region'].isnull()) & (~ctyregs_df['region'].isnull()))]
ctyregs_df.shape

longdf_cp = pd.read_stata(cpy_filepath)


# ===============================================================================================
exp_codes = set(longdf_cp['exporter'])
iso_codes = set(ctyregs_df['alpha-3'])
cty_codes = set(exp_codes) & set(iso_codes)

# Subsetting and cleaning
longdf_cp = longdf_cp[longdf_cp["exporter"].isin(cty_codes) & 
                      (longdf_cp["export_value"] > 0) & 
                      (longdf_cp["commoditycode"] != ".")].drop(["inatlas",
                                                                 "oppval",
                                                                 "oppgain",
                                                                 "distance",
                                                                 "import_value"], 1)




# ===============================================================================================
# products that are in all years
pcodes2keep = longdf_cp[longdf_cp['year']=='2015'].commoditycode.unique()

# countries that are in all years
ccodes2keep = longdf_cp[longdf_cp['year']=='2015'].exporter.unique()

for yr in longdf_cp.year.unique():
    pcodesyr = longdf_cp[longdf_cp['year']==yr].commoditycode.unique()
    pcodes2keep = np.array(list(set(pcodes2keep) & set(pcodesyr)))

    ccodesyr = longdf_cp[longdf_cp['year']==yr].exporter.unique()
    ccodes2keep = np.array(list(set(ccodes2keep) & set(ccodesyr)))

#longdf_cp = longdf_cp[list(map(lambda x: not(x in CountriesToRemove), longdf_cp['exporter']))]
#longdf_cp = longdf_cp[list(map(lambda x: x in pcodes2keep, longdf_cp['commoditycode']))]
longdf_cp = longdf_cp[~longdf_cp["exporter"].isin(CountriesToRemove) & 
                      longdf_cp["exporter"].isin(ccodes2keep) & 
                      longdf_cp["commoditycode"].isin(pcodes2keep)]
#print((len(longdf_cp.year.unique()), len(longdf_cp.exporter.unique()), len(longdf_cp.commoditycode.unique())))


# ===============================================================================================
# Creating totals
longdf_cp = longdf_cp.merge(longdf_cp.groupby(by = ["year", "exporter"])[["export_value"]].sum(), left_on = ["year", "exporter"], right_index = True, suffixes = ("", " exporter tot"))
longdf_cp = longdf_cp.merge(longdf_cp.groupby(by = ["year", "commoditycode"])[["export_value"]].sum(), left_on = ["year", "commoditycode"], right_index = True, suffixes = ("", " commoditycode tot"))

# Running regression
longdf_cp["constant"] = 1
longdf_cp[["log_export_value", "log_export_value exporter tot", "log_export_value commoditycode tot"]] = np.log(longdf_cp[["export_value", "export_value exporter tot", "export_value commoditycode tot"]])

rca_regression = OLS(longdf_cp["log_export_value"], exog = longdf_cp[["log_export_value exporter tot", "log_export_value commoditycode tot", "constant"]], hasconst = True).fit()
#longdf_cp[rcaval] = (rca_regression.resid > 0).astype(int)  
longdf_cp["mcpfromOLS"] = (rca_regression.resid > 0).astype(int)  


# ===============================================================================================
ctryskept = longdf_cp.exporter.unique()
prodskept = longdf_cp.commoditycode.unique()
yearskept = longdf_cp.year.unique()
print((len(yearskept), len(ctryskept), len(prodskept)))


# ===============================================================================================
# Looping through the years
## 1) Load the data for the year
## 2) Create the Mcp and Mcpi
## 3) Compute the Ccc_cp and Ccc_cpi
## 4) Compute the left and right eigenvectors
## 5) Compute the distance to the center
## 6) Append distance to a long list of countries and years.

dists_allyears_df = longdf_cp[['year', 'exporter', 'population', 'eci']].drop_duplicates(subset=['year', 'exporter'], keep='last')
dists_allyears_df['leftDistance2Origin_cp'] = 0
dists_allyears_df['leftDistance2Origin_ccp'] = 0

dists_allyears_df['leftSimDirection2MostComplex_cp'] = 0
dists_allyears_df['leftSimDirection2MostComplex_ccp'] = 0

dists_allyears_df['leftDistance2MostComplex_cp'] = 0
dists_allyears_df['leftDistance2MostComplex_ccp'] = 0

dists_allyears_df['leftEntropyCommunities_cp'] = 0
dists_allyears_df['leftEntropyCommunities_ccp'] = 0

dists_allyears_df['rightDistance2Origin_cp'] = 0
dists_allyears_df['rightDistance2Origin_ccp'] = 0

dists_allyears_df['rightSimDirection2MostComplex_cp'] = 0
dists_allyears_df['rightSimDirection2MostComplex_ccp'] = 0

dists_allyears_df['rightDistance2MostComplex_cp'] = 0
dists_allyears_df['rightDistance2MostComplex_ccp'] = 0

dists_allyears_df['rightEntropyCommunities_cp'] = 0
dists_allyears_df['rightEntropyCommunities_ccp'] = 0


# ===============================================================================================
years = np.arange(1962, 2016)

for year in years:
    #year = 1970
    print("Starting the year {}".format(year))

    # 1) Load the data cpi for the year
    longdf_ccpy = pd.read_csv(ccpy_filepath.format(yr=year))
    longdf_cpy = longdf_cp[longdf_cp['year']==str(year)][[rowvarstring, colvarstring, rcaval]]

    # zero-padding numeric commodity codes, some modifications
    longdf_ccpy[colprod] = longdf_ccpy[colprod].apply(lambda x: str(x).zfill(4))
    longdf_ccpy = longdf_ccpy[~longdf_ccpy[rowvarstring].isnull() & ~longdf_ccpy["importer"].isnull()] # removing NAs
    longdf_ccpy = longdf_ccpy[longdf_ccpy[rowvarstring].isin(ctryskept) & 
                          longdf_ccpy["importer"].isin(ctryskept) &
                          longdf_ccpy[colprod].isin(prodskept) & 
                  (longdf_ccpy[valvarstring] > 0) & 
                  (longdf_ccpy[colprod] != ".")]



    # compute discrete rca
    longdf_ccpy = EComm_0001_complexities.RCA_from_longOLS_insingleyear(longdf_ccpy,
                                                                      dim_vars=[rowvarstring, colimporter, colprod],
                                                                      value_var=valvarstring, 
                                                                      new_var=rcaval,
                                                                      discrete=True)
    longdf_ccpy[newcol] = longdf_ccpy[colimporter]+longdf_ccpy[colprod] 

    # 2) Create the Mcp and Mcpi
    print("Ready to start to create the mcps")
    Mccp_widedf = longdf_ccpy[[rowvarstring, newcol, rcaval]].pivot(index=rowvarstring, 
                                                                   columns=newcol, 
                                                                   values=rcaval).fillna(0.0)    
    print("The shape of the Mccp is {}".format(Mccp_widedf.shape))
    # From long to wide format
    Mcp_widedf = longdf_cpy[[rowvarstring, colvarstring, rcaval]].pivot(index=rowvarstring, 
                                                                       columns=colvarstring, 
                                                                       values=rcaval).fillna(0.0)
    print("The shape of the Mcp is {}".format(Mcp_widedf.shape))

    # 3) Compute the Ccc_cp and Ccc_cpi and 4) Compute the left and right eigenvectors
    (Mc2c_ccp, Dc_ccp, leftVc_ccp, rightVc_ccp) = EComm_0001_complexities.ECeigenvecs(Mccp_widedf)
    (Mc2c_cp, Dc_cp, leftVc_cp, rightVc_cp) = EComm_0001_complexities.ECeigenvecs(Mcp_widedf)
        
    # 5) Compute the complexity measures
    print("Ready to compute the distances")
    kccp = 4
    kcp  = 3
    # 5.1) distance to the origin in the right-eigenspace
    ddf_ccp = EComm_0001_complexities.distance_to_center(Mccp_widedf, kcomm=kccp, Mc2c=Mc2c_ccp, Dc=Dc_ccp, leftVc=leftVc_ccp, rightVc=rightVc_ccp)
    ddf_cp = EComm_0001_complexities.distance_to_center(Mcp_widedf, kcomm=kcp, Mc2c=Mc2c_cp, Dc=Dc_cp, leftVc=leftVc_cp, rightVc=rightVc_cp)
    mostcomplexix_ccp = ddf_ccp.sort_values(by=["rightDist2Origin"], ascending=False).index[0]
    mostcomplexix_cp = ddf_cp.sort_values(by=["rightDist2Origin"], ascending=False).index[0]

    # 5.2) CosineSimilarity to most distant (1) country in left-eigenspace
    sims_ccp = EComm_0001_complexities.cosine_angle_to_target(Mccp_widedf, mostcomplexix_ccp, kcomm=kccp, Mc2c=Mc2c_ccp, Dc=Dc_ccp, leftVc=leftVc_ccp, rightVc=rightVc_ccp)
    sims_cp = EComm_0001_complexities.cosine_angle_to_target(Mcp_widedf, mostcomplexix_cp, kcomm=kcp, Mc2c=Mc2c_cp, Dc=Dc_cp, leftVc=leftVc_cp, rightVc=rightVc_cp)

    # 5.3) distance to the most distant (1) country in the right-eigenspace
    d2t_ccp = EComm_0001_complexities.distance_to_target(Mccp_widedf, mostcomplexix_ccp, kcomm=kccp, Mc2c=Mc2c_ccp, Dc=Dc_ccp, leftVc=leftVc_ccp, rightVc=rightVc_ccp)
    d2t_cp = EComm_0001_complexities.distance_to_target(Mcp_widedf, mostcomplexix_cp, kcomm=kcp, Mc2c=Mc2c_cp, Dc=Dc_cp, leftVc=leftVc_cp, rightVc=rightVc_cp)


    # 6) Append distance to a long list of countries and years.
    print("Generating the dataframe")

    ddf_ccp['year'] = str(year)
    ddf_cp['year'] = str(year)
    sims_ccp['year'] = str(year)
    sims_cp['year'] = str(year)
    d2t_ccp['year'] = str(year)
    d2t_cp['year'] = str(year)

    # 7) Merging
    print("Merging")

    # First merging the ccp information
    ctysorderyear = dists_allyears_df[dists_allyears_df['year']==str(year)]['exporter'].values

    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftDistance2Origin_ccp'] = ddf_ccp.loc[ctysorderyear]['leftDist2Origin'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftDistance2Origin_cp'] = ddf_cp.loc[ctysorderyear]['leftDist2Origin'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftSimDirection2MostComplex_ccp'] = sims_ccp.loc[ctysorderyear]['leftAngle2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftSimDirection2MostComplex_cp'] = sims_cp.loc[ctysorderyear]['leftAngle2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftDistance2MostComplex_ccp'] = d2t_ccp.loc[ctysorderyear]['leftDistance2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'leftDistance2MostComplex_cp'] = d2t_cp.loc[ctysorderyear]['leftDistance2target'].values

    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightDistance2Origin_ccp'] = ddf_ccp.loc[ctysorderyear]['rightDist2Origin'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightDistance2Origin_cp'] = ddf_cp.loc[ctysorderyear]['rightDist2Origin'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightSimDirection2MostComplex_ccp'] = sims_ccp.loc[ctysorderyear]['rightAngle2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightSimDirection2MostComplex_cp'] = sims_cp.loc[ctysorderyear]['rightAngle2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightDistance2MostComplex_ccp'] = d2t_ccp.loc[ctysorderyear]['rightDistance2target'].values
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'rightDistance2MostComplex_cp'] = d2t_cp.loc[ctysorderyear]['rightDistance2target'].values


    # Attaching measure of the number of communities
    Dc_ccp[Dc_ccp==0] = 1.0
    Dc_cp[Dc_cp==0] = 1.0

    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'EntropyCommunities_ccp'] = -np.dot(Dc_ccp.real, np.log(Dc_ccp.real))
    dists_allyears_df.loc[dists_allyears_df['year']==str(year),'EntropyCommunities_cp'] = -np.dot(Dc_cp.real, np.log(Dc_cp.real))

    # take the values that of interest from the clashes

    print("end of cycle")
    print("")

    


# ===============================================================================================
save2file = True
if(save2file):
    # name of file
    folder = "~\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\outputdata\\"
    finalfilename = "EComm_0003_constructingdataallyears.csv"

    # To csv file
    dists_allyears_df.to_csv(folder + finalfilename, index=False, encoding='utf-8')


# ===============================================================================================



# ===============================================================================================







