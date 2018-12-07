import sys
import numpy as np
import pandas as pd
from scipy.linalg import eig
from statsmodels.regression.linear_model import OLS

M_ccp = pd.read_csv("sitc4_ccp_2013.csv", sep = "\t")

# Adds a column with the total value exported by the exporter
M_ccp = M_ccp.merge(M_ccp.groupby(by = "exporter")[["dollars"]].sum(), left_on = "exporter", right_index = True, suffixes = ("", " exporter tot"))
# Adds a column with the total value exported to the importer
M_ccp = M_ccp.merge(M_ccp.groupby(by = "importer")[["dollars"]].sum(), left_on = "importer", right_index = True, suffixes = ("", " importer tot"))
# Adds a column with the total value exported of the particular sitc product
M_ccp = M_ccp.merge(M_ccp.groupby(by = "sitc")[["dollars"]].sum(), left_on = "sitc", right_index = True, suffixes = ("", " sitc tot"))

# We need a constant for the regression
M_ccp["constant"] = 1
# We want to log all values so that we can exploit log's nice numerical properties
M_ccp[["dollars", "dollars exporter tot", "dollars importer tot", "dollars sitc tot"]] = np.log(M_ccp[["dollars", "dollars exporter tot", "dollars importer tot", "dollars sitc tot"]])

# The OLS creates an expectation of how much an exporter should export a product to an importer. This is the denominator of the RCA
rca_regression = OLS(M_ccp["dollars"], exog = M_ccp[["dollars exporter tot", "dollars importer tot", "dollars sitc tot", "constant"]], hasconst = True).fit()
# If we exponentiate the residuals of the regression, we have an idea of how much more than expected the exporter-importer-product connection was. This is logically equivalent to RCA.
M_ccp["rca"] = np.exp(rca_regression.resid) > 1



### This block calculates the exporter and the importer-product complexities ###

M_e_ip = pd.pivot_table(M_ccp, index = "exporter", columns = ("importer", "sitc"), values = "rca").fillna(False)

ubiquity = M_e_ip.sum(axis = 0) # Sums the number of exporters exporting the product with RCA > 1
diversity = M_e_ip.sum(axis = 1) # Sums the number of products exported by the exporter with RCA > 1
ubiquity[ubiquity == 0] = 1

Q = M_e_ip / ubiquity # Column normalized M_cp
R = (M_e_ip.T / diversity).T # Row normalized M_cp
#S_ip_ip = np.dot(R.T, Q).astype(float) # Square product-product matrix
S_e_e = np.dot(Q.values, R.T.values).astype(float) # Square exporter-exporter matrix

S_e_e_eigen = eig(S_e_e, left = True)
idx = S_e_e_eigen[0].argsort()  # Numpy returns eigenvectors in random order, so we need to sort them so that we are sure we're picking the second largest
eci = np.real(S_e_e_eigen[1][:,idx][:,-2])

if np.corrcoef(eci, diversity)[0,1] < 0:
   eci *= -1
######


eci = M_e_ip.reset_index().merge(pd.DataFrame(eci), left_index = True, right_index = True).iloc[:,[0,-1]]
eci.columns = ("exporter", "eci")
eci["eci"] = (eci["eci"] - eci["eci"].mean()) / (eci["eci"].std())
eci.to_csv("sitc4_2013_meci.csv", index = False, sep = "\t")
