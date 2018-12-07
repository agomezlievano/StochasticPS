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

# We need to make sure that our data table will have all pairs of exporter-importer, so we generate them manually
complete_index = pd.MultiIndex.from_tuples([(e,i) for i in np.sort(M_ccp["importer"].unique()) for e in np.sort(M_ccp["exporter"].unique())], names = ("exporter", "importer"))
# This generates a table exporter+importer X products table
M_ccp = pd.pivot_table(M_ccp, index = ["exporter", "importer"], columns = "sitc", values = "rca").reindex(index = complete_index).fillna(False)
# This makes M_ccp a 3-dimensional matrix product x importer x exporter
M_ccp = M_ccp.unstack().unstack().values.reshape((M_ccp.columns.shape[0], M_ccp.index.levels[1].shape[0], M_ccp.index.levels[0].shape[0]))

kp = []
ki = []
ke = []

kp.append(M_ccp.sum(axis = 2).sum(axis = 1).astype(float)) # How many exporter-importer pairs involve the product?
ki.append(M_ccp.sum(axis = 0).sum(axis = 1).astype(float)) # How many exporter-product pairs involve the importer?
ke.append(M_ccp.sum(axis = 1).sum(axis = 0).astype(float)) # How many product-importer pairs involve the exporter?

for n in range(1, 19):
   kp.append((M_ccp * (ki[n - 1][:,None] * ke[n - 1][:,None].T)).sum(axis = 2).sum(axis = 1)) # Misses the normalizing term
   ki.append((np.swapaxes(M_ccp, 0, 1) * (kp[n - 1][:,None] * ke[n - 1][:,None].T)).sum(axis = 2).sum(axis = 1)) # Misses the normalizing term
   ke.append((np.swapaxes(np.swapaxes(M_ccp, 0, 1), 0, 2) * (kp[n - 1][:,None] * ki[n - 1][:,None].T)).sum(axis = 2).sum(axis = 1)) # Misses the normalizing term

