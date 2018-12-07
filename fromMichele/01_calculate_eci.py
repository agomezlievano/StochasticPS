# Hidalgo-Hausmann classical ECI implementation (ignores the importers)

import sys
import numpy as np
import pandas as pd
from scipy.linalg import eig

M_ccp = pd.read_csv("sitc4_ccp_2013.csv", sep = "\t")
M_cp = M_ccp.groupby(by = ["exporter", "sitc"]).sum().reset_index() # Collapse the matrix over the importers

# Adds a column with the total value exported by the exporter
M_cp = M_cp.merge(M_cp.groupby(by = "exporter")[["dollars"]].sum(), left_on = "exporter", right_index = True, suffixes = ("", " exporter tot"))
# Adds a column with the total value exported of the particular sitc product
M_cp = M_cp.merge(M_cp.groupby(by = "sitc")[["dollars"]].sum(), left_on = "sitc", right_index = True, suffixes = ("", " sitc tot"))
# Calculates the revealed comparative advantage (RCA), and binarizes it -- filling True only for the exporters which exported significant quantities of the product (RCA > 1)
M_cp["rca"] = (M_cp["dollars"] / M_cp["dollars exporter tot"]) / (M_cp["dollars sitc tot"] / M_cp["dollars"].sum()) > 1
# Transforms the data table into a dense M_cp matrix
M_cp = pd.pivot_table(M_cp, index = "exporter", columns = "sitc", values = "rca").fillna(False)

ubiquity = M_cp.sum(axis = 0) # Sums the number of exporters exporting the product with RCA > 1
diversity = M_cp.sum(axis = 1) # Sums the number of products exported by the exporter with RCA > 1

Q = M_cp / ubiquity # Column normalized M_cp
R = (M_cp.T / diversity).T # Row normalized M_cp
Spp = np.dot(R.T, Q).astype(float) # Square product-product matrix
Scc = np.dot(Q, R.T).astype(float) # Square exporter-exporter matrix

### This block calculates the eigenvector of the square product-product matrix and takes the second largest. If it doesn't correlate with ubiquity, multiplies by -1 ###

Spp_eigen = eig(Spp, left = True)
idx = Spp_eigen[0].argsort() # Numpy returns eigenvectors in random order, so we need to sort them so that we are sure we're picking the second largest
pci = np.real(Spp_eigen[1][:,idx][:,-2])

if np.corrcoef(pci, ubiquity)[0,1] > 0:
   pci *= -1

######

### This block calculates the eigenvector of the square exporter-exporter matrix and takes the second largest. If it doesn't correlate with ubiquity, multiplies by -1 ###

Scc_eigen = eig(Scc, left = True)
idx = Scc_eigen[0].argsort()  # Numpy returns eigenvectors in random order, so we need to sort them so that we are sure we're picking the second largest
eci = np.real(Scc_eigen[1][:,idx][:,-2])

if np.corrcoef(eci, diversity)[0,1] < 0:
   eci *= -1

######

# Merges pci with M_cp, so that each product gets assigned to its pci score
pci = M_cp.T.reset_index().merge(pd.DataFrame(pci), left_index = True, right_index = True)[["sitc", 0]].rename(columns = {0: "pci"})
# Merges eci with M_cp, so that each exporter gets assigned to its eci score
eci = M_cp.reset_index().merge(pd.DataFrame(eci), left_index = True, right_index = True)[["exporter", 0]].rename(columns = {0: "eci"})

# pci and eci gets standardized with the same strategy described in the Atlas of Economic Complexity
pci["pci"] = (pci["pci"] - pci["pci"].mean()) / (pci["pci"].std())
eci["eci"] = (eci["eci"] - eci["eci"].mean()) / (eci["eci"].std())

# Write tables to file
pci.to_csv("sitc4_2013_pci.csv", index = False, sep = "\t")
eci.to_csv("sitc4_2013_eci.csv", index = False, sep = "\t")

