import sys
import numpy as np
import pandas as pd
from scipy.linalg import eig
from sklearn.metrics.pairwise import euclidean_distances
from statsmodels.regression.linear_model import OLS

cpy = pd.read_stata("S2_final_cpy_all.dta")
countries_to_keep = cpy[(cpy["year"] == 2015)][["exporter", "population"]].drop_duplicates()
countries_to_keep = set(countries_to_keep[countries_to_keep["population"] >= 330815.0]["exporter"])
cpy = cpy[cpy["exporter"].isin(countries_to_keep)]

for year in range(1962, 2016):
   longdf = cpy[(cpy["year"] == year) & (cpy["export_value"] > 0) & (cpy["commoditycode"] != ".")][["exporter", "commoditycode", "export_value"]]
   longdf = longdf.merge(longdf.groupby(by = "exporter")[["export_value"]].sum(), left_on = "exporter", right_index = True, suffixes = ("", " exporter tot"))
   longdf = longdf.merge(longdf.groupby(by = "commoditycode")[["export_value"]].sum(), left_on = "commoditycode", right_index = True, suffixes = ("", " commoditycode tot"))
   longdf["constant"] = 1
   longdf[["export_value", "export_value exporter tot", "export_value commoditycode tot"]] = np.log(longdf[["export_value", "export_value exporter tot", "export_value commoditycode tot"]])
   rca_regression = OLS(longdf["export_value"], exog = longdf[["export_value exporter tot", "export_value commoditycode tot", "constant"]], hasconst = True).fit()
   longdf["rca"] = rca_regression.resid > 0
   M_cp = pd.pivot_table(longdf, index = "exporter", columns = "commoditycode", values = "rca").fillna(False)
   M_cp = M_cp.loc[:, (M_cp != 0).any(axis = 0)]
   M_cp = M_cp.loc[(M_cp != 0).any(axis = 1), :]
   ubiquity = M_cp.sum(axis = 0) # Sums the number of exporters exporting the product with RCA > 1
   diversity = M_cp.sum(axis = 1) # Sums the number of products exported by the exporter with RCA > 1
   ubiquity[ubiquity == 0] = 1
   Q = M_cp / ubiquity # Column normalized M_cp
   R = (M_cp.T / diversity).T # Row normalized M_cp
   Mc2c = np.dot(Q.values, R.T.values).astype(float) # Square exporter-exporter matrix
   Dc, Vc, rightVc = eig(Mc2c, left = True)
   idx = Dc.argsort()[::-1]  # Numpy returns eigenvectors in random order, so we need to sort them so that we are sure we're picking the second largest
   Dc = Dc[idx]
   Vc = Vc[:,idx]
   rightVc = rightVc[:,idx]
   ncomms = (1 - (Dc.real * range(1, Dc.shape[0] + 1)) >= -.1).sum()
   sys.stderr.write("%s (%s)...\n" % (year, ncomms))
   Fc = np.ones(M_cp.shape[0]).astype(float)
   Qp = np.ones(M_cp.shape[1]).astype(float)
   for i in xrange(20):
      Fcnew = M_cp.dot(Qp)
      Fcnew[Fcnew == 0.0] = 1.0
      denom = M_cp.T.dot(1.0 / Fc)
      denom[denom == 0.0] = 0.1
      Qpnew = 1.0 / denom
      Fc = Fcnew / np.mean(Fcnew)
      Qp = Qpnew / np.mean(Qpnew)
   complexities = pd.DataFrame(diversity, columns = ("diversity",))
   complexities["eci"] = Vc.T[1]
   complexities["fitness"] = Fc
   complexities["redfc fix"] = euclidean_distances([[0] * rightVc.shape[0],], Y = rightVc)[0]
   complexities["redfc var"] = euclidean_distances([[0] * ncomms,], Y = rightVc[:,:ncomms])[0]
   complexities.reset_index().to_csv("complexities_cp/%s.csv" % year, sep = "\t", index = False)
