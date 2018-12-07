import numpy as np
import pandas as pd
from scipy.stats import linregress

def calculate_growth(row, gdp):
   if row["year"] < 2006:
      timeline = gdp[(gdp["exporter"] == row["exporter"]) & (gdp["year"] > row["year"]) & (gdp["year"] <= row["year"] + 10)]["gdppc"]
      #return linregress(range(timeline.shape[0]), np.log(timeline.values)).slope
      return timeline.values[-1] / timeline.values[0]
   else:
      return None

gdp = pd.read_csv("gdp_percapita_2010const.csv").drop("Country Name", 1).rename(columns = {"Country Code": "exporter"}).set_index("exporter").unstack().reset_index().dropna()
gdp.columns = ("year", "exporter", "gdppc")
gdp["year"] = gdp["year"].astype(int)
gdp["gdppc"] = gdp["gdppc"].astype(float)

natural_resources = pd.read_csv("natural_resources.csv").drop("Country Name", 1).rename(columns = {"Country Code": "exporter"}).set_index("exporter")
natural_resources = natural_resources.mean(axis = 1).reset_index().fillna(0)
drop_countries = set(natural_resources[natural_resources[0] > 8]["exporter"])

gdp = gdp[~gdp["exporter"].isin(drop_countries)]
gdp["growth"] = gdp.apply(calculate_growth, axis = 1, args = (gdp,))

complexities = pd.DataFrame()
for year in range(1962, 2016):
   complexity = pd.read_csv("complexities/%s.csv" % year, sep = "\t")
   complexity["year"] = year
   complexities = pd.concat([complexities, complexity])

table = complexities.merge(gdp.dropna(), on = ("exporter", "year"))
table.to_csv("eci_growth_regression_table.csv", sep = "\t", index = False)
