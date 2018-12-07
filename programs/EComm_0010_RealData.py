
# coding: utf-8

# In[1]:

get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
get_ipython().magic('matplotlib inline')

import os
import sys


# In[2]:

# Allow imports from parent directory 
# https://stackoverflow.com/questions/34478398/import-local-function-from-a-module-housed-in-another-directory-with-relative-im/35273613#35273613
#module_path = os.path.abspath(os.path.join(os.pardir))
module_path = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\programs\\'
if module_path not in sys.path:
    sys.path.append(module_path)
sys.path


# In[3]:

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

LETTERS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]


# In[27]:

import EComm_0001_complexities 

# Paths
path_fig = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\figures\\'
path_data = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\data\\'
path_outputdata = 'C:\\Users\\agomez\\Dropbox\\Harvard\\LittleProjects\\StochasticPS\\outputdata\\'

# format of figures
figformat = "pdf"
save2file = True


# # Loading SITC 2015 data 
# from: https://intl-atlas-downloads.s3.amazonaws.com/index.html

# In[5]:

ccpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CCPY/S2_final_2015.dta"
cpy_filepath = "https://intl-atlas-downloads.s3.amazonaws.com/CPY/S2_final_cpy_all.dta"
ctyregions_filepath = "https://raw.githubusercontent.com/lukes/ISO-3166-Countries-with-Regional-Codes/master/all/all.csv"


# In[6]:

rowvarstring = 'exporter'
colvarstring = 'commoditycode'
valvarstring = 'mcp'
year = 2015


# In[7]:

ctyregs_df = pd.read_csv(ctyregions_filepath)
ctyregs_df.head()


# In[8]:

ctyregs_df = ctyregs_df[~ctyregs_df['sub-region'].isnull()]
ctyregs_df.shape


# In[9]:

###### LOADING THE DATA #######
longdf = pd.read_stata(cpy_filepath)


# In[10]:

# subsetting to the year
longdf = longdf[longdf['year']==year]
longdf.shape


# ## Getting the codes and checking consistency

# In[11]:

exp_codes = np.sort(longdf['exporter'].unique())
prod_codes = np.sort(longdf['commoditycode'].unique())
iso_codes = np.sort(ctyregs_df['alpha-3'].unique())

cty_codes = np.array(list(set(exp_codes) & set(iso_codes)))
print(len(cty_codes))

longdf = longdf[longdf[rowvarstring].isin(cty_codes)]
print(longdf.shape)


# ## Converting long to wide

# In[35]:

# ---------------------
# From long to wide format
Mcp_widedf = longdf[longdf['year']==year][[rowvarstring, colvarstring, valvarstring]].pivot(index=rowvarstring, 
                                                                       columns=colvarstring, 
                                                                       values=valvarstring).fillna(0.0)


# In[ ]:




# # -------------------------------------------------------------------------------------

# # EXPLORING THE ECI WITH THE REAL DATA

# In[13]:

kcomm = 5
mycolors = sns.color_palette("Set1", n_colors=kcomm, desat=.5)


# In[36]:

# Getting rid of the rows/columns of pure zeros
Mcp_widedf = Mcp_widedf.loc[:, (Mcp_widedf != 0).any(axis=0)]
Mcp_widedf = Mcp_widedf.loc[(Mcp_widedf != 0).any(axis=1), :]
Mcp_mat = np.array(Mcp_widedf)
print(Mcp_widedf.shape)


# In[28]:

nmat, ncP, npP = EComm_0001_complexities.ReorderingMatrix(Mcp_widedf)

fig = plt.figure(figsize=(8,5))

ax1 = fig.add_subplot(111)
ax1.spy(Mcp_widedf, aspect='auto')
ax1.set_xlabel('Products')
ax1.set_ylabel('Countries')

#ax2 = fig.add_subplot(122)
#ax2.spy(nmat, aspect='auto')
#ax2.set_xlabel('Products')
#ax2.set_ylabel('Countries')

plt.show()

#save2file=True
if(save2file):
    fig.savefig(path_fig + "EComm_0010_real_Mcp_matrix.{ff}".format(ff=figformat), bbox_inches='tight')


# In[16]:

# Calculating the c2c and p2p matrices, eigenvalues and left-eigenvectors
(Mc2c, Dc, Vc), (Mp2p, Dp, Vp) = EComm_0001_complexities.ECeigenvecs(Mcp_widedf)
minsize = min(Mcp_widedf.shape)
print(minsize)

# Calculating the right-eigenvectors
rightDc, rightVc = np.linalg.eig(Mc2c)


# In[29]:

fig = plt.figure(figsize=(5,5))

ax1 = fig.add_subplot(111)
ax1.imshow(Mc2c, aspect='auto', interpolation='nearest')
ax1.set_xlabel('Countries')
ax1.set_ylabel('Countries')

#ax2 = fig.add_subplot(122)
#ax2.spy(nmat, aspect='auto')
#ax2.set_xlabel('Products')
#ax2.set_ylabel('Countries')

plt.show()

#save2file=True
if(save2file):
    fig.savefig(path_fig + "EComm_0010_real_C_matrix.{ff}".format(ff=figformat), bbox_inches='tight')


# In[18]:

communitycolumn = 'region'

# left-eigenvalue data frame
Vc_df = pd.DataFrame(Vc, index=Mc2c.index)
Vc_df['Community'] = [ctyregs_df[ctyregs_df['alpha-3']==cname][communitycolumn].values[0] for cname in Mc2c.index.values]

# right-eigenvalue data frame
rightVc_df = pd.DataFrame(rightVc, index=Mc2c.index)
rightVc_df['Community'] = [ctyregs_df[ctyregs_df['alpha-3']==cname][communitycolumn].values[0] for cname in Mc2c.index.values]

print(Vc_df.head())
print(rightVc_df.head())


# In[30]:

fig = plt.figure(figsize=(10,4))
plt.subplots_adjust(wspace=0.2)

ax1 = fig.add_subplot(1,2,1)
# histogram of country eigenvalues
sns.distplot(Dc, bins=50, rug=True, kde=False, norm_hist=True, ax=ax1, axlabel="Eigenvalues")
sns.kdeplot(Dc, bw=.005, ax=ax1)
ax1.set_ylabel("Statistical Frequency")

ax2 = fig.add_subplot(1,2,2)
# histogram of country ECI's
sns.distplot(Vc[:,1], bins=40, rug=True, kde=False, norm_hist=True, ax=ax2, axlabel="2nd left-eigenvector (aka, ECI)")
sns.kdeplot(Vc[:,1], bw=.01, ax=ax2)
plt.show()

#save2file=True
if(save2file):
    fig.savefig(path_fig + "EComm_0010_real_Frequencies.{ff}".format(ff=figformat), bbox_inches='tight')


# In[31]:

fig = plt.figure(figsize=(10,4))
plt.subplots_adjust(wspace=0.2)

ax1 = fig.add_subplot(1,2,1)
# histogram of country eigenvalues
ax1.scatter(Dc.real, Dc.imag, marker="+", c="red", s=50)
#ax1.kdeplot(Dc, bw=.005, ax=ax1)
#ax1.set_ylabel("Statistical Frequency")

ax2 = fig.add_subplot(1,2,2)
# histogram of country ECI's
sns.distplot(Vc[:,1], bins=40, rug=True, kde=False, norm_hist=True, ax=ax2, axlabel="2nd left-eigenvector (aka, ECI)")
sns.kdeplot(Vc[:,1], bw=.01, ax=ax2)
plt.show()

#save2file=True
if(save2file):
    fig.savefig(path_fig + "EComm_0010_real_dots_and_Frequencies.{ff}".format(ff=figformat), bbox_inches='tight')


# In[21]:

print("The number of clusters is: {}".format(np.sum(Dc>0.25)))


# In[22]:

realcomms = np.unique(Vc_df['Community'].values)
kcomm = len(realcomms)
mycolors = sns.color_palette("Set1", n_colors=kcomm, desat=.5)

communities_vec = realcomms
numcommunities = len(communities_vec)
rncomm = np.arange(numcommunities)
mycolors = sns.color_palette("Set1", n_colors=numcommunities, desat=.5)

(communities_vec, numcommunities, rncomm)


# In[23]:

cty_marker_sizes = (Mcp_widedf.sum(axis=1).values/(0.8*np.min(Mcp_widedf.sum(axis=1).values)))**1.1
(np.min(cty_marker_sizes), np.max(cty_marker_sizes))


# In[32]:

fig = plt.figure(figsize=(15,10))

#ax1 = fig.add_subplot(131)
#ax2 = fig.add_subplot(132)
#ax3 = fig.add_subplot(133)


#########################################################################
# FIRST ROW: LEFT-EIGENVECTORS
ax1 = fig.add_subplot(2,3,1)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax1.scatter(Vc_df[Vc_df.Community==target_name][1], Vc_df[Vc_df.Community==target_name][2],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax1.set_xlabel('2nd left-eigenvector (aka, ECI)', fontsize=16)
ax1.set_ylabel('3rd left-eigenvector', fontsize=16)
#plt.legend(loc="best", shadow=False, scatterpoints=1)

ax2 = fig.add_subplot(2,3,2)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax2.scatter(Vc_df[Vc_df.Community==target_name][1], Vc_df[Vc_df.Community==target_name][3],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax2.set_xlabel('2nd left-eigenvector (aka, ECI)', fontsize=16)
ax2.set_ylabel('4th left-eigenvector', fontsize=16)
#plt.legend(loc="upper center", shadow=False, scatterpoints=1, bbox_to_anchor=(0.5, -0.2), ncol=5)

ax3 = fig.add_subplot(2,3,3)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax3.scatter(Vc_df[Vc_df.Community==target_name][2], Vc_df[Vc_df.Community==target_name][3],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax3.set_xlabel('3rd left-eigenvector', fontsize=16)
ax3.set_ylabel('4th left-eigenvector', fontsize=16)
#plt.legend(loc="center left", shadow=False, scatterpoints=1, bbox_to_anchor=(1.1, 0.5))



#########################################################################
# SECOND ROW: RIGHT-EIGENVECTORS
ax4 = fig.add_subplot(2,3,4)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax4.scatter(-rightVc_df[rightVc_df.Community==target_name][1], -rightVc_df[rightVc_df.Community==target_name][2],
                c=color, label=target_name, 
                s=cty_marker_sizes[rightVc_df.Community==target_name])
ax4.set_xlabel('2nd right-eigenvector', fontsize=16)
ax4.set_ylabel('3rd right-eigenvector', fontsize=16)
#plt.legend(loc="best", shadow=False, scatterpoints=1)

ax5 = fig.add_subplot(2,3,5)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax5.scatter(-rightVc_df[rightVc_df.Community==target_name][1], rightVc_df[rightVc_df.Community==target_name][3],
                c=color, label=target_name, 
                s=cty_marker_sizes[rightVc_df.Community==target_name])
ax5.set_xlabel('2nd right-eigenvector', fontsize=16)
ax5.set_ylabel('4th right-eigenvector', fontsize=16)

# LEGEND
ax5legend = ax5.legend(loc="upper center", shadow=False, scatterpoints=1, bbox_to_anchor=(0.5, -0.2), ncol=6,
          title=r'$\bf{Communities}$', fontsize=16, frameon=True, fancybox=True, markerscale=1)
plt.setp(ax5legend.get_title(),fontsize=16)

ax6 = fig.add_subplot(2,3,6)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax6.scatter(-rightVc_df[rightVc_df.Community==target_name][2], rightVc_df[rightVc_df.Community==target_name][3],
                c=color, label=target_name, 
                s=cty_marker_sizes[rightVc_df.Community==target_name])
ax6.set_xlabel('3rd right-eigenvector', fontsize=16)
ax6.set_ylabel('4th right-eigenvector', fontsize=16)





plt.subplots_adjust(wspace=0.4, hspace=0.4)

#plt.axis([-2, 3, -3, 3])
plt.show()

#save2file=True
if(save2file):
    fig.savefig(path_fig + "EComm_0010_real_eigenvectors.{ff}".format(ff=figformat), bbox_inches='tight')


# ### Principal Component Analysis on the Mcp

# In[33]:

pca = PCA(n_components = 5, whiten = True)
X_pca = pca.fit_transform(Mcp_widedf)


# In[34]:

fig = plt.figure(figsize=(15,5))

ax1 = fig.add_subplot(1,3,1)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax1.scatter(X_pca[Vc_df.Community==target_name, 0], X_pca[Vc_df.Community==target_name, 1],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax1.set_xlabel('PC0', fontsize=16)
ax1.set_ylabel('PC1', fontsize=16)
#plt.legend(loc="best", shadow=False, scatterpoints=1)

ax2 = fig.add_subplot(1,3,2)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax2.scatter(X_pca[Vc_df.Community==target_name, 0], X_pca[Vc_df.Community==target_name, 2],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax2.set_xlabel('PC0', fontsize=16)
ax2.set_ylabel('PC2', fontsize=16)
plt.legend(loc="upper center", shadow=False, scatterpoints=1, bbox_to_anchor=(0.5, -0.2), ncol=5)

ax3 = fig.add_subplot(1,3,3)
for color, i, target_name in zip(mycolors, rncomm, np.array(communities_vec)):
    ax3.scatter(X_pca[Vc_df.Community==target_name, 1], X_pca[Vc_df.Community==target_name, 2],
                c=color, label=target_name, 
                s=cty_marker_sizes[Vc_df.Community==target_name])
ax3.set_xlabel('PC1', fontsize=16)
ax3.set_ylabel('PC2', fontsize=16)
#plt.legend(loc="center", shadow=False, scatterpoints=1, bbox_to_anchor=(0.5, -0.05), ncol=5)


#plt.axis([-2, 3, -3, 3])
plt.show()


# In[ ]:



