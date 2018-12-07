# EComm_0001_complexities.py

# #########################################################################
# Functions to compute complexity indices
# #########################################################################


# =========================================================================
# LIBRARIES
# -------------------------------------------------------------------------
import pandas as pd
import numpy as np
import numpy.random as npr
from numpy.linalg import matrix_power as matpower

#import numpy.linalg as linalg
from scipy.linalg import eig
from operator import itemgetter

from sklearn.preprocessing import StandardScaler, RobustScaler
from statsmodels.regression.linear_model import OLS

from numpy.random import exponential, negative_binomial, randint, choice, binomial, lognormal
from random import shuffle

from scipy.spatial import distance


# =========================================================================
# MY FUNCTIONS
# -------------------------------------------------------------------------
def RCA_from_Xcp(mat_df, discrete=False):
    if(discrete==True):
        return(np.floor(mat_df.sum().sum()*(mat_df.T/mat_df.sum(axis=1)).T/mat_df.sum(axis=0)).astype(bool).astype(int))
    else:
        return(mat_df.sum().sum()*(mat_df.T/mat_df.sum(axis=1)).T/mat_df.sum(axis=0))


def RCA_from_longOLS_insingleyear(longdf_cp, dim_vars=["exporter", "commoditycode"], value_var="export_value", new_var="rca_ols", discrete=False):
    log_var_names = [None]*(len(dim_vars)+1)
    var_names = [None]*(len(dim_vars))
    for vari in range(len(dim_vars)):
        var = dim_vars[vari]
        longdf_cp = longdf_cp.merge(longdf_cp.groupby(by = var)[[value_var]].sum(), left_on = var, right_index = True, suffixes = ("", "_{}_tot".format(var)))
        log_var_names[vari] = "log_" + value_var + "_{}_tot".format(var)
        var_names[vari] = value_var + "_{}_tot".format(var)

    longdf_cp["constant"] = 1
    log_var_names[len(dim_vars)] = "constant"
    
    longdf_cp[log_var_names[:len(dim_vars)]] = np.log(longdf_cp[var_names])

    rca_regression = OLS(longdf_cp[log_var_names[0]], exog = longdf_cp[log_var_names[1:]], hasconst = True).fit()
    if(discrete==True):
        longdf_cp[new_var] = (rca_regression.resid > 0).astype(int)        
        return(longdf_cp[list(np.concatenate((dim_vars, [new_var])))])
    else:
        longdf_cp[new_var] = np.exp(rca_regression.resid)       
        return(longdf_cp[list(np.concatenate((dim_vars, [new_var])))])

        
#function equivalent to Position[vec, object] in Mathematica
def Position(vec, object):
    positions = [j for j,x in enumerate(vec) if np.all(x==object)]
    return(positions)

#function that fills a matrix with 1's with probability q
def randomchoice(q,Na,Np):
    #matrix Na \times Np
    mat = np.array([[ floor(npr.uniform()+q)+0.0 for i in range(Np)] for j in range(Na)])
    return mat


#Reordering of a matrix from the most diverse city to the least, and from the most ubiquitous
#product to the least.
def ReorderingMatrix(mat, cityPos=None, prodPos=None):
    #cityPos[i]=j: current j-th city will now be in position i.
    #prodPos[i]=j: current j-th capability will now be in position i.

    #new matrix
    nmat = np.array(mat)

    #new positions for cities and products
    ncP = cityPos
    npP = prodPos

    Nc = len(nmat)   # number of cities
    Np = len(nmat.T) # number of products

    if(cityPos==None):
        diversification = np.sum(nmat, axis=1)
        if(Nc!=len(diversification)): print("mat has wrong dimensions!")
        
        ncP = np.argsort(-diversification)

    if(prodPos==None):
        ubiquity = np.sum(nmat, axis=0)
        if(Np!=len(ubiquity)): print("mat has wrong dimensions!")
        
        npP = np.argsort(-ubiquity)

    nmatfinal = np.array([[ nmat[ncP[c]][npP[p]] for p in range(Np)] for c in range(Nc)])
    return(nmatfinal, ncP, npP)

	
def Replace(vector, value, replacement):
    vec = list(vector)
    for index, item in enumerate(vec):
        if(item==value):
            vec[index] = replacement
    return(np.array(vec))

# modified RCA
def modRCA_func(rcaval, specifym=False, mspecified=25):
    """
    The idea behind modRCA is to transform the RCA-type measure
    into another measure that is normally distributed, but 
    which respects the distribution of mass
    between 0 and 1, and from 1 and above, that RCA has. The idea is also
    to keep modRCA=0 if RCA=0. 
    Hence, the interpretations of RCA and modRCA are
    the same.

    Inputs
    ------
    rcaval : numpy array of non-negative floats.
        These are the original RCA values.

    specifym : Binary. (optional)
        If 'False', then the function will ignore the
        parameter 'mspecified', and will calculate
        it so that minimum positive value of RCA
        is equal to the minimum positive value of 
        modRCA.
        If 'True', the user must pass a value for
        'mspecified', which will be used to 
        displace and scale the log(RCA).

    mspecified : float. (optional)
        This values is valid of specifym is True, 
        otherwise, it is ignored.

    Returns
    -------
    modRCA: numpy array of non-negative floats.
        This is the 'modified' RCA.

    Examples
    --------
    # Simple example:
    modrca_vec = modRCA_func(rca_vec)

    """

    # Generate a copy
    modRCA = np.array(rcaval)

    # Take only positive values
    xp = modRCA[modRCA>0]

    # Transform the positive values
    if(specifym):
        modRCA[modRCA>0] = (np.log(xp) + mspecified)/mspecified
    else:
        xmin = np.min(xp)
        m = np.log(xmin)/(xmin - 1.0)
        modRCA[modRCA>0] = (np.log(xp) + m)/m
        

    return modRCA
# end def modRCA_func

def weightedproduct(vec, weights):
    #this function is equivalent to 'weightedgeommean', but apparently faster
    if(any(x==0 for x in vec)): 
        return(0.0)
    else:
        logave = np.average(np.log(vec), weights=weights)
        return(np.exp(logave))


#==============================================================================================
# Community discovery
# def laplacian_community(A, k):
    # # Finding 
    
        
        
#==============================================================================================
# Production function
def Leontief(cvec, pvec):
    if(len(cvec)!=len(pvec)): print("Number of capabilities is different in cities than in products!")

    #positions to index the capabilities of the product defined by pvec
    caps = [i for i,x in enumerate(pvec) if x>0]
        
    #I subset the capabilities of the city given by caps to get the factors of production available to produce product pvec.
    factorsofprod = np.array(cvec)[caps]

    #I take the minimum among the factors of production.
    return(min(factorsofprod))


def CobbDouglas(cvec, pvec):
    #positions to index the capabilities of the product defined by pvec
    caps = [i for i,x in enumerate(pvec) if x>0]

    #I subset the capabilities of the city given by caps to get the factors of production available to produce product pvec.
    factorsofprod = np.array(cvec)[caps]

    #weights for each capability could be different
    weightsvec = np.ones(len(factorsofprod))

    #I take the corresponding geometric average
    return(weightedproduct(factorsofprod, weightsvec))


def ArithmeticPF(cvec, pvec):
    #positions to index the capabilities of the product defined by pvec
    caps = [i for i,x in enumerate(pvec) if x>0]

    #I subset the capabilities of the city given by caps to get the factors of production available to produce product pvec.
    factorsofprod = np.array(cvec)[caps]

    #weights for each capability could be different
    weights = np.ones(len(factorsofprod))

    return(np.average(factorsofprod, weights=weights))



def ProductionMatrix(Cmat, Pmat, option):
    # In Cmat, rows are countries and capabilities are columns
    # In Pmat, rows are capabilities and columns are products
    
    #option 1: Leontief production function: minimum of the production factors: min(x1,x2,...,xn)
    #option 2: Cobb-Douglas production function: geometric average between the production factors: (x1 x2 ... xn)^(1/n) 
    #option 3: Arithmetic average production function.
    #option 4: ...
    Xcp=None

    if(option==1):
        #Xcp = [[float(Leontief(c, p)) for p in np.transpose(Pmat)] for c in Cmat]
        pcaps = np.sum(Pmat, axis=0)
        pcaps = Replace(pcaps, 0.0, 1.0)
        newPmat = Pmat/pcaps
        Xcp = np.floor(np.dot(Cmat, newPmat))
    elif(option==2):
        Xcp = [[float(CobbDouglas(c, p)) for p in np.transpose(Pmat)] for c in Cmat]
    elif(option==3):
        Xcp = [[float(ArithmeticPF(c, p)) for p in np.transpose(Pmat)] for c in Cmat]
    elif(option==4):
        print("No option=4 yet!")
    else:
        print("No option input!")

    return(np.array(Xcp))



#==============================================================================================
# Measuring complexity 

def ECeigenvecs(McpDF):
    """
    
    """ 
    #this function assumes McpMatrix is already an np.array object
    McpMatrix = np.array(McpDF)

    Nc = McpDF.shape[0]   #number of cities
    Np = McpDF.shape[1] #number of products

    cities = McpDF.index.values
    products = McpDF.columns.values

    #DIVERSIFICATION kc0 AND UBIQUITY kp0
    kc0 = np.sum(McpMatrix, axis=1) + 0.0
    kp0 = np.sum(McpMatrix, axis=0) + 0.0

    kc0 = Replace(kc0, 0.0, 1.0)
    kp0 = Replace(kp0, 0.0, 1.0)
    # the above replacements are so that I don't divide by zero.
    # if any element of kc0 or kp0 are 0, then any "corresponding"
    # element in Mcp is also 0 anyway... hence this division 0/0 should be 0 = 0/1.
    
    Mc2p = (McpDF.T/kc0).T # right-stochastic matrix
    Mp2c = McpDF/kp0 # left-stochastic matrix
    
    # left-stochastic matrices
    Mc2c = Mp2c.dot(Mc2p.T) # country to country
    Mp2p = Mc2p.T.dot(Mp2c) # product to product
    
    Dc, leftVc, rightVc = eig(Mc2c, left=True) #eigenvalues Dc (Dc[i]) and eigenvectors Vc (as columns Vc[:,i])
    Dp, leftVp, rightVp = eig(Mp2p, left=True) #eigenvalues Dp (Dp[i]) and eigenvectors Vp (as columns Vp[:,i])
    
    return((Mc2c, Dc, leftVc, rightVc), (Mp2p, Dp, leftVp, rightVp))


def ECIcalculate(McpDF):
    #this function assumes McpDF is either an np.array object or a pandas data frame
    McpMatrix = np.array(McpDF)

    Nc = len(McpMatrix)   #number of cities
    Np = len(McpMatrix.T) #number of products

    cities = range(Nc)
    products = range(Np)

    #DIVERSIFICATION kc0 AND UBIQUITY kp0
    kc0 = np.sum(McpMatrix, axis=1) + 0.0
    kp0 = np.sum(McpMatrix, axis=0) + 0.0

    Mcc = [ [ 0.0 for cj in cities ] for ci in cities ]

    #Matrix for ECI
    kc0copy = Replace(kc0, 0.0, 1.0)
    kp0copy = Replace(kp0, 0.0, 1.0)
    # the above replacements are so that I don't divide by zero.
    # if any element of kc0 or kp0 are 0, then any "corresponding"
    # element in Mcp is also 0 anyway... hence this division 0/0 should be 0 = 0/1.

    for ci in cities:
        for cj in cities:
            Mcc[ci][cj] = np.dot(McpMatrix[ci], McpMatrix[cj]/kp0copy)/kc0copy[ci]

    Mcc = np.array(Mcc)
    Dc, Vc = linalg.eig(Mcc)  #eigenvalues Dc (Dc[i]) and eigenvectors Vc (as columns Vc[:,i])
    csorted=Dc.real.argsort()

    Kc = Vc.T[csorted[len(csorted)-2]].real #eigenvector corresponding the second largest eigenvalue
    ECI = (Kc - np.average(Kc))/np.std(Kc)  #standardization

    # The values of the eigenvector may be switching positive <=> negative.
    # The location-complexity should be POSITIVELY related to the ubiquity.
    #print(np.corrcoef(kc0,ECI)[0][1])
    if (np.corrcoef(kc0,ECI)[0][1]<0.0):
        ECI = -1.0*np.array(ECI)

    
    return(ECI)

def FITNESScalculate(McpDF, numiter=100):
    #this function assumes McpDF is either an np.array object or a pandas data frame
    McpMatrix = np.array(McpDF)

    Nc = McpMatrix.shape[0]   #number of cities
    Np = McpMatrix.shape[1] #number of products


    # Initialization of values
    Fc = np.ones(Nc) + 0.0
    Qp = np.ones(Np) + 0.0

    # For a given iteration
    for i in np.arange(numiter):
        Fcnew = np.dot(McpMatrix + 0.0, Qp)
        Fcnew = Replace(Fcnew, 0.0, 1.0)
        denom = np.dot(np.transpose(McpMatrix + 0.0), 1.0/Fc)
        if(any(denom==0)):
            denom = Replace(denom, 0.0, 0.1)
        Qpnew = 1.0/denom
        
        #normalize after step
        Fc = Fcnew/np.mean(Fcnew)
        Qp = Qpnew/np.mean(Qpnew)

    return(Fc)

def ECIPLUScalculate(McpDF, numiter=100, meanmethod='geometric'):
    # this function assumes XcpMatrix is already an np.array object
    # IMPORTANTLY, XcpMatrix is the matrix of total export values 
    # per country per product.
    XcpMatrix = np.array(McpDF)

    Nc = len(XcpMatrix)   #number of cities
    Np = len(XcpMatrix.T) #number of products

    # Initialize value
    Xc = np.sum(XcpMatrix + 0.0, axis=1)
    Xp = np.sum(XcpMatrix + 0.0, axis=0)

    # normalize by geometric mean
    Xc = Xc/np.exp(np.mean(np.log(Xc)))

    # For a given iteration
    for i in np.arange(numiter):
        Xc = Replace(Xc, 0.0, 1.0)
        denom = np.dot(XcpMatrix.T + 0.0, 1.0/Xc)
        if(any(denom==0)):
            denom = Replace(denom, 0.0, 0.1)
        difficultyp = 1.0/denom
        Xcnew = np.dot(XcpMatrix + 0.0, difficultyp)
        
        if(meanmethod=='geometric'):
            #normalize after step by the geometric mean
            Xcnew = Replace(Xcnew, 0.0, 0.1)
            Xc = np.exp(np.log(Xcnew)- np.mean(np.log(Xcnew)))
        elif(meanmethod=='arithmetic'):
            #normalize after step by the arithmetic mean
            Xc = Xcnew/np.mean(Xcnew)
        else:
            Xc = Xcnew

    #pdb.set_trace()
    Xp = Replace(Xp, 0.0, 1.0)
    ECIplus = np.log(Xc) - np.log(np.dot(XcpMatrix + 0.0, 1.0/Xp))

    return(Xc)



def create_complexity_series(dflong):
    # This is a function specifically for the SITC-type data
    Xcp_df = pd.pivot_table(dflong,
                    index = ['country_code','year'],
                    columns = 'product_code',
                    values = 'export_value',
                    fill_value=0.0)
    Xcp_df.name = 'None'
    Xcp_df.filename = 'None'

    allyears = dflong.year.unique()

    #Initializing the data frame
    newDF = pd.DataFrame() #creates a new dataframe that's empty

    for yr in allyears:
        indices = list(map(itemgetter(0), Xcp_df.query('year=={:d}'.format(yr)).index.values))
        cols = list(Xcp_df.query('year=={:d}'.format(yr)).columns.values)
        
        newXcp = pd.DataFrame(data=Xcp_df.query('year=={:d}'.format(yr)).values, index=list(indices), columns=cols) + 0.0
        newXcp.index.name = 'country_code'
        newXcp.columns.name = 'product_code'
        
        mcp = McpCalculate(np.array(RCAcalculate(newXcp+0.0)), threshold=1.0)
        
        eci = ECIcalculate(mcp+0.0)
        fitness = FITNESScalculate(mcp+0.0, numiter=100)
        eciplus = ECIPLUScalculate(newXcp.values+0.0, numiter=100)
        
        complexities = np.array([[cty,yr,eci[i],fitness[i],eciplus[i]] for i,cty in enumerate(indices)])
        
        tempdf = pd.DataFrame(complexities, columns=['country_code', 'year', 'ECI', 'FITNESS', 'ECIplus'])
        
        # Concatenate/Append/rbind
        #pd.concat([df1, df2])
        newDF = newDF.append(tempdf, ignore_index = True)

    return(newDF)
#
# ------------------------------------

# ==============================================================================================
# Constructing matrices

def create_single_community(Nc1, Np1, K, cty_typeletter="A", cap_typeletter="A", prod_typeletter="A"):
    # ---------------------
    # Names of countries
    names1 = np.array(["c" + cty_typeletter + str(i).zfill(np.int(np.ceil(np.log10(Nc1*K+1)))) for i in np.arange(1, Nc1+1)])
    # The distribution of complexities (i.e., number of capabilities per country) will be lognormal
    cty_compl1 = sorted(np.ceil(1.0+lognormal(mean=np.log(200.0), sigma=1.5, size=Nc1)), reverse=True)

    # ---------------------
    # Names of products
    prod_names1 = np.array(["p" + prod_typeletter + str(i).zfill(np.int(np.ceil(np.log10(Np1*K+1)))) for i in np.arange(1, Np1+1)])
    # The distribution of product complexities (i.e., number of capabilities required per product) will be binomial
    average_number_of_capabilities_required = 3.0
    prod_compl1 = sorted(np.array(1.0 + binomial(p=average_number_of_capabilities_required/max(cty_compl1), 
                                                    n=np.int(np.ceil(max(cty_compl1))), 
                                                    size=Np1)), reverse=True)

    # ---------------------
    # Number of capabilities
    numcapabilities1 = np.int(max(cty_compl1)+1)
    # Names of capabilities
    cap_names1 = np.array(["a" + cap_typeletter + str(i).zfill(np.int(np.ceil(np.log10(2*numcapabilities1*K+1)))) for i in np.arange(1, numcapabilities1+1)])

    # Each capability has a different likelihood of being used. Some are rare to find, others are very common.
    # This value per capability will serve as a weight to sample it.
    weight_cap_gen1 = 1.0+lognormal(mean=np.log(10.0), size=numcapabilities1)**2
    weight_cap_gen1 = sorted(weight_cap_gen1/weight_cap_gen1.sum(), reverse=True)

    # ---------------------
    # Creating the Cca

    # Creating the list of capabilities that each country has. Each country will thus be represented
    # as a tuple of capabilities.
    cty_endowments1 = {names1[country]: tuple(choice(cap_names1,
                                                    size = np.int(cty_compl1[country]),
                                                    replace = False,
                                                    p = weight_cap_gen1))
                                                    for country in range(0,Nc1) }
    # Creating the edge list
    ctycap_edgelist1 = np.array([[country, capability] for country in names1 for capability in cty_endowments1[country]])


    # ---------------------
    # Creating the Ppa

    # Creating the list of capabilities that each product requires. Each product will thus be represented
    # as a tuple of capabilities.
    prod_requirements1 = {prod_names1[product]: tuple(choice(cap_names1,
                                                            size = np.int(prod_compl1[product]), 
                                                            replace = False, 
                                                            p = weight_cap_gen1)) 
                                                            for product in range(0,Np1) }
    # Creating the edge list
    prodcap_edgelist1 = np.array([[product, capability] for product in prod_names1 for capability in prod_requirements1[product]])

    # ---------------------
    # OUTPUT
    return(ctycap_edgelist1, prodcap_edgelist1)

def create_toy_Cca_Ppa(Nc, Np, K=1):
    LETTERS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
    Nc_i = np.int(np.floor(Nc/(K+0.0)))
    Np_i = np.int(np.floor(Np/(K+0.0)))

    # ---------------------
    # First community
    ki = 0
    letter_i = LETTERS[ki]
    ctycap_edgelist_all, prodcap_edgelist_all = create_single_community(Nc_i, 
                                                                    2*Np_i, # the first community we will make it very large 
                                                                    K, 
                                                                    cty_typeletter=letter_i,
                                                                    cap_typeletter=letter_i, 
                                                                    prod_typeletter=letter_i)

    # ---------------------
    # Looping over communities
    for ki in range(1, int(K)):
        letter_i = LETTERS[ki]
        ctycap_edgelist_i, prodcap_edgelist_i = create_single_community(Nc_i, 
                                                                        Np_i, 
                                                                        K, 
                                                                        cty_typeletter=letter_i,
                                                                        cap_typeletter=letter_i, 
                                                                        prod_typeletter=letter_i)
        # Appending
        ctycap_edgelist_all = np.concatenate((ctycap_edgelist_all, ctycap_edgelist_i), axis=0)
        prodcap_edgelist_all = np.concatenate((prodcap_edgelist_all, prodcap_edgelist_i), axis=0)
    # end of loop for creating communities

    # ---------------------
    # Creating off-diagonal elements (i.e., adding some white noise everywhere in the Cca)
    allcty_names = np.sort(np.unique(ctycap_edgelist_all[:,0]))
    allprod_names = np.sort(np.unique(prodcap_edgelist_all[:,0]))
    allnewcap_names = np.sort(np.unique(ctycap_edgelist_all[:,1]))

    newNa = len(allnewcap_names)
    fraction2fill = 0.2
    ncross = np.int(fraction2fill*Nc*newNa)
    cty_subsampl = np.array(choice(allcty_names, size=ncross))
    cap_subsampl = np.array(choice(allnewcap_names, size=ncross))
    offdiag_edgelist = np.array([[cty_subsampl[i], cap_subsampl[i]] for i in range(0, ncross)])

    bipartite_long = np.concatenate((ctycap_edgelist_all, offdiag_edgelist), axis=0)

    # ---------------------
    # Some edges might have duplicates
    bipartite_long = np.unique(bipartite_long, axis=0)
    prodcap_edgelist_all = np.unique(prodcap_edgelist_all, axis=0)

    # ---------------------
    # Convert the edgelist arrays into a DataFrames
    Cca_longdf = pd.DataFrame(bipartite_long, columns=['country_code', 'capability_code'])
    Cca_longdf['Edge'] = 1.0

    Ppa_longdf = pd.DataFrame(prodcap_edgelist_all, columns=['product_code', 'capability_code'])
    Ppa_longdf['Edge'] = 1.0

    # ---------------------
    # From long to wide format
    Cca_df_wide = Cca_longdf.pivot(index='country_code', columns='capability_code', values='Edge').fillna(0.0)
    Ppa_df_wide = Ppa_longdf.pivot(index='product_code', columns='capability_code', values='Edge').fillna(0.0)

    # The edge list of products may not have the same capabilities, so we need to re-index
    Ppa_df_wide = Ppa_df_wide.reindex(columns = Cca_df_wide.columns.values, fill_value=0.0)

    # ---------------------
    # OUTPUT
    return(Cca_df_wide, Ppa_df_wide)


def create_singleMcp(Nc1, Np1, fill, cty_typeletter="A", prod_typeletter="A"):
    # ---------------------
    # Names of countries
    names1 = np.array(["c" + cty_typeletter + str(i).zfill(np.int(np.ceil(np.log10(Nc1+1)))) for i in np.arange(1, Nc1+1)])

    # ---------------------
    # Names of products
    prod_names1 = np.array(["p" + prod_typeletter + str(i).zfill(np.int(np.ceil(np.log10(Np1+1)))) for i in np.arange(1, Np1+1)])

    # ---------------------
    # Creating the Mcp

    # Creating the list of products that each country has. Each country will thus be represented
    # as a tuple of products.
    cty_products1 = {names1[country]: tuple(choice(prod_names1,
                                                    size = np.int(Np1*fill),
                                                    replace = False))
                                                    for country in range(0,Nc1)}
    # Creating the edge list
    ctyprod_edgelist1 = np.array([[country, product] for country in names1 for product in cty_products1[country]])


    # ---------------------
    # OUTPUT
    return(ctyprod_edgelist1)

def create_toyMcp(Nc, Np, K=1, withinfill=0.8, betweenfill=0.2):
    LETTERS = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
    Nc_i = np.int(np.floor(Nc/(K+0.0)))
    Np_i = np.int(np.floor(Np/(K+0.0)))

    # ---------------------
    # First community
    ki = 0
    letter_i = LETTERS[ki]
    ctyprod_edgelist_all = create_singleMcp(Nc_i, 
                                            Np_i, 
                                            fill=withinfill, 
                                            cty_typeletter=letter_i,
                                            prod_typeletter=letter_i)

    # ---------------------
    # Looping over communities
    for ki in range(1, int(K)):
        letter_i = LETTERS[ki]
        ctyprod_edgelist_i= create_singleMcp(Nc_i, 
                                            Np_i, 
                                            fill=withinfill, 
                                            cty_typeletter=letter_i,
                                            prod_typeletter=letter_i)
        # Appending
        ctyprod_edgelist_all = np.concatenate((ctyprod_edgelist_all, ctyprod_edgelist_i), axis=0)
    # end of loop for creating communities

    # ---------------------
    # Creating off-diagonal elements (i.e., adding some white noise everywhere in the Mcp)
    allcty_names = np.sort(np.unique(ctyprod_edgelist_all[:,0]))
    allprod_names = np.sort(np.unique(ctyprod_edgelist_all[:,1]))

    ncross = np.int(betweenfill*Nc*Np)
    cty_subsampl = np.array(choice(allcty_names, size=ncross))
    prod_subsampl = np.array(choice(allprod_names, size=ncross))
    offdiag_edgelist = np.array([[cty_subsampl[i], prod_subsampl[i]] for i in range(0, ncross)])

    bipartite_long = np.concatenate((ctyprod_edgelist_all, offdiag_edgelist), axis=0)

    # ---------------------
    # Some edges might have duplicates
    bipartite_long = np.unique(bipartite_long, axis=0)

    # ---------------------
    # Convert the edgelist arrays into a DataFrames
    Mcp_longdf = pd.DataFrame(bipartite_long, columns=['country_code', 'product_code'])
    Mcp_longdf['Edge'] = 1.0

    # ---------------------
    # From long to wide format
    Mcp_df_wide = Mcp_longdf.pivot(index='country_code', columns='product_code', values='Edge').fillna(0.0)


    # ---------------------
    # OUTPUT
    return(Mcp_df_wide)
    
    
# ===============================================================================================
# Alternatives to complexity
def distance_to_center(Mcp_widedf, kcomm=None, Mc2c=None, Dc=None, leftVc=None, rightVc=None):
    # Calculating the c2c and p2p matrices, eigenvalues and left-eigenvectors
    #(Mc2c, Dc, leftVc, rightVc), (Mp2p, Dp, leftVp, rightVp) = ECeigenvecs(Mcp_widedf)
    if(all((Mc2c is None, Dc is None, leftVc is None, rightVc is None))):
        (Mc2c, Dc, leftVc, rightVc) = ECeigenvecs(Mcp_widedf)
    

    # Compute the number of communities, if not passed
    if(kcomm is None):
        kcomm = np.sum(Dc.real>0.22)
        print("The number of clusters identified are: {}".format(kcomm))
    else:
        print("The number of clusters identified are: {}".format(np.sum(Dc.real>0.22)))
        print("The number of clusters chosen: {}".format(kcomm))


    # Take only the firts k eigenvectors. Hence, take the diversity plus (k-1) next eigen vectors
    leftX_mat = leftVc[:,1:kcomm]
    rightX_mat = rightVc #[:,0:kcomm]

    # Compute the center
    leftDistances = np.array(list(map(lambda coordinates_vec: np.sqrt(np.dot(10.0*coordinates_vec, 10.0*coordinates_vec)), leftX_mat)))
    rightDistances = np.array(list(map(lambda coordinates_vec: np.sqrt(np.dot(10.0*coordinates_vec, 10.0*coordinates_vec)), rightX_mat)))

    dist_df = pd.DataFrame(np.transpose([leftDistances,rightDistances]), index=Mc2c.index, columns=["leftDist2Origin","rightDist2Origin"])

    return(dist_df)    
    
def cosine_angle_to_target(Mcp_widedf, targetindex=None, kcomm=None, Mc2c=None, Dc=None, leftVc=None, rightVc=None):
    # Constructing the left-eigenspace, and computing the cosine similarity
    # between every observation and a target observation.
    if(all((Mc2c is None, Dc is None, leftVc is None, rightVc is None))):
        (Mc2c, Dc, leftVc, rightVc) = ECeigenvecs(Mcp_widedf)


    # Compute the number of communities, if not passed
    if(kcomm is None):
        kcomm = np.sum(Dc.real>0.22)
        print("The number of clusters identified are: {}".format(kcomm))
    else:
        print("The number of clusters identified are: {}".format(np.sum(Dc.real>0.22)))
        print("The number of clusters chosen: {}".format(kcomm))

    
    # Take only the firts k eigenvectors. Discard, however, the first. 
    # Hence, take the indices from 1 to k+1.
    leftX_mat = leftVc[:,1:kcomm]
    rightX_mat = rightVc#[:,0:kcomm]

    # I need the distances to the origin. 
    leftDistances = np.array(list(map(lambda coordinates_vec: np.sqrt(np.dot(10.0*coordinates_vec, 10.0*coordinates_vec)), leftX_mat)))
    rightDistances = np.array(list(map(lambda coordinates_vec: np.sqrt(np.dot(10.0*coordinates_vec, 10.0*coordinates_vec)), rightX_mat)))
    leftDistances[leftDistances==0.0] = 1.0
    rightDistances[rightDistances==0.0] = 1.0

    # Matrix of Similarities
    leftSim_df = pd.DataFrame(np.dot(np.dot(np.diag(1.0/leftDistances), np.dot(leftX_mat, leftX_mat.T)), np.diag(1.0/leftDistances)), index=Mc2c.index, columns=Mc2c.index)
    rightSim_df = pd.DataFrame(np.dot(np.dot(np.diag(1.0/rightDistances), np.dot(rightX_mat, rightX_mat.T)), np.diag(1.0/rightDistances)), index=Mc2c.index, columns=Mc2c.index)

    return(pd.DataFrame(np.transpose([leftSim_df[targetindex].values, rightSim_df[targetindex].values]), index=Mc2c.index, columns=["leftAngle2target", "rightAngle2target"]))   

    
def distance_to_target(Mcp_widedf, targetindex=None, kcomm=None, Mc2c=None, Dc=None, leftVc=None, rightVc=None):
    # Constructing the left-eigenspace, and computing the cosine similarity
    # between every observation and a target observation.
    if(all((Mc2c is None, Dc is None, leftVc is None, rightVc is None))):
        (Mc2c, Dc, leftVc, rightVc) = ECeigenvecs(Mcp_widedf)
    

    # Compute the number of communities, if not passed
    if(kcomm is None):
        kcomm = np.sum(Dc.real>0.22)
        print("The number of clusters identified are: {}".format(kcomm))
    else:
        print("The number of clusters identified are: {}".format(np.sum(Dc.real>0.22)))
        print("The number of clusters chosen: {}".format(kcomm))

    # Matrix of coordinates
    leftX_mat = leftVc[:,1:kcomm]
    rightX_mat = rightVc#[:,0:kcomm]

    # I need the distances to the origin. 
    leftPairDistsDf = pd.DataFrame(distance.squareform(distance.pdist(leftX_mat, 'euclidean')), index = Mcp_widedf.index, columns = Mcp_widedf.index)
    rightPairDistsDf = pd.DataFrame(distance.squareform(distance.pdist(rightX_mat, 'euclidean')), index = Mcp_widedf.index, columns = Mcp_widedf.index)
    
    return(pd.DataFrame(np.transpose([leftPairDistsDf[targetindex].values, rightPairDistsDf[targetindex].values]), index=Mcp_widedf.index, columns=["leftDistance2target", "rightDistance2target"]))   

# def generate_Cca(countryendowments):
# countryendowments = cty_endowments1
# countrynames = names1
# capabilities_vec = np.unique([capability for country in countrynames for capability in countryendowments[country]])

# Xi_mat = np.zeros((len(capabilities_vec), len(Pnames_vec)))
# for f in range(Xi_mat.shape[0]):
    # for p in range(Xi_mat.shape[1]):
        # ps = Position(combinations2produceP_mat[p], features_vec[f])
        # if(features_vec[f] in combinations2produceP_mat[p]):
            # Xi_mat[f][p] = coefficients2produceP_vec[p][ps[0]]
# # ADDING A CONSTANT
# # Assumptions: 
# # - the CONSTANT should be negative (so that probabilities are below 0.5 when not producing anything).
# # - the abs(CONSTANT) should be a measure of how difficult the product is to invent.
# # - On the other hand, when producing the product, a country should keep producing it.
# # - Hence: beta0 < 0 and beta0 + beta1 > 0
# # - I will assume beta0 = -beta1 + param
# param = np.log(prob/(1.0-prob))
# Xi_df_noconstant = pd.DataFrame(Xi_mat, index=features_vec, columns=Pnames_vec)
# firstrow = pd.DataFrame([[-coefficients2produceP_vec[p][0]+param for p in Pnames_vec]], index=["constant"], columns=Pnames_vec)
# Xi_df = pd.concat([firstrow, Xi_df_noconstant])

# return(Xi_df)


# ==============================================================================================
def eci_transform(Mcp_widedf, kdim=3):
    # Automatic embedding
    
    # Taking the eigenvectors
    (Mc2c, Dc, leftVc, rightVc) = ECeigenvecs(Mcp_widedf)
    
    
    










