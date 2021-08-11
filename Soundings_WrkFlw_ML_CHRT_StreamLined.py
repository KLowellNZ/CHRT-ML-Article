# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 14:58:11 2018

@author: XXXXX
"""

#%%
############################################################################
# This program is the outcome of a large body of work to find the best way
# to classify .las soundings as Bathy or NotBathy. Ultimately, the only
# sounding data needed is depth -- not any of the sounding meta-data
# explored earlier. Input to this program is a .las file, and output from
# CHRT (aug_hypos and chrt_out but neighbours for each hypothesis are not
# needed). (CHRT is run outside this program to produce the necessary
# data outputs.) Output from this program is 1)a final .csv Bathy/NotBathy
# classification for all soundings, 2)a .txt cross-tabulation of the CHRT-ML
# and NOAA Bathy/NotBathy classification including accuracy/summary
# statistics, a .pdf providing descriptive information on all soundings in
# the .las, and two .csv files giving summaries of all CHRT hypotheses
# formed for all estimation nodes (ENs)), and 5)CHRT hypothesis information
# summarised by EN.
# 
# This program 1)screens out outlier CHRT estimation nodes (ENs), 
# 2)uses k-means clustering on information from CHRT ENs to define
# three clusters representing ocean surface, ocean bottom, and
# somewhere in between, 3)applies rules to the clusters to define a
# "bathymetry depth interval" (BDI) for soundings, 4)reads soundings from
# a .las file and classifies them as Bathy/NotBathy based on the BDI,
# and 5)outputs information aobut the clustering process and the
# accuracy of the classificaiton relative to the NOAA reference
# classification.
#
# NOTE ON "CLIPPED" .LAS INPUT: If any tiles cover area(s) that is/are above
# Mean Sea Level, the .las file that is input to this program should have been
# clipped in "Step 0" to eliminate those areas.  CHRT must then be run on
# the CLIPPED .las file. Clipping can be achieved a number of ways using a
# variety of software, but the fundamental steps are:
# 1)Create an empty shapefile and digitize a geo-referenced polygon that
#   identifies the area(s) above MSL. If there are multiple non-contiguous
#   areas, they should be digitized as a single MultiPolygon (a shapefile
#   data type).
# 2)Clip the .las file using software such as publicly avaialble QGIS,
#   commercially available ArcMap, R libraries that address .las data.
#   (This author found that the R routines did not create correct 
#   shapefile metadata.  This can be resolved by writing a utility program
#   that replaces the metadata of the clipped .las file with the metadata
#   of the original .las file.)
#
############ CALC_WEIGHTED_STATS (from ClusterClassifyPulses) ##############
# This function calculates a weighted average and weighted standard
# deviation of a column from frequencies.
def calc_weighted_stats(df,x,freq):
    wghtavg=(df[x]*df[freq]).sum()/df[freq].sum()
    wghtSD=math.sqrt((((df[x]-wghtavg)**2*df[freq]).sum())/df[freq].sum())
    return wghtavg, wghtSD
####### CONF_LIMITS_3CLSTRS_ALLOVRLAP (from ClusterClassifyPulses) #########
# conf_limits_3clstrs_allovrlap calculates the surface and deep "bathy depth
# interval" limits for winning hypotheses. The df passed to the function is
# a groupby object. This function provides for overlap between mid-depth and
# deep and mid-depth and shallow.  Whether the mid-
# depth "belongs" to the surface or deep cluster is determined by the location
# of its mean. If its mean is within the prop_ratio (distance from Deep and
# Shallow), it belongs to both and the sounding classification CI is based on
# the CIs for shallow, mid, and deep. If it is within the prop_ratio for
# deep, the CI is based on the CI for mid and Deep.
# NOTE: Mean values are positive here -- i.e., a depth of approx +22 is surface.
def conf_limits_3clstrs_allovrlap(dfclass,crosstabfile,conf_level_deep,
                        conf_level_surface,overlap_tol,mtoadd,minENs,
                        minovrlap):
# Create empty dataframe and then fill out with the means, std. devs., and n
# of dpthwin of each cluster.
    dftemp=pd.DataFrame()
    dftemp['mean']=dfclass['dpthwin']['mean'].tolist()
    dftemp['std']=dfclass['dpthwin']['std'].tolist()
    dftemp['count']=dfclass['dpthwin']['count'].tolist()
# Sort it by ascending mean values. 
    dftemp.sort_values('mean',axis=0, inplace=True)
    dftemp.reset_index(drop=True,inplace=True)
# Set up a comparable df for disambiguaiton rules.
    dfwinorno=pd.DataFrame()
    dfwinorno['dpthavg_win']=dfclass['dpthwin']['mean'].tolist()
    dfwinorno['dpthavg_nowin']=dfclass['dpthavg_nowin']['mean'].tolist()
    dfwinorno.sort_values(by='dpthavg_win',inplace=True,ignore_index=True)
# Determine CI_limits for each cluster.             
# A 95% lower/surface limit is employed if the CI limits of the shallow cluster
# are used. Similarly, at times the 95% LL of mid depth is used.
    LL_shall95=dftemp.loc[0,'mean'] - \
             abs(t.ppf((1-0.95)/2,dftemp.loc[0,'count']-1))*dftemp.loc[0,'std']
    LL_mid95=dftemp.loc[1,'mean'] - \
            abs(t.ppf((1-0.95)/2,dftemp.loc[1,'count']-1))*dftemp.loc[1,'std']
# 99.9% CIs for evaluation of overlap.
    LL_shall=dftemp.loc[0,'mean'] - \
             abs(t.ppf((1-conf_level_surface)/2,dftemp.loc[0,'count']-1))*dftemp.loc[0,'std']
    UL_shall=dftemp.loc[0,'mean'] + \
        abs(t.ppf((1-conf_level_deep)/2,dftemp.loc[0,'count']-1))*dftemp.loc[0,'std']
    LL_mid=dftemp.loc[1,'mean'] - \
            abs(t.ppf((1-conf_level_surface)/2,dftemp.loc[1,'count']-1))*dftemp.loc[1,'std']
    UL_mid=dftemp.loc[1,'mean'] + \
        abs(t.ppf((1-conf_level_deep)/2,dftemp.loc[1,'count']-1))*dftemp.loc[1,'std']
    LL_deep=dftemp.loc[2,'mean'] - \
        abs(t.ppf((1-conf_level_surface)/2,dftemp.loc[2,'count']-1))*dftemp.loc[2,'std']
    UL_deep=dftemp.loc[2,'mean'] + \
        abs(t.ppf((1-conf_level_deep)/2,dftemp.loc[2,'count']-1))*dftemp.loc[2,'std']
# To ensure all deeper soundings are identified as Bathy, create a UL_Deep
# that is mtoadd m deeper (mtoadd is recommended to be 1 m).
    UL_bathy=UL_deep+mtoadd
    print('\n*** THE SURFACE CI LIMIT FOR SHALLOW HERE IS 99.9% (NOT 95%) ***',
          file=crosstabfile)
    print('Cluster CIs:  Mean   StdDev   n    ',
              round(conf_level_deep*100,1),'% (LL UL)',
          '\n   Shallow ',str(round(dftemp.loc[0,'mean'],3)).rjust(7),
              str(round(dftemp.loc[0,'std'],3)).rjust(7),
              str(int(dftemp.loc[0,'count'])).rjust(5),
              str(round(LL_shall,3)).rjust(8),
              str(round(UL_shall,3)).rjust(7),
          '\n   MidDpth ',str(round(dftemp.loc[1,'mean'],3)).rjust(7),
              str(round(dftemp.loc[1,'std'],3)).rjust(7),
              str(int(dftemp.loc[1,'count'])).rjust(5),
              str(round(LL_mid,3)).rjust(8),
              str(round(UL_mid,3)).rjust(7),
          '\n   Deepest ',str(round(dftemp.loc[2,'mean'],3)).rjust(7),
              str(round(dftemp.loc[2,'std'],3)).rjust(7),
              str(int(dftemp.loc[2,'count'])).rjust(5),
              str(round(LL_deep,3)).rjust(8),
              str(round(UL_deep,3)).rjust(7),'\n',file=crosstabfile)
# If the number of ENs in a cluster is less than minENs, make n=minENs --
# i.e., where changes in t values from n-1 to n start to taper off --
# so that the CI does not widen excessively.
    if dftemp.loc[0,'count'] < minENs:
        LL_shall95=dftemp.loc[0,'mean'] - \
             abs(t.ppf((1-0.95)/2,minENs-1))*dftemp.loc[0,'std']
        print('   LOW DEGREES OF FREEDOM: LL Shallow now based on n=',minENs,
              file=crosstabfile)
    if dftemp.loc[1,'count'] < minENs:
        LL_mid95=dftemp.loc[1,'mean'] - \
             abs(t.ppf((1-0.95)/2,minENs-1))*dftemp.loc[1,'std']
        print('   LOW DEGREES OF FREEDOM: LL MidDepth now based on n=',
              minENs,file=crosstabfile)
    if dftemp.loc[2,'count'] < minENs:
        LL_deep15=dftemp.loc[2,'mean'] - \
             abs(t.ppf((1-conf_level_deep)/2,minENs-1))*dftemp.loc[2,'std']
        print('   LOW DEGREES OF FREEDOM: LL Deep now based on n=',minENs,
              file=crosstabfile)
# SPECIAL CASE: IF THE AMOUNT OF DATA AND NUMBER OF ESTIMATION NODES ARE
# LIMITED, CONFIDENCE LIMITS WILL NOT CALCULATE PROPERLY. CREATE CIS FROM
# AVAILABLE INFO. THIS IS NOTABLY FOR CHUNK 4 TILE 429000E_27170000N
    if pd.isnull(LL_shall) or pd.isnull(UL_shall) or \
        pd.isnull(LL_mid) or pd.isnull(UL_mid) or \
        pd.isnull(LL_shall) or pd.isnull(UL_shall):
# Redefine CIs using the total number of observations for all three clusters
# as the degrees of freedom.
        fake_n=dftemp.loc[0,'count']+dftemp.loc[1,'count']+dftemp.loc[2,'count']
        LL_shall95=dftemp.loc[0,'mean'] - \
            abs(t.ppf((1-0.95)/2,fake_n-1))*dftemp.loc[0,'std']
        LL_shall=dftemp.loc[0,'mean'] - \
            abs(t.ppf((1-conf_level_surface)/2,fake_n-1))*dftemp.loc[0,'std']
        UL_shall=dftemp.loc[0,'mean'] + \
            abs(t.ppf((1-conf_level_deep)/2,fake_n-1))*dftemp.loc[0,'std']
        LL_mid=dftemp.loc[1,'mean'] - \
            abs(t.ppf((1-conf_level_surface)/2,fake_n-1))*dftemp.loc[1,'std']
        LL_mid95=dftemp.loc[1,'mean'] - \
            abs(t.ppf((1-0.95)/2,fake_n-1))*dftemp.loc[1,'std']
        UL_mid=dftemp.loc[1,'mean'] + \
            abs(t.ppf((1-conf_level_deep)/2,fake_n-1))*dftemp.loc[1,'std']
        LL_deep=dftemp.loc[2,'mean'] - \
            abs(t.ppf((1-conf_level_surface)/2,fake_n-1))*dftemp.loc[2,'std']
        UL_deep=dftemp.loc[2,'mean'] + \
            abs(t.ppf((1-conf_level_deep)/2,fake_n-1))*dftemp.loc[2,'std']
        print('   **** TOO FEW OBSERVATIONS AND ENs & CIs RECALCULATED ****',
            '\n   **** WITH',fake_n,'DEGREES OF FREEDOM ***********************',
              file=crosstabfile)
# Define CI as maximum of upper limits and minimum of lower limits without
# considering nan values.
        UL_bathy=np.nanmax([UL_shall,UL_mid,UL_deep + 1.])
        LL_bathy=np.nanmin([LL_shall,LL_mid,LL_deep + 1.])
# Now start "real" Bathy Depth Interval (BDI) rules -- i.e., based on a
# sufficient amount of data. Determine if Deep & Mid overlap by at
# least minovrlp (usually set to 2 cm)).
    elif LL_deep > (UL_mid-minovrlap):
        print('   No',round(minovrlap*100,0),
              'cm overlap between Deep and Mid',file=crosstabfile)
# NO OVERLAP DEEP WITH MID-DEPTH. Check if LL_deep < LL_mid which indicates
# Deep has high variance. If so, base the sounding CI on the deeper LL_mid.
# LATER: I THINK THE FOLLOWING IS IMPOSSIBLE. IF THERE IS NO DEEP-MID
#        OVERLAP LL_DEEP CANNOT BE < LL MID
#************** CONDITION A: This is the easiest condition. Deep does not
#************** overlap Mid or Shallow. LL_Bathy = 99% LL_Deep.
        print('   CI uses Deep only: 99% LL Deep to UL Deep +',
                  int(mtoadd),'m',file=crosstabfile)
    # If there are < 15 ENs in the cluster, use a CI based on n=15.
        if dftemp.loc[2,'count'] < minENs:
                LL_bathy=LL_deep15
        else:
                LL_bathy=LL_deep
# Determine if Mid-depth and Deep overlap by at minovrlap (usually 2 cm),
# but Mid-Depth and Shallow do not.
    elif LL_deep < (UL_mid-minovrlap) and LL_mid > (UL_shall-minovrlap):
        print('   No',round(minovrlap*100,0),
              'cm overlap between Mid and Shallow',file=crosstabfile)
# Deep and Mid overlap but not Mid and Shallow.
# Check if the difference between mean dpthwin and mean nodpthwin for Mid-depth
# is less than the user-specified threshold for mid_depth. If so, mid_depth
# does not contain Bathy and LL is LL Deep. NOTE: Mid-depth is also not used
# if dpthavg_nowin is deeper than dpthwin (-- i.e., nowin contains Bathy).
#************** CONDITION B: Deep overlaps Mid but Mid has narrow dbn
#************** possibly indicating noise near ocean bottom. 
#************** LL_Bathy = 99% LL_Deep.
        if dfwinorno.loc[1,'dpthavg_win']-dfwinorno.loc[1,'dpthavg_nowin'] \
            <= overlap_tol:
# Mid depth does not contain Bathy
            LL_bathy=LL_deep
            print('   Mid and Deep overlap -- but not Shallow. Mid-depth',
                '\n   has no Bathy. CI uses 99% LL of Deep and UL Deep +',
                        int(mtoadd),'m', file=crosstabfile)
# Mid-depth contains Bathy. However, if the LL_mid is less than LL_shall,
# use LL_deep as the lower limit of the CI.
        else:
# LL_mid is larger than LL_shallow. LL_mid defines LL of CI. Calculate 
# and use the LL of the 95% CI on mid-depth.
#************** CONDITION C: Deep overlaps Mid (but not Shallow) and Mid is
#************** broad enough to have Bathy. Bathy contained in two clusters.
#************** LL_Bathy = 95% LL_Mid.
            LL_bathy=LL_mid95
            print('   Mid and Deep overlap -- but not Shallow. Mid-depth',
                '\n   contains Bathy. CI uses 95% LL Mid-depth to UL Deep +',
                int(mtoadd),'m',file=crosstabfile)
# 99.9% CIs for all three clusters overlap. First check if mid-depth
# contains Bathy. (If not, shallow cannot contain Bathy either.)
    else:
# Check if Mid-depth contains Bathy.
        if dfwinorno.loc[1,'dpthavg_win']-dfwinorno.loc[1,'dpthavg_nowin'] \
            <= overlap_tol:
# Mid-depth does not have Bathy. Check if mid's LL is shallower than Deep's LL.
# If so, Deep has high variability and LL_mid is the Bathy limit.
#************** CONDITION D: All three -- Deep, Mid, & Shallow -- overlap
#************** but Deep is very broad (highly variable/low number of ENs)
#************** and Mid does not contain Bathy (i.e., dbn is narrow).
#************** LL_Bathy = 95% LL_Mid.
                if LL_mid>LL_deep:
                    LL_bathy=LL_mid95
                    print('   All 3 clusters overlap, but Mid-depth has no Bathy.',
                        '\n   But LL Mid > LL Deep indicating high variability',
                        '\n   in Deep CI uses 95% LL of Mid-depth and UL Deep +',
                        int(mtoadd),'m',file=crosstabfile)
# Mid-depth does not have Bathy and Deep is not highly variable. LL_Deep
# is the limit of Bathy.
#************** CONDITION E: All three -- Deep, Mid, & Shallow -- overlap
#************** but Mid does not have Bathy. LL_Bathy = 99% LL_Deep.
                else:
                    LL_bathy=LL_deep
                    print('   All 3 clusters overlap, but Mid-depth has no Bathy.',
                        '\n   CI uses 99% LL of Deep and UL Deep +',
                        int(mtoadd),'m',file=crosstabfile)
# Mid-depth has Bathy. Check if Shallow does. NOTE: Shallow is also not used
# if dpthavg_nowin is deeper than dpthwin (-- i.e., nowin contains Bathy).
        elif dfwinorno.loc[0,'dpthavg_win']-dfwinorno.loc[0,'dpthavg_nowin'] <= \
            overlap_tol:
# Shallow does not have Bathy. However, if LL_mid is shallower than LL_shall,
# base CI on Deep only.
#************** CONDITION F: All three -- Deep, Mid, & Shallow -- overlap
#************** but Shallow does not contain Bathy and Mid is very broad
#************** (highly variable/low number of ENs). LL_Bathy = 99% LL_Deep.
            if LL_shall>LL_mid:
                LL_bathy=LL_deep
                print('   All 3 clusters overlap, Mid-depth has Bathy but not',
                   '\n   Shallow. But LL Shall > LL Mid indicating high variability',
                   '\n   in Mid. CI uses 99% LL of Deep and UL Deep +',
                   int(mtoadd),'m',file=crosstabfile)
#************** CONDITION G: All three -- Deep, Mid, & Shallow -- overlap
#************** but Shallow does not contain Bathy. Mid is not overly broad
#************** (highly variable/low number of ENs). LL_Bathy = 95% LL_Mid.
            else:
                LL_bathy=LL_mid95
                print('   All 3 clusters overlap, but Shallow has no Bathy.',
                        '\n   CI uses 95% LL of Mid-depth and UL Bathy is UL Deep +',
                        int(mtoadd),'m',file=crosstabfile)
# Shallow contains Bathy. Use a 95% CI for the LL/surface limit of Bathy.
#************** CONDITION H: All three -- Deep, Mid, & Shallow -- overlap
#************** and all three contain Bathy. LL_Bathy = 95% LL_Shallow.
        else:
                LL_bathy=LL_shall95
                print('   All 3 clusters overlap, and Shallow has Bathy.',
                        '\n   CI uses 95% LL of Shallow and UL Bathy is UL Deep +',
                        int(mtoadd),'m',file=crosstabfile)
# Print a CR/LF for ease of reading.
    print('',file=crosstabfile)
    return LL_bathy, UL_bathy
#########################  DEFSTATS (from las_read)  ######################
# This function returns the min, max, mean, std dev and CV of a column in
# a dataframe.
def getstats(df,colname):
    minval=float(df[colname].min())
    maxval=float(df[colname].max())
    meanval=float(df[colname].mean())
    stdevval=float(df[colname].std())
    try:
        CVval=float(stdevval/meanval*100)
    except:
        CVval=-9999.9
    return minval, maxval, meanval,stdevval,CVval
################ JOIN_LIST_TO_DF (from ClusterClassifyPulses) ##############
# join_list_to_df adds a list of length x to a df of same length as a 
# column while maintaining the indexing of the df.
def join_list_to_df(df,thelist,colname):
    indicesdf=df.index.tolist()
    dflist=pd.DataFrame(data=thelist,columns=[colname],index=indicesdf)
    df=pd.concat([df,dflist],axis=1)
    return df,dflist
################# MAHAL_DISTS (from ClusterClassifyPulses) #################
# mahal_dists calculates Mahalanobis distances for each estimation node.
# Input is a dataframe with a single column normalised between 0.0 and 1.0
def mahal_dists(df):
# Convert data frame to matrix and get inverse of the covariance matrix.
    arrmahal=df.values
    covmatrix=np.cov(arrmahal,rowvar=False)
    invcovmatrix=np.linalg.inv(covmatrix)
# Get means.
    meanmahal=[]
    for i in range(arrmahal.shape[0]):
        meanmahal.append(list(arrmahal.mean(axis=0)))
# Subtract means and get Mahalanobis distances.
    diff=arrmahal-meanmahal
    mahaldists=[]
    for i in range(len(diff)):
        mahaldists.append(np.sqrt(diff[i].dot(invcovmatrix).dot(diff[i])))
    return mahaldists
########### PRINT_CLUSTER_INFO (from ClusterClassifyPulses) ################
# print_cluster_info prints out clustering information
def print_cluster_info(dfclass,clustercols,LL,UL,crosstabfile):
    varlist=['dpthwin','dpthavg_nowin']
# Print header
    header= '      Cluster'+dfclass['dpthwin'].columns[0].center(7)+' '+ \
        dfclass['dpthwin'].columns[1].center(7)+'  '+ \
        dfclass['dpthwin'].columns[2].center(7)+'  '+ \
        dfclass['dpthwin'].columns[3].center(6)+'  '+ \
        dfclass['dpthwin'].columns[7].center(6)
# Print info for each variable.
    for var in (varlist):
        print('   INFO ON:',var,'(Sorted by ascending dpthwin)',file=crosstabfile)
        print(header,file=crosstabfile)
# Create a non-multi-indexed dataframe for this variable. NOTE: The new
# df will have the previous index -- which is the cluster -- as a column
# labelled "index".
        try:
            dftemp=dfclass[var]
            dftemp.reset_index(inplace=True)
# If processing dpthwin, sort and proceed. If processing dpthavg_nowin,
# create a dataframe using the same order as the sorted dpthwin.
            if var == 'dpthwin':
# Sort on mean. (Order of dpthwin and nodpthwin may be different.)
                dftemp.sort_values('mean', inplace=True, ignore_index=True)
                dforder=pd.DataFrame(columns=['Cluster'])
                dforder['Cluster']=dftemp['Cluster']
# If processing dpthavg_nowin (i.e., NOT dpthwin) join to dataframe.
            else:
                dftemp=pd.merge(left=dforder,right=dftemp,how='left',
                                left_on='Cluster', right_on='Cluster')
        except:
            print('         VARIABLE',var,'WAS REMOVED FROM CLUSTERING DUE TO MISSING DATA.',
                      '\n       AS A CONSEQUENCE, NO STATISTICS ARE AVAILABLE.',
                      file=crosstabfile)
            break
# Print descriptive information about the clusters.
        for i in range(dfclass.shape[0]):
            print('         ',dftemp.loc[i,'Cluster'],
                  str(int(dftemp.loc[i,'count'])).rjust(7),
                  str(round(dftemp.loc[i,'mean'],3)).center(8),
                  str(round(dftemp.loc[i,'std'],3)).center(7), 
                  str(round(dftemp.loc[i,'min'],3)).center(7),
                  str(round(dftemp.loc[i,'max'],3)).center(7),
                  file=crosstabfile)
    print('\nBathy-based interval: Lower limit:',round(LL,3),
          '  Upper limit:',round(UL,3),'\n',file=crosstabfile)
############### PRINT_CROSSTAB (from ClusterClassifyPulses) ###############
# print_crosstab prints a crosstab table properly formatted.
def print_crosstab(dfcross,label1,crosstabfile):
    print(label1,file=crosstabfile)
# Get row and column totals.
    total=sum(dfcross.sum(0))
    sumcols=dfcross.sum(0)
    sumrows=dfcross.sum(1)
# Get row and column names.
    cols_indices=sumcols.index
    rows_indices=sumrows.index
# Print row by row.
    print('\n             NOAA:',end=' ',file=crosstabfile)
    for colname in cols_indices:
        print(colname.rjust(8),end=' ',file=crosstabfile)
    print('   TOTAL',file=crosstabfile)
    print('      CHRT/Kim:        ',file=crosstabfile)
    for rowname in rows_indices:
        print('         ',str(rowname).ljust(8),end=' ',file=crosstabfile)
        for colname in cols_indices:
            print(str(dfcross.loc[rowname,colname]).rjust(8),end=' ',file=crosstabfile)
        print(str(int(sumrows[rowname])).rjust(8),file=crosstabfile)
    print('          TOTAL   ', end=' ',file=crosstabfile)
    for colname in cols_indices:
        print(str(int(sumcols[colname])).rjust(8),end=' ',file=crosstabfile)
    print(str(int(total)).rjust(8),file=crosstabfile)
    return cols_indices
################ PRINT_STATS (from ClusterClassifyPulses) ##################
# print_stats prints summary statistics about a symmetrical confusion matrix.
def print_stats(dfcross,crosstabfile):
    if dfcross.shape[0] == dfcross.shape[1]:
# Global accuracy.
        glblacc=np.trace(dfcross)/dfcross.values.sum()
        print("      Global Accuracy:",round(glblacc,3),
              "\n\n\t           True   False  User's  Producer's",
              file=crosstabfile)
# Calculate various true/false positive/negative rates.
        for i in range(dfcross.shape[0]):
            rowsum=dfcross.sum(axis=1)
            colsum=dfcross.sum(axis=0)
            tp=dfcross.iloc[i,i]/colsum[i]
            fp=(colsum[i]-dfcross.iloc[i,i])/colsum[i]
            users=dfcross.iloc[i,i]/rowsum[i]
            producers=tp
            print('\t',dfcross.columns[i].ljust(8),str(round(tp,3)).rjust(6),str(round(fp,3)).rjust(6),
                  str(round(users,3)).rjust(6),str(round(producers,3)).rjust(9),
                  file=crosstabfile)
    else:
        print('      Confusion matrix not symmetrical.\n      No summary statistics printed.',
              file=crosstabfile)
    print(' ',file=crosstabfile)
################ SORT_CONF_MATRIX (from ClusterClassifyPulses) ############
# sort_conf_matrix arranges colums and rows in the desired order -- i.e.,
# (combined) column representing NotBathy first and (combined) column
# representing Bathy second.
def sort_conf_matrix(df):
    df.sort_index(axis=1, ascending=False, inplace=True)
    df.sort_index(axis=0, ascending=False, inplace=True)
    return df
########################  MAIN PROGRAM  ###################################
########################### LIBRARIES #####################################
import pandas as pd
import numpy as np
import math
from sklearn.cluster import KMeans
from scipy.stats import t
from laspy.file import File
import matplotlib.pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
import time
# Suppress warnings
import warnings
warnings.filterwarnings('ignore')
########################### FILE INFO #####################################
# Naming conventions are important throughout. Moreover,
# the inpath must accommodate clipped (Y:) and unclipped .las files that may
# reside in different places. CHRT outputs chrt_aughypos and chrt_out are
# produced by running CHRT external to this Python program. Output filenames
# of this program are determined after a .las file is read and are based on
# the .las file name.
################## CHUNK 0 -- ORIGINAL FOUR FILES PLUS ONE #################
# Chunk 0 takes 3 hours to process.
chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk0/Chunk0_CHRT/'
laspath='Y:'
laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk0/Chunk0_Clipped/'
outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk0/Chunk0_Step3b_CHRTafter/'
########### FILES: ORIGINAL 4 PLUS CHUNK 1/FILE 3 W/ CHANNELS #############
# filenames=['2016_428000e_2719500n_CLIPPED.las','2016_420500e_2728500n.las',
#             '2016_426000e_2708000n.las','2016_430000e_2707500n.las',
#             '2016_423000e_2727000n.las']
# filenames=['2016_428000e_2719500n_CLIPPED.las']
# filenames=['2016_428000e_2719500n_CLIPPED.las','2016_430000e_2707500n.las']
filenames=['2016_430000e_2707500n.las']
###################### CHUNK 1 FILENAMES ##################################
# Chunk 1 takes 18 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk1/Chunk1_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk1/Chunk1_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk1/Chunk1_Step3B_CHRTafter/'
############### FILES 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ########################
# filenames=['2016_415000e_2730000n.las','2016_415500e_2719000n.las',
#            '2016_423000e_2727000n.las','2016_427500e_2720500n.las',
#             '2016_423000e_2716000n_CLIPPED.las','2016_423500e_2720000n.las',
#             '2016_427000e_2730500n.las','2016_429500e_2715500n.las',
#             '2016_422000e_2728000n.las','2016_417000e_2727500n.las']
###################### CHUNK 2 FILENAMES ##################################
# Chunk 2 takes 15 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk2/Chunk2_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk2/Chunk2_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk2/Chunk2_Step3B_CHRTafter/'
############### FILES 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 ############
# filenames=['2016_417500e_2716000n_CLIPPED.las',
#             '2016_421500e_2725000n_CLIPPED.las',
#             '2016_421500e_2719000n_CLIPPED.las',
#             '2016_418000e_2727000n.las','2016_416000e_2719000n.las',
#             '2016_423000e_2715500n.las','2016_426000e_2719500n.las',
#             '2016_419000e_2723000n.las','2016_424500e_2730500n.las',
#             '2016_425500e_2734500n.las']
######################### CHUNK 3 FILENAMES #################################
# Chunk 3 takes 23 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk3/Chunk3_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk3/Chunk3_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk3/Chunk3_Step3B_CHRTafter/'
############## FILES 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 ###############
# filenames=['2016_424500e_2733000n.las','2016_421000e_2707500n.las',
#             '2016_424000e_2726000n.las','2016_423500e_2721000n.las',
#             '2016_421500e_2728000n.las','2016_426000e_2716000n.las',
#             '2016_421000e_2727000n.las','2016_419000e_2729000n.las',
#             '2016_425000e_2728500n.las','2016_420500e_2717500n_CLIPPED.las']
###################### CHUNK 4 FILENAMES ###################################
# Chunk 4 takes 24 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk4/Chunk4_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk4/Chunk4_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk4/Chunk4_Step3B_CHRTafter/'
############# FILES 31, 32, 33, 34, 35, 36, 37, 38, 39, 40 ################
# filenames=['2016_417000e_2726500n.las','2016_426500e_2728000n.las',
#             '2016_419000e_2732000n.las','2016_429000e_2717000n.las',
#             '2016_424000e_2732000n.las','2016_417500e_2720000n.las',
#             '2016_422500e_2731000n.las','2016_428000e_2722500n.las',
#             '2016_420000e_2726500n.las','2016_418000e_2721500n.las']
###################### CHUNK 5 FILENAMES ###################################
# Chunk 5 takes 9.5 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk5/Chunk5_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk5/Chunk5_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk5/Chunk5_Step3B_CHRTafter/'
############# FILES 41, 42, 43, 44, 45, 46, 47, 48, 49, 50 ################
# filenames=['2016_419500e_2726000n.las','2016_430500e_2709000n.las',
#             '2016_420000e_2723500n.las','2016_422500e_2730000n.las',
#             '2016_424000e_2707500n.las','2016_423500e_2726500n.las',
#             '2016_416000e_2725500n.las','2016_418500e_2706500n.las',
#             '2016_420500e_2707500n.las','2016_420000e_2707500n.las']
###################### CHUNK 6 FILENAMES ###################################
# Chunk 6 takes 19 hours to process.
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk6/Chunk6_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk6/Chunk6_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk6/Chunk6_Step3B_CHRTafter/'
############# FILES 51, 52, 53, 54, 55, 56, 57, 58, 59, 60 #################
# filenames=['2016_416000e_2722500n.las','2016_427500e_2723000n_CLIPPED.las',
#             '2016_427000e_2729500n.las','2016_423000e_2708500n.las',
#             '2016_420500e_2725000n.las','2016_420500e_2721500n.las',
#             '2016_426500e_2723000n_CLIPPED.las',
#             '2016_418000e_2715000n_CLIPPED.las',
#             '2016_425500e_2719500n_CLIPPED.las',
#             '2016_415500e_2726500n.las']
############################ CHUNK 7 ######################################
# Chunk 7 takes 17 hours to process
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk7/Chunk7_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk7/Chunk7_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk7/Chunk7_Step3B_CHRTafter/'
############# FILES 61, 62, 63, 64, 65, 66, 67, 68, 69, 70 #################
# filenames=['2016_418500e_2725500n.las', '2016_423000e_2719500n.las',
#             '2016_427500e_2727000n.las','2016_419000e_2717500n_CLIPPED.las',
#             '2016_422000e_2723000n.las',
#             '2016_426000e_2708000n.las', '2016_415500e_2724500n.las',
#             '2016_429000e_2715500n.las', '2016_425000e_2730500n.las',
#             '2016_415500e_2728000n.las']
###################### CHUNK 8 FILENAMES ###################################
# 11.5 hours to run this chunk (primarily because of the first three files.)
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk8/Chunk8_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk8/Chunk8_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk8/Chunk8_Step3B_CHRTafter/'
############### FILES 71, 72, 73, 74, 75, 76, 77, 78, 79 & 80 ##############
# filenames=['2016_424500e_2722500n.las', '2016_425000e_2725500n.las',
#             '2016_425000e_2714500n.las', '2016_422500e_2729000n.las',
#             '2016_428000e_2724000n.las', '2016_416000e_2718500n.las',
#             '2016_420000e_2731000n.las', '2016_423000e_2731500n.las',
#             '2016_416000e_2718000n.las','2016_421000e_2725000n_CLIPPED.las']
###################### CHUNK 9 FILENAMES ###################################
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk9/Chunk9_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk9/Chunk9_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk9/Chunk9_Step3B_CHRTafter/'
############ FILES 81, 82, 83, 84, 85, 86, 87, 88, 89, & 90 #################
# filenames=['2016_419000e_2715000n_CLIPPED.las','2016_417000e_2722000n.las',
#             '2016_424500e_2722000n.las','2016_417000e_2716500n.las',
#             '2016_426000e_2726000n_CLIPPED.las','2016_427000e_2724000n.las',
#             '2016_426000e_2710000n.las','2016_419000e_2706500n.las',
#             '2016_427500e_2716500n_CLIPPED.las',
#             '2016_419500e_2715500n_CLIPPED.las']
########################## CHUNK 10 FILENAMES ##############################
# chrtpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk10/Chunk10_CHRT/'
# laspath='Y:'
# laspathclip='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk10/Chunk10_Clipped/'
# outpath='C:/LAS_Kim/ML_CHRT/ML_CHRT_Data/Chunk10/Chunk10_Step3B_CHRTafter/'
############################# FILES: 91 - 100 ##############################
# filenames=['2016_423500e_2727000n_CLIPPED.las','2016_426500e_2734500n.las',
#             '2016_418000e_2707000n.las','2016_427000e_2709000n.las',
#             '2016_418000e_2717500n_CLIPPED.las',
#             '2016_421000e_2715500n_CLIPPED.las',
#             '2016_427500e_2718000n_CLIPPED.las',
#             '2016_428500e_2710000n.las','2016_423500e_2731500n.las',
#             '2016_424500e_2726500n.las']
############################################################################
################## CLUSTERCLASSIFYPULSES HYPERPARAMETERS ###################
# Column names for various files.
outcols=['row','col','x','y','numbhyps','totpnts','numwinpnts','numnowinpnts',
                 'pntsperhypo_all','pntsperhypo_nowin','pntsperhypo_win',
                'dpthwin','dpthwin_SD','dpthavg_all','dpthSD_all','dpthavg_nowin',
                'dpthSD_nowin','winstrength']
# NOAAclassdict remaps the NOAA classes input to this program to Bathy
# or NotBathy.
NOAAclassdict={'Un':'NotBathy','Gr':'NotBathy','LP-Nz':'NotBathy',
               'Bth':'Bathy','TBWS':'NotBathy'}
# clusterdict maps descriptive information to the column "clstr" in the
# df that is output to the node summary file.
clusterdict={'1.0':'1','0.0':0,'nan':'Outlier'}
# clustercols are the columns on which clustering will be based.
clustercols=['dpthwin','dpthavg_nowin']
# maxpendepth and mindepth are gross error traps.  maxpendepth is the
# maximum depth (below the ellipsoid) to which it can reasonably be
# expected that lidar will penetrate. mindepth is the minimum depth at
# which a sounding might still be possible given tides, etc.
# NOTE: A value of +/- 20 is the ocean surface.  Therefore a maxpendepth
# value of 40 (for example) indicates an ocean depth of 20 m; a mindepth
# of 18 suggests an ocean "depth" 2m above the surface.
# CHANGE: mindepth is now determined by the VDATUMvalue for each tile.
maxpendepth=40
mindepth=20
# Mahal_vars identifies the variables that will be used in screening out
# estimation nodes using Mahalanobis distances.
Mahal_vars=['dpthwin', 'dpthwin_SD', 'dpthavg_all', 'dpthSD_all',
       'dpthavg_nowin', 'dpthSD_nowin', 'winstrength']
# sparse_Mahal_vars are used if there are very few soundings and therefore
# very few estimation nodes which then makes the Mahalanobis matrix non-
# singular. Thing variables to the bare minimum based on depth.
sparse_Mahal_vars=['dpthwin','dpthavg_all','dpthavg_nowin']
# minwinsoundings is the minimum number of soundings that a most likely
# depth hypothesis must have for an EN to be used in clustering. This
# was added because despite running CHRT in LIDAR mode which should
# select the hypothesis with the largest number of soundings as the
# MLD/winner, this was not happening for certain ENs in 1)shallower
# areas where 2)there were relatively few ENs.
minwinsoundings=2
# minovrlap determines the minimum distance in m that 99.9% confidence
# intervals between clusters must overlap to confirm that two CIs do
# overlap.
minovrlap=0.02
# mtoadd is the value in metres that can be added to the depth of the
# sounding confidence interval to "ensure" that deeper bathymetry are
# found.
mtoadd=1.
# minENs defines the minimum number of estimation nodes that will be used
# to calculate confidence intervals. If the actual number of ENs in a
# cluster are below minENs, the degrees of freedom used to calculate
# the CI is minENs-1.
minENs=15
# outlier_conf is the percent of the Mahalanobis Distance distribution that
# will be retained for clustering.
outlier_conf=0.999
# conf_levels are the probability widths of the confidence intervals for
# determining which estimation nodes are Bathy. conf_level_surface extends
# the confidence ellipse towards the surface and is big to ensure a bigger
# data capture. conf_level_deep defines how deep to take points, and since
# we use a cluster indicating the deepest depth, we do not want to take
# as many "deep" points.
# NOTE: IF THE CIs ARE TOO NARROW, 3-CLUSTER CLUSTERING MAY NOT FIND OVERLAP
# BETWEEN TWO CLUSTERS THAT REPRESENT BATHYMETRY. LEAVE THEM AT 0.999.
conf_level_surface=0.999
conf_level_deep = 0.999
# overlap_tol is the difference between the MLD depth and the average depth
# of the non-MLD hypotheses for a cluster to determine if the cluster
# contains Bathy.
overlap_tol=0.41
# subset_label describes what "option" is being explored so that one can keep
# track of all the options. 2- vs- 3-clusters, agglomerated vs. non-agglomerated
# CIs, only labelling MLD hypotheses as Bathy or all hypotheses that are in
# the CI range (regardless of whether they are MLD hypos or not). The print
# statement for this is about 40 lines below.
subset_label='\n\t***** 3-clusters/non-agglomerated/sounding class/'+\
    '\n\t      Clustering on windpth and nowindpth *****'
######################### LAS_READ HYPERPARAMETERS #########################
# Columns to keep -- las files read.
cols=['PulseNo','x','y','z','NOAAclass']
# Conversion of NOAA class number in .las file to text.
classnamedict={'1':'Un','2':'Gr','7':'LP-Nz','25':'WC','26':'Bth',
               '27':'TBWS','29':'SubObj','30':'IHOS57Obj'}
# statvars is for the summary statistics for the .las tile being analysed.
# The 1st element in each is the print label, the
# second element is the variable name, and the third is the number of
# bins to use in the histogram. (NOAA class is not in the list because
# it is a special case.)
statvars=[['X coord ','x', 200],['Y coord ','y',200],['Depth(m)','z',200]]
############################ START PROGRAM ################################
############# START BY PROCESSING CHRT ESTIMATION NODE INFORMATION ########
################## START PROCESSING EACH TILE ############################
# prog_tic is to record time for processing all tiles in this set.
prog_tic=time.time()
for ii,indataname in enumerate(filenames):
    tile_tic=time.time()
# tempname strips the "2016_" year from indata_name.
    tempname=indataname.replace('2016_','')
    print('\n\n**** TILE: ',indataname,'(',ii+1,'of',len(filenames),') ****')
# Set up input and output file names.
    inhypo_name=tempname.replace('.las','_chrt_aughypos.txt')
    inchrtout_name=tempname.replace('.las','_chrt_out.txt')
    individ_hypo_name=tempname.replace('.las','_individ_hypos_postCHRT.csv')
    hyposumm_by_EN_name=tempname.replace('.las','_nodehypo_summary_postCHRT.csv')
    crosstab_name=tempname.replace('.las','_crosstab_postCHRT.txt')
    dataout_name=tempname.replace('.las','_Classed_postCHRT.csv')
    pdffilename=tempname.replace('.las','_postChrt.pdf')
# Open file where cross-tabulation of NOAA and CHRT classifications
# will be output.
    crosstabfile= open(outpath+crosstab_name,'w')
    print('***** TILE/.LAS FILE:',indataname,'*****',file=crosstabfile)
# cross_label is for the cross tabulation output. 
    cross_label = '   ***** ' + inhypo_name  + ' *****\n'+ \
                  '\tNOAA vs CHRT Classification Cross-tab'+subset_label
    print(cross_label,file=crosstabfile)
    print('\n   ****************** TUNING PARAMETERS ******************',
          '\n   Min/Max sounding depth (in m below ellipsoid):',
                round(mindepth,1),'/',round(maxpendepth,1),
          '\n   Required overlap betw. two clustered distributions (cm):',
                round(minovrlap*100,1),
          '\n   Added to deep end of sounding confidence interval (m):',
                round(mtoadd,1),
          '\n   Min. deg of freedom used to calculate sounding CI:',minENs,
          '\n   p value for Mahalanobis CI to eliminate outliers:',
                round(outlier_conf,3),
          '\n   p values for confidence limits for clusters:',
          '\n      Surface:',round(conf_level_surface,3),
          '\n      Deep:   ',round(conf_level_deep,3),
          '\n   Min. separation betw. win/nowin depth for Bathy (m)',
                round(overlap_tol,3),
          '\n   Variables used for Mahalanobis distance calculation:\n\t',
          Mahal_vars,file=crosstabfile)
############################# read aughypos file #########################
# Read chrt_aughypos file -- information about all hypos for each estimation
# node -- as a series of lines and strip off the CR/LF characters.
    with open(chrtpath+inhypo_name,'r') as hypo:
        hypolines=[lines.strip() for lines in hypo]
############################## read chrt_out file #########################
# Read chrt_out file -- identifies bounding box, (rows,cols) of estimation
# nodes, and most likely depth hypothesis (plus some ephemera) -- as a
# series of lines and strip off the CR/LF
    with open(chrtpath+inchrtout_name,'r') as chrtout:
        chrtlines=[lines.strip() for lines in chrtout]
# Prepare dfout -- the output dataframe that summarises individual hypos
# for each EN -- and dfallhyposout that will have information on each of
# the hypos for each of the ENs.
    dfout=pd.DataFrame(columns=outcols)
    dfallhyposout=pd.DataFrame(columns=['row','col','x','y','hypoID','depth',
                                        'dpthSD','n_hypo'])
# From header line of the chrt_out file, get number of rows, number of
# columns, and calculate estimation node spacing and (UTM) coordinates
# of lower-left of grid -- i.e., row,col of [0,0].
    chrt1list=chrtlines[0].split()
# NOTE: The header is in (row,col) order for the chrtout file but (col,row)
# order for the hypo file.
    nrows=int(chrt1list[0])
    ncols=int(chrt1list[1])
    xmin=float(chrt1list[2])
    ymin=float(chrt1list[3])
    xmax=float(chrt1list[4])
    ymax=float(chrt1list[5])
    spacing=(ymax-ymin)/nrows
    xorig=xmin+spacing/2
    yorig=ymin+spacing/2
# Set up parameters to use appropriate lines for each estimation node.
# hypo_ind is for chrt_out and the first line is skipped. chrt_ind is
# the index for chrt_out and the first two lines are skipped. tot_records
# is for estimating total processing time.
    hypo_ind=1
    chrt_ind=2
    tot_records=2
    oldrow=0
    print('   Summarising CHRT estimation nodes....\n      Working on row',
          oldrow,'of',nrows)
# Process and collect information on each estimation node.
    time_tic=time.time()
    while tot_records < len(chrtlines):
# Identify most likely depth (MLD) hypothesis from the chrt_output file.
# (MLD was previously referred to as the "winning" hypothesis.)
        winline=chrtlines[chrt_ind].split()
        winrow=int(winline[0])
        wincol=int(winline[1])
        newrow=winrow
# Print processing progress and prepare for new row if this is the last
# column in this row.
        if newrow != oldrow and newrow%25 ==0:
            print('      Working on row',newrow,'of',nrows)
            oldrow=newrow
            time_toc=time.time()
            elapsedtime=(time_toc-time_tic)
            estimatedtime=len(chrtlines)/tot_records*elapsedtime
            print('         Elapsed/estimated CHRT summary time (mins):',
                  round(elapsedtime/60,1),'/', round(estimatedtime/60,1))
        windepth=float(winline[2])
        winvar=float(winline[3])
        numbhypos=int(winline[4])
        winstrength=float(winline[5])
        winobs=int(winline[6])
# All needed info extracted from the chrt_out file. Get the number of
# hypotheses associated with this EN.  Update index to skip/account for
# the variable density nodes that are established around each EN and 
# account for the two lines that are headers for each EN.
        numbsubhypos=int(chrtlines[chrt_ind+1].split()[5])
        chrt_ind = chrt_ind+int(numbsubhypos)+2
# Determine coordinates of the EN.
        xcentre=xorig+wincol*spacing
        ycentre=yorig+winrow*spacing
# Read the hypotheses for this EN into a dataframe from aug_hypos. First
# create empty dataframe that will be used to calculate summary statistics
# for the hypos for the current EN.
        dfhypo=pd.DataFrame(columns=['depth','dpthSD','n_hypo','hypoID'])
        for i in range(numbhypos):
            hypoline=hypolines[hypo_ind+i].split()
            dftemp=pd.DataFrame({'depth':[float(hypoline[0])],
                'dpthSD':[float(hypoline[1])],'n_hypo':[float(hypoline[2])],
                'hypoID':[int(round(float(hypoline[4]),0))]})
            dfhypo=dfhypo.append(dftemp)
# Hypotheses for each EN will be summarised after removing outliers. To
# retain all hypotheses for output, a separate df -- dfhold -- must be
# created and have additional identifier info added to it.
        allhypos=list(dfhypo['hypoID'])
        dfhold=dfhypo.copy()
        dfhold['row']=winrow
        dfhold['col']=wincol
        dfhold['x']=xcentre
        dfhold['y']=ycentre
# Get rid of any hypotheses for this EN that are obvious outliers -- i.e.,
# beyond the depth to which lidar can reasonably be expected to penetrate
# and shallower than mindepth (usually 2m) above the surface (as a buffer
# for tide).This only eliminates the outlier hypothesis from the calculation
# of statistics. Outlier ENs must be eliminated subsequently if they are
# MLD hypotheses or the outlier hypotheses will corrupt the clustering.
        dfhypo=dfhypo[(dfhypo['depth'] < maxpendepth) &
                      (dfhypo['depth'] > mindepth)]
        dfhypo.reset_index(drop=True,inplace=True)
        if dfhypo.shape[0] > 0:
# Calculate various summary stats for this EN.
            tothypopnts=int(dfhypo['n_hypo'].sum())
            numnowinpnts=float(tothypopnts-winobs)
            pntsperhypo_all=float(tothypopnts)/float(numbhypos)
# It is possible to have only 1 hypothesis per supercell -- i.e., 0 losers.
# If this occurs, division by zero results. Avoid this.
            try:
                pntsperhypo_nowin=float(numnowinpnts)/(numbhypos-1)
            except:
                pntsperhypo_nowin=0
            pntsperhypo_win=float(winobs)/1
# Calculate weighted average and standard deviation for all pulse returns.
            wghtavg,wghtSD= calc_weighted_stats(dfhypo,'depth','n_hypo')
# Create dataframe without winning hypothesis and get weighted stats of
# non-winning hypotheses. (Get hypo with depth closest to winning depth.)
            winindex=dfhypo['depth'].sub(windepth).abs().idxmin()
            dfhypo=dfhypo.drop(winindex)
            wghtavg_nowin, wghtSD_nowin = calc_weighted_stats(dfhypo,'depth','n_hypo')
# Set up output dataframe. Then hold df and append to the output dataframe.
            nowinhypos=list(dfhypo['hypoID'])
            winninghypo=list(set(allhypos).difference(nowinhypos))[0]
            dfhold['winhypo']=np.where(dfhold['hypoID']==winninghypo,1,0)
        else:
            dfhold['winhypo']=np.where(dfhold['depth']==windepth,1,0)
# REMINDER: dfallhyposout contains the information on individual hypos
# for each node.
        dfallhyposout=dfallhyposout.append(dfhold,sort=False)
# If there is one hypothesis only for the current EN, make all stats =
# winstats and nowin stats are set to missing values.
        if numbhypos == 1:
            wghtavg=windepth
            wghtSD=winvar**0.5
            wghtavg_nowin=None
            wghtSD_nowin=None
# Create a one-line df for appending to larger df. Then append.
# REMINDER: dfout summarises information for hypos belonging to each EN.
        dfouttemp=pd.DataFrame({'row':[winrow],'col':[wincol],'x':[xcentre],
                    'y':[ycentre],'numbhyps':[float(numbhypos)],
                    'totpnts':[float(tothypopnts)],'numwinpnts':[float(winobs)],
                    'numnowinpnts':numnowinpnts,'pntsperhypo_all':pntsperhypo_all,
             'pntsperhypo_win':pntsperhypo_win,'pntsperhypo_nowin':pntsperhypo_nowin,
                    'dpthwin':[windepth],'dpthwin_SD':winvar**0.5,
                    'dpthavg_all':[wghtavg],'dpthSD_all':[wghtSD],
                    'dpthavg_nowin':[wghtavg_nowin],'dpthSD_nowin':[wghtSD_nowin],
                    'winstrength':[winstrength]})
        dfout=dfout.append(dfouttemp,sort=False)
# Update hypo_ind and tot_records
        hypo_ind += numbhypos+1
        tot_records += numbsubhypos+2
# Reset the index, and output the dataframe that summarises individual
# hypos for each EN that does not include estimation nodes with no data
# -- i.e., no soundings.
    dfout.reset_index(drop=True, inplace=True)
    dfout=dfout[outcols]
# Information needed for clustering now available. Begin clustering
# process.
    cluster_tic=time.time()
    print('\n   Total number of estimation nodes:',dfout.shape[0],file=crosstabfile)
# Eliminate variables with missing data from Mahalanobis screening. Update
# list of variables to be used in Mahalanobis screening. (This is necessary
# when a .las file has very few points -- like Chunk 4 File 34 -- or very
# few points for an estimation node which can occur for tiles that are not
# completely filled and/or have "rogue" individual points.
    dfnomiss=dfout.dropna(axis='rows')
    print('      Estimation nodes having missing data:',
          dfout.shape[0]-dfnomiss.shape[0],
          '\n        (Includes missing data due to dropping too\n         deep/shallow hypotheses)',
          file=crosstabfile)
# Calculate Mahalanobis distances for estimation nodes.
    print('   ENs with missing data dropped\n',
          '  Eliminating aberrant and outlier ENs')
# Subset variables to prep a df for clustering that will define the depth
# limits for which most likely depths (MLDs) are considered Bathy. First
# calculate normalised Mahalnobis distance to eliminate aberrant
# estimation nodes. (Do not use the variable pntsperhypo_win because it
# is perfectly correlated with numwinpnts.)
    dfmahal=dfnomiss[Mahal_vars].copy()
    nomisstot=dfnomiss.shape[0]
# A copy of nomisstot will be needed to screen for CHRT ENs where the
# MLD only has one associated sounding. Modify it in the same way as dfmahal.
    dfnomiss=dfnomiss[dfnomiss['dpthwin']<maxpendepth]    
    dfmahal=dfmahal[(dfmahal['dpthwin']<maxpendepth)]
    toodeeprows=nomisstot-dfnomiss.shape[0]
    preshallow=dfnomiss.shape[0]
    dfnomiss=dfnomiss[dfnomiss['dpthwin']>mindepth]
    dfmahal=dfmahal[dfmahal['dpthwin']>mindepth]
    tooshallowrows=preshallow-dfnomiss.shape[0]
# Get numwinpnts as a list to join it after Mahalanobis distances
# are calculated.
    numwinlist=dfnomiss['numwinpnts'].to_list()
# In some cases where soundings have weird spatial arrangements,
# dpthavg_nowin has missing data which means it is read as a string which
# in turn eventually converts everything to a string. Make sure all 
# variables are float. Then normalise.
    dfmahal=dfmahal.astype('float64')
    dfmahalnorm=(dfmahal-dfmahal.min())/(dfmahal.max()-dfmahal.min())
# If there are very few soundings (and therefore estimation nodes), the
# Mahalanobis matrix will be non-singular. Reduce the number of variables
# considered.
    if (len(Mahal_vars)>=dfmahalnorm.shape[0]):
        dfmahalnorm=dfmahalnorm[sparse_Mahal_vars]
# mahal_dists returns a list of Mahalanobis distances for the df passed
# to it. Join the list of Mahal dists to the not normalised Mahal df.
    mahaldists=mahal_dists(dfmahalnorm)
    dfmahal, dfmahallist=join_list_to_df(dfmahal,mahaldists,'MahalDist')
    dfmahal,dfnumwinlist=join_list_to_df(dfmahal,numwinlist,'numwinpnts')
# Drop outlier ENs -- anything greater than the outlier_conf p value
# that is set initially at 0.999 -- i.e., 99.9% of all ENs will be
# retained.
    UL_mahal=dfmahal['MahalDist'].mean()+ \
        abs(t.ppf((1-outlier_conf)/2,len(mahaldists)))*dfmahal['MahalDist'].std()
    LL_mahal=dfmahal['MahalDist'].mean()- \
        abs(t.ppf((1-outlier_conf)/2,len(mahaldists)))*dfmahal['MahalDist'].std()
    preoutlrows=dfmahal.shape[0]
# Get df containing observations w/o missing data and that are not outliers.
    dfmahal=dfmahal[(dfmahal['MahalDist']<=UL_mahal) &
                    (dfmahal['MahalDist']>=LL_mahal)]
    outlrrows=preoutlrows-dfmahal.shape[0]
# Now eliminate rows where numwinpnts is less than minimum specified.
    prenwineq1=dfmahal.shape[0]
    dfmahal=dfmahal[dfmahal['numwinpnts']>=minwinsoundings]
    minwindrops=prenwineq1-dfmahal.shape[0]
# Drop obvious dpthwin outliers from the data that will be used for
# clustering. If this is not done, the standard deviation of the cluster
# that will contain Bathy ENs inflates.
    print('      "Too deep" most likely EN depths dropped:',toodeeprows, \
          '\n      "Too shallow" most likely EN depths dropped:',tooshallowrows, \
          '\n      ENs whose MLD has <',minwinsoundings,'soundings dropped:',
                   minwindrops,
          '\n      Outlier ENs (Mahal Dist) dropped:',outlrrows, \
          '\n   Total ENs for clustering:',dfmahal.shape[0],
          '\n   Clustering ENs on Most Likely Depth hypotheses')
    print('      "Too deep" most likely EN depths dropped:',toodeeprows, \
          '\n      "Too shallow" most likely EN depths dropped:',tooshallowrows, \
          '\n      ENs whose MLD has <',minwinsoundings,'soundings dropped:',
                   minwindrops,
          '\n      Outlier ENs (Mahal Dist) dropped:',outlrrows, \
          '\n   Total ENs for clustering:',dfmahal.shape[0],
          '\n\n   Clustering on:',clustercols, file=crosstabfile)
# Renormalise the sanitised -- non-outlier non-deptherrs -- data.
    dfmahalnorm2=(dfmahal-dfmahal.min())/(dfmahal.max()-dfmahal.min())
# Cluster -- variables on which to cluster are specified
# in Hyperparamaters -- to define limits of Bathymetry.
    clusterdata2=dfmahalnorm2[clustercols].values
    kmeans2=KMeans(n_clusters=3, init='random').fit(clusterdata2)
    dfmahal, dfmahallist=join_list_to_df(dfmahal,kmeans2.labels_,'Cluster')
# Summarise by cluster.
    dfclass=dfmahal.groupby('Cluster').describe()
# conf_limits_3clstrs_allovrlap examines overlap between the deep and
# mid-depth clusters and the mid-depth to shallow clusters to calculate
#  confidence intervals.
    LL_bathy, UL_bathy = conf_limits_3clstrs_allovrlap(dfclass,crosstabfile,
                              conf_level_deep,conf_level_surface,
                              overlap_tol,mtoadd,minENs,minovrlap)
# Write clustering info for dpthwin and dpthavg_nowin to .txt crosstab
# file.
    print_cluster_info(dfclass,clustercols,LL_bathy,UL_bathy,crosstabfile)
# Assign each estimation node -- including missing data and outliers -- to
# Bathy or NotBathy. Being classed as a Bathy node means the MLD hypo
# is Bathy -- but non-MLD hypos are always assigned to NotBathy.
    dfout['BathyOrNo']=np.where((dfout['dpthwin']>=LL_bathy) &  \
                            (dfout['dpthwin']<=UL_bathy),'Bathy','NotBathy')
# Assign each hypothesis -- including missing and outliers -- as Bathy or
# NotBathy based on the confidence limits.
####### THE FOLLOWING ONLY DENOTES A HYPOTHESIS AS BATHY IF THE EN IS A
####### BATHY EN AND IF IT IS THE MLD HYPOTHESIS. AN ALTERNATIVE CONSIDERED
# WAS ASSIGNING A HYPO TO BATHY IF ITS EN IS BATHY AND ITS DEPTH IS IN THE
# CI RANGE.
    dfallhyposout['BathyOrNo']=np.where((dfallhyposout['depth']>=LL_bathy) &  \
                            (dfallhyposout['depth']<=UL_bathy) & \
                            (dfallhyposout['winhypo']==1),'Bathy','NotBathy')
# Assign cluster number to each EN using an index to index match. Convert
# the column to string, and change missing values to "Outlier".
    dfout['clstr']=dfmahal['Cluster']
    dfout['clstr']=dfout['clstr'].astype(str)
    dfout['clstr']=dfout['clstr'].map(clusterdict)
    dfout.to_csv(outpath+hyposumm_by_EN_name, index=False)
    dfallhyposout.reset_index(drop=True, inplace=True)
    dfallhyposout.to_csv(outpath+individ_hypo_name, index=False)
    print('      Time to complete clustering (min):',
          round((time.time()-cluster_tic)/60,1))
# Free up memory
    del dfallhyposout
    del dftemp
    del dfhold
    del dfouttemp
    del dfnomiss
    del dfmahal
    del dfmahalnorm
    del dfmahalnorm2
    del dfmahallist
############################################################################
################# SOUNDINGS-BASED CLASSIFICATION STARTS HERE ###############
############################################################################
# At this point, the limits of the Bathymetry Depth Interval have been
# determined.  Use these to classify individual pulses as Bathy/NotBathy,
# summarise, and output the df. 
# First read the .las file and map NOAA's classes to Bathy/NotBathy.
    print('\n   Reading .las file....')
    if 'CLIPPED' in indataname:
        inlaspath=laspathclip
    else:
        inlaspath=laspath
    las_tic=time.time()
    lasfile = File(inlaspath+indataname, mode = "r")
    print('   File read. Getting variables....')
# Scale and offset factors are grabbed for output in the metadata file but
# are not used because the assignment of lasfile.n does the scaling and
# offset.
    offsets=lasfile.header.offset
    scalefactors=lasfile.header.scale
# Start filling a dataframe.
    dfpulse=pd.DataFrame(columns=cols)
    dfpulse['x']=lasfile.x
    dfpulse['PulseNo']=dfpulse.index
    dfpulse['y']=lasfile.y
    dfpulse['z']=lasfile.z
    dfpulse['NOAAclass']=lasfile.Classification.astype(str)
    dfpulse['NOAAclass']=dfpulse['NOAAclass'].map(classnamedict)
    tot_rtrns=dfpulse.shape[0] 
# Get observations by NOAA class for subsequent output.
    dfNOAAclass=dfpulse.groupby(['NOAAclass']).size().reset_index(name='count')
    dfNOAAclass['pct']=dfNOAAclass['count']/tot_rtrns*100
# NOw map to Bathy/NotBathy
    dfpulse['NOAAclass']=dfpulse['NOAAclass'].map(NOAAclassdict)
# Classify soundings.
    print('   Assigning Bathy/NotBathy to individual soundings....')
    dfpulse['CHRTclass']=np.where((dfpulse['z']>=-UL_bathy) &  \
                            (dfpulse['z']<-LL_bathy),'Bathy','NotBathy')
# Print summary stats to be able to evaluate initial classification.
    print('   Calculating summary statistics....')
    classy=dfpulse['NOAAclass'].values
    chrty=dfpulse['CHRTclass'].values
    dfcross=pd.crosstab(chrty,classy,rownames=[ 'chrty'],colnames=['classy'])
    dfcross=sort_conf_matrix(dfcross)
    print('******************************************'.center(46),
          file=crosstabfile)
    label1='   ***** CHRT Class vs. NOAA Class ******'
    cols_indices=print_crosstab(dfcross,label1,crosstabfile)
# Write summary stats to the crosstab .txt file.
    print_stats(dfcross,crosstabfile)
    print('   Outputting pulse file....')
    dfpulse.to_csv(outpath+dataout_name,index=False)
# Variables have been extracted and are in a dataframe. Now print metadata
# including spatial information.
    minx,maxx,meanx,stdevx,CVx= getstats(dfpulse,'x')
    miny,maxy,meany,stdevy,CVy= getstats(dfpulse,'y')
# Projection parameters must be extracted from vlrs[1] and projection
# description must be extracted from vlrs[2].
    prj_params=lasfile.header.vlrs[1].parsed_body
    prj_string=str(lasfile.header.vlrs[2].parsed_body)
    false_easting=prj_params[1]
    false_northing=prj_params[2]
    cntrl_meridian=prj_params[3]
    scale_factor=prj_params[4]
    lat_origin=prj_params[0]
    linear_unit=prj_params[5]
# Now get some specific elements. Split the string and then loop through
# the elements looking for specific things.
    prj_elements=prj_string.split('|')
    for element in prj_elements:
        equ_radius=prj_params[6]
        recip_flattening=prj_params[7]
        rad2deg=prj_params[9]
        if 'projection' in element:
            projname=element.split(':')[1]
        if 'PCS' in element:
            projcoorsys=element.split('=')[1]
        if 'GCS' in element:
            geogcoorsys=element.split('=')[1]
        if 'Datum' in element:
            datum=element.split('=')[1]
        if 'Ellipsoid =' in element:
            ellipsoid=element.split('=')[1]
        if 'Primem' in element:
            prmmeridian=element.split('=')[1]
# Print geographic information to file.
    print('\n***** .LAS FILE DESCRIPTIVES (SEE ALSO .PDF FILE) *****',
          '\n   las file:',indataname,'\n\n   GEOGRAPHIC INFORMATION',
      '\n     Geographic Coordinate System:',geogcoorsys,'\n     Datum:',datum,
      '\n       Prime Meridian:',prmmeridian,
      '\n       Angular unit: Degree (Radians to degree conversion (',round(rad2deg,6),'rad = 1 degree))',
      '\n\n   PROJECTED COORDINATE SYSTEM:',projcoorsys,'\n     Projection:',projname,
      '\n     Ellipsoid:',ellipsoid,'\n       Equatorial radius (m):',equ_radius,
      '\n       Flattening reciprocal:',recip_flattening,
      '\n     False Easting=',false_easting,'\n     False Northing=',false_northing,
      '\n     Central Meridian:',cntrl_meridian,'\n     Scale Factor:',scale_factor,
      '\n     Latitude of Origin:',lat_origin,'     Linear Unit: Meter (',linear_unit,')',
      '\n\n   BOUNDING AREA:\n     Min x, min y:',minx,',',miny,'\n     Max x, max y:',maxx,',',maxy,
      '\n\n   LIDAR INFORMATION:\n     Total points:',tot_rtrns,
      '\n     X, Y, Z offsets:',offsets,'\n     Scale factors:',scalefactors,'\n',
      file=crosstabfile)
# The following returns a df rather than a groupby object. The variable name
# is 'count.'
    print('   *** NOAA: [Bth] = Bathy[Un, Gr, Lp-NZ]=NotBathy ***',
         '\n            [Un, Gr, Lp-NZ]=NotBathy',file=crosstabfile)
    print('\n   NOAA Class Frequency\n     Class    Count     %',file=crosstabfile)
    for row in range(dfNOAAclass.shape[0]):
        print('     {} {:7d}   {:4.1f}'.format(dfNOAAclass.loc[row,'NOAAclass'].ljust(6),
                        int(dfNOAAclass.loc[row,'count']),
                        round(dfNOAAclass.loc[row,'pct'],1)),
                        file=crosstabfile)
    print('     {} {:7d}'.format('Total '.ljust(6),tot_rtrns), file=crosstabfile)
# Get stats for each variable. Certain things will be of minimal interest,
# but even the mean of X,Y cooredinates, for example will indicate if the
# pulse returns are evenly distributed. Also print distributions.
    print('\n   Variable        Min.          Max.          Mean',
      '   Std Dev    CV(%)',file=crosstabfile)
# Create the histogram file.
    outpdf= PdfPages(outpath+pdffilename)
# Set it up to plot the 3 frequency distributions -- x, y, and depth -- on
# the same page.
    plot.rc('font',size=12)
    fig,axs=plot.subplots(nrows=3,ncols=1,figsize=(8,11))
    fig.tight_layout(pad=3.0)
    plotrow=-1
    for var in statvars:
        plotrow += 1
        minval,maxval,meanval,stdvval,CVval= getstats(dfpulse,var[1])
        print('   {} {:13.3f} {:13.3f} {:13.3f} {:8.3f} {:8.3f}'.
              format(var[0],minval, maxval, meanval, stdvval, CVval),
              file=crosstabfile)
# Plot histograms and save to pdf. Width of histogram is the minimum
# and maximum values except for depth which is constrained to the 
# "reasonable limits" of depth of lidar.
        historange=[minval,maxval]
# For the range for the plot of depth, use a clipped 99% confidence interval.
        if var[1] == 'z':
            historange=[-36,-20]
        ax=axs[plotrow]
        ax.hist(dfpulse[var[1]],bins=var[2],range=historange)
        ax.set_title(indataname+' freq.: '+var[0])
    outpdf.savefig()
    plot.show()
    print('   Time to process tile (min):',round((time.time()-tile_tic)/60,1))
    print('\n***** Time to process tile (min):',round((time.time()-tile_tic)/60,1),
          '*****', file=crosstabfile)
# Write output files that have info about the pulses in this file.
    crosstabfile.close()
    outpdf.close()
print('\n***** Total time to process all tiles (min):',
           round((time.time()-prog_tic)/60,1))
