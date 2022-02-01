import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import LocalOutlierFactor
import os
"""
Main Functions

* ExtractFeatures(...)
* KmerFeatures(...)
* BarPLot(...)
* Zoom(...)
* ScatterPlot_2Feat(...)
* ScatterPlot_3Feat(...)
"""
def ExtractFeatures(df) :      
    """
    ExtractFeatures function extracts features from JACUSA2 CALL2 output as a dataframe columns: 
    'Ref','Pos','Base','strand','Cov1','Cov2','Min_cov','Mis','Ins','Del','dA','dC',
    'dG','dT','A1','C1','G1','T1','A2','C2','G2','T2','Ins11','Ins21','Del11','Del21'

    :param df: output of JACUSA2 Call2
    :return: a dataframe of columns representing features: Reference, Base, coverage from both conditions,
    min coverage, Base frequencies, Mismatch, deletion and insertion scores and counts
    """ 
    # select features from the JACUSA2 CALL 2 output
    df=df[df['strand']=='+']
    df2= df.rename(columns={"#contig": "Ref", "end": "Pos", "score":"Mis", 'ref':'Base'})
    df2 =  df2[df2['Mis']!='*']
    df2['Mis'] =  df2['Mis'].astype(float)
    tmp=np.zeros((df2.shape[0],6))
    df2[['A2','C2','G2','T2']] = df['bases21'].str.split(',',expand=True).astype(float)
    df2["Cov2"]= df2["A2"] +df2["C2"]+df2["G2"]+df2["T2"]
    df2[['A1','C1','G1','T1']] = df['bases11'].str.split(',',expand=True).astype(float)
    df2["Cov1"]= df2["A1"] +df2["C1"]+df2["G1"]+df2["T1"]
    df2["Min_cov"]=pd.DataFrame([df2['Cov1'], df2['Cov2']]).min()
    df2[["A1","C1","G1","T1"]]=df2[["A1","C1","G1","T1"]].div(df2.Cov1, axis=0)
    df2[["A2","C2","G2","T2"]]=df2[["A2","C2","G2","T2"]].div(df2.Cov2, axis=0)
    df2[["dA","dC","dG","dT"]] = abs(df2[["A1","C1","G1","T1"]].values - df2[["A2","C2","G2","T2"]].values)

    # extract features from info field
    for idx in range(df2.shape[0]):
        strg = df2['info'].iloc[idx].replace(';resetP=0','')
        if strg.find("deletion_score") != -1:
                interm = strg[strg.find("deletion_score=")+len("deletion_score="):]
                if interm.find(";") != -1:
                    tmp[idx][0] = float(interm[0 : interm.find(";")])
                else:
                    tmp[idx][0] = float(interm)
        if strg.find("insertion_score") != -1:
                interm = strg[strg.find("insertion_score=")+len("insertion_score="):]
                if interm.find(";") != -1:
                    tmp[idx][1] = float(interm[0 : interm.find(";")])
                else:
                    tmp[idx][1] = float(interm)
        if strg.find("ins11") != -1:
                interm = strg[strg.find("ins11=")+len("ins11="):]
                tmp[idx][2] = float(interm[0 : interm.find(",")])

        if strg.find("ins21") != -1:
                interm = strg[strg.find("ins21=")+len("ins21="):]
                tmp[idx][3] = float(interm[0 : interm.find(",")])

        if strg.find("del11") != -1:
                interm = strg[strg.find("del11=")+len("del11="):]
                tmp[idx][4] = float(interm[0 : interm.find(",")])

        if strg.find("del21") != -1:
                interm = strg[strg.find("del21=")+len("del21="):]
                tmp[idx][5] = float(interm[0 : interm.find(",")])

    df2.loc[:,['Del','Ins','Ins11','Ins21','Del11','Del21']]=np.nan_to_num(tmp)    
    df = df2[['Ref', 'Pos', 'Base', 'strand','Cov1','Cov2','Min_cov',"Mis","Ins","Del","dA","dC","dG","dT","A1","C1","G1","T1","A2","C2","G2","T2",'Ins11','Ins21','Del11','Del21']]
    return df

def KmerFeatures(df,K=5):
    """
    KmerFeatures function extracts features from JACUSA2 CALL2 output in Kmer context as a dataframe columns: 
    'Mis_x_$k','Ins_x_$k','Del_x_$k','5mer'

    :param df: JACUSA2 Call2 extracted features 
    :param K: mer size
    :return: a dataframe of features in Kmer context
    """     
    size = df.shape[0]
    df['5mer'] =np.zeros((size,1))
    for k in range(K):
        if k != int(K/2):
            df[['Mis_x_'+str(k),'Ins_x_'+str(k),'Del_x_'+str(k)]] =np.zeros((size,3))

    # Add Mis, Ins and Del in 5mer context
    for i in range(size):
        dff = df[df['Ref'] == df.loc[i,'Ref']]
        inter=list('NN'+ df.loc[i,'Base'] + 'NN')

        for k in range(K):
            if k != int(K/2):
                    iii = dff[dff['Pos']==dff.loc[i,'Pos']-2+k].index.values
                    if iii.size > 0 :
                        df.loc[i,'Mis_x_'+str(k)] = df.loc[iii[0],'Mis']
                        df.loc[i,'Ins_x_'+str(k)] = df.loc[iii[0],'Ins']
                        df.loc[i,'Del_x_'+str(k)] = df.loc[iii[0],'Del']
                        inter[k] = df.loc[iii[0],'Base']
        df.loc[i,'5mer']= "".join(inter)

    # Add sum of Mis, Ins, Del   
    df[['SumMis','SumIns','SumDel']]  = df[['Mis','Ins','Del']]
    for k in range(K):
        if k != int(K/2):
            df['SumMis'] = df['SumMis'] + df['Mis_x_'+str(k)]        
            df['SumIns'] = df['SumIns'] + df['Ins_x_'+str(k)]         
            df['SumDel'] = df['SumDel'] + df['Del_x_'+str(k)]         

    return df

def BarPlot(df_,feat,title,target,  indx=[],outlier = 0.001,neigh=20, path ='',method = 'IQR'):
    """
    BarPlot function perform bar plot of one features(score) from JACUSA2 CALL2 output

    :param df_: a dataframe of features including the score to plot  
    :param feat: feature's name
    :param title: title of the barplot
    :param indx: outliers indices, if provided, no need to do the outlier detection
    :param outlier: outlier contamination for the LocalOutlierFactor method
    :param zoom: provided if we want to zoom in a specific rRNA squence region 
    :param path: path to save plots as eps/pdf format
    :param method: method used to detect outliers : IQR or LOF
    
    """     
    val = df_[feat]
    pos =df_['Pos']
    # outliers indices
    if len(indx)>0:
        # outliers are given
        indices = indx
    else: 
        if method == 'LOF':
           # perform outlier detection with LocalOutlierFactor method
            lof = LocalOutlierFactor(contamination =outlier, n_neighbors = neigh)
            ohat = lof.fit_predict(val.values.reshape(-1, 1))
            indices = np.where(ohat == -1)
            outl = df_.iloc[indices[0],:]
        else: 
        # Using IQR method, extract the upper and lower quantiles
            pokemon_HP_lq = df_[feat].quantile(0.25)
            pokemon_HP_uq = df_[feat].quantile(0.75)
            #extract the inter quartile range
            pokemon_HP_iqr = pokemon_HP_uq - pokemon_HP_lq#get the upper and lower bounds
            lower_bound = pokemon_HP_lq - 1.5*pokemon_HP_iqr
            upper_bound = pokemon_HP_uq + 1.5*pokemon_HP_iqr#extract values outside these bounds 
            outl = df_[(df_[feat] <= lower_bound) | (df_[feat] >= upper_bound)]
                
    # barplots
    _, ax = plt.subplots(figsize=(20,10), dpi = 300)
    ax.stem(pos,val, 'silver',markerfmt=" ", label = 'Inliers')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Highlight each classe of outliers with a diffrent color : modfied, neighbors and unmodified positions
    ind = outl[outl['ModStatus' ] != 'Unm'].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt='.', label = 'Outlier & Modified')
                plt.setp([m, s], color='#2CBDFE')
    ind = outl[outl['Status' ].isna()].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt=".", label = 'Outlier & Neighbor')
                plt.setp([m, s], color='#47DBCD')
    ind = outl[outl['Status' ] == 'Unm'].index
    if len(ind) != 0:
                m, s, b = ax.stem(pos[ind],val[ind], markerfmt=".", label = 'Outlier & Not modified')             
                plt.setp([m, s], color='#F5B10C')
    for t in target:
                ax.annotate(str(t), (t, val[pos == t].values[0]), xytext=(0, 1), textcoords="offset points", fontweight='bold',fontstyle='italic', size = 6)          
          
    ax.legend()
    ax.set_xlabel("Position (nt)",fontsize = 16)
    ax.set_ylabel(feat,fontsize = 16)
    ax.set_title(title, fontsize = 16)
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(14)    
    ax.margins(x=0.01, y = 0.05)
    plt.savefig(path+ title+'_'+feat+'_'+ str(outlier) + '.eps') 
    plt.savefig(path+ title+'_'+feat+'_'+ str(outlier) + '.pdf') 
    
def Zoom(df_,feat, pos,outliers, s =50):
    """
    Zoom function allows to plot barplot of a specific sequence region

    :param df_: a dataframe of features including the score to plot  
    :param feat: feature's name
    :param pos: centeral position
    :param outliers: outliers indices
    :param s: the length of the sequence from both sides of the central position 
    
    """
    indx = np.where(df_['Pos'] == pos)
    df_ = df_.loc[indx[0][0]-s:indx[0][0]+s,:].reset_index()
    indx = np.searchsorted(df_['Pos'], outliers)
    BarPlot( df_,feat,title, indx,'zoom'+str(pos) )
    
    
def ScatterPlot_2Feat(df,label1,label2,feat,ref,outlier=0.001,neigh=20, path =''):
    """
    ScatterPlot_2Cond function perform scatter plot with two features (scores) from two JACUSA2 CALL2 outputs

    :param df: a dataframe of two features (scores) and position columns
    :param label1: label of the first JACUSA2 CALL2 output 
    :param label2: label of the second JACUSA2 CALL2 output   
    :param feat: feature's name
    :param ref: rRNA type 18S/28S
    :param outlier: outlier contamination for the LocalOutlierFactor method    
    :param path: path to save plots in eps/pdf format
    
    """     
    
    # select condition features
    X = df.loc[:,feat+'_x'].values.reshape(-1, 1)  
    Y = df.loc[:,feat+'_y'].values.reshape(-1, 1)  
    
    # perform linear regression
    linear_regressor = LinearRegression()  
    linear_regressor.fit(X, Y)  
    Y_pred = linear_regressor.predict(X)

    # identify outliers
    lof = LocalOutlierFactor(contamination = outlier, n_neighbors = neigh) 
    ohat = lof.fit_predict(df.loc[:,[feat+'_x',feat+'_y']])  
    indices = np.where(ohat == -1)    
    zz = np.zeros(len(X)).astype(int)
    zz[indices[0]]=int(1)
    
    # scatter plot
    fig, ax = plt.subplots(dpi = 1200)
    colors = np.array(['silver', 'lime'])
    sp = ax.scatter(X, Y,color=colors[zz])   
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.plot(X, Y_pred, color='red')
    
    # add labels to outliers
    for ii in indices[0]:
        ax.annotate(str(df.loc[ii,'Pos_x']), (X[ii],Y[ii]), xytext=(0, 1), textcoords="offset points",fontweight='bold',size= 8)

    ax.set_title(ref + '_'+feat)
    ax.set_xlabel(label1)
    ax.set_ylabel(label2)
    plt.savefig(path +ref + '_'+feat+'_'+ str(outlier) +'.eps') 
    plt.savefig(path+ref + '_'+feat+'_'+ str(outlier) +'.pdf')     
    
    
def ScatterPlot_3Feat(df,KO, feat,ref,target,Dtype, outlier=0.001,neigh=20, path =''):
    """
    ScatterPlot_3Cond function perform scatter plot with three features (scores) from three JACUSA2 CALL2 outputs

    :param df: a dataframe of two features (scores) and position columns
    :param KO: KO label 
    :param feat: feature's name
    :param ref: rRNA type 18S/28S
    :param target: KO target position  
    :param Dtype: Flongle/MinION    
    :param outlier: outlier contamination for the LocalOutlierFactor method    
    :param path: path to save plots in eps/pdf format
    
    """     
    
    # select condition features
    X = df.loc[:,feat+'_x'].values.reshape(-1, 1)  
    Y = df.loc[:,feat+'_y'].values.reshape(-1, 1)  
    Z = df[feat]    
    
    # perform linear regression
    linear_regressor = LinearRegression()  
    linear_regressor.fit(X, Y)  
    Y_pred = linear_regressor.predict(X)

    # identify outliers
    lof = LocalOutlierFactor(contamination = outlier, n_neighbors = neigh) 
    ohat = lof.fit_predict(df.loc[:,[feat+'_x',feat+'_y', feat]])  
    idx_bool = ohat == -1
    indices = np.where(ohat == -1)    
    zz = np.zeros(len(X)).astype(int)
    zz[indices[0]]=int(1)
    lofscores = pd.DataFrame({'pos': df['Pos'],'scores':lof.negative_outlier_factor_, 'outlier':idx_bool})
    lofscores.to_csv(path+feat+'_lof_scores.csv', index = False)
    
    # scatter plot
    color_map = plt.cm.get_cmap('Spectral')
    reversed_color_map = color_map.reversed()    
    fig, ax = plt.subplots(dpi = 1200)

    zz = np.zeros(len(X) + len(target)).astype(int)
    for sz in range(len(target)):
        dff = df[df['Pos'] == target[sz]].index.values
        X= np.append(X,X[dff[0]])
        Y=np.append(Y,Y[dff[0]])
        Z=np.append(Z,Z[dff[0]])
        Y_pred=np.append(Y_pred,Y_pred[dff[0]])
        zz[-sz-1]=int(1)
    
    colors = np.array(['None', 'lime'])
    sp = ax.scatter(X, Y, c= Z, cmap=reversed_color_map, edgecolor=colors[zz])   
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.colorbar(sp, label = ' rRNA WTvsKO', shrink=0.5, fraction=0.046, pad = 0.04)
    ax.plot(X, Y_pred, color='red')

    
    # add labels to outliers
    dff = df[df['Pos'].isin(target)].index.values
    for ii in indices[0]:
        if ii in dff :
            ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", fontweight='bold',fontstyle='italic', size = 6)          
        else:
            if abs(ii - dff[0]) < 3:
                ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", fontweight='bold', size = 6)
            else:
                ax.annotate(str(df.loc[ii,'Pos']), (X[ii], Y[ii]), xytext=(0, 1), textcoords="offset points", size = 6)
            
    for position in  dff:
        if (position in indices[0]) == False :
            ax.annotate(str(df.loc[position,'Pos']), (X[position], Y[position]),  xytext=(0, 1), textcoords="offset points",fontstyle='italic',size = 6)     
        
    ax.set_title(KO + '_'+Dtype+'_' + feat)
    ax.set_xlabel(' rRNA WTvsIVT')
    ax.set_ylabel(' rRNA KOvsIVT')
    plt.savefig(path+ KO+ '_'+ ref+'_'+Dtype+'_'+feat+'_'+ str(outlier) + '.eps',  bbox_inches = "tight") 
    plt.savefig(path+ KO+ '_'+ ref+'_'+Dtype+'_'+feat+'_'+ str(outlier) + '.pdf',  bbox_inches = "tight") 
   