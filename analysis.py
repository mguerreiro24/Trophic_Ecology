#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                     #
#                     Miguel Guerreiro                                #
#                          2020                                       #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
##from network_TW import *
import trophic_WEBS
from copy import deepcopy
from random_sampler import *
from random import randint,shuffle
import numpy as np
import matplotlib.pyplot as plt
from statistics import stdev as sd


def SE(l):
    """standard error of list
"""
    try:
        return sd(l)/(len(l)**(.5))
    except:
        return 0


def confidence_interval_95(l):
    return SE(l)*1.96


def select_binary(lst, indices):
    return [lst[ii] for ii,i in enumerate(indices) if i==1]


def make_binary(indices,lst):
    l = len(lst)
    g = [0]*l
    for i in indices:
        g[i]=1
    return g


def isEmpty(deleturn):
    for i in deleturn:
        if len(i)>0:
            return False
    return True


#----------------------------------------------------------------------------------------
def load_data(filename,Type):
    """Loading files to graph objects
requires: str(filename) - name title of file(tsv format, ) without its extension (.txt)
          str(Type)     - type of file: 'EL' edgelist,
                                        'AM+' Adjacency matrix:
                                                'AMU' unweighted,
                                                'AMW' weighted
ensures: graph object
"""
    if Type=='EL':
        nodes_names,edges_w,edges_coords = open_EL(filename)
    elif Type[:2]=="AM":
        nodes_names,MA = load_AM(filename,h=False,sep="\t")
        if Type[-1]=="U":
            RA = MA
        elif Type[-1]=="W":
            RA = absolute_AM_2_relative_AM(MA)
        edges_w,edges_coords=process_AM(RA)
    g = from_out_O_fmt_files(nodes_names,edges_w,edges_coords)
    return g


#--------------------------------------------------
def linear_curve(filename,maxi):
    #https://matplotlib.org/3.2.0/gallery/lines_bars_and_markers/fill_between_demo.html
    rang = list(range(1,maxi))
    X,x,xx,Y,Z,W = [[],[],[],[],[],[]]
    Am = []
    AE = []
    AS = []
    Em = []
    EE = []
    ES = []
    for i,f in [("{}_{}.txt".format(filename,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        A = sorted([float(i[0]) for i in j])
        E = sorted([float(i[1]) for i in j])
        s = round((len(j)-1)*.025)
        e = round((len(j)-1)*.975)
        Am.append(sum(A)/len(A))
        AS.append(A[s])
        AE.append(A[e])
        
        Em.append(sum(E)/len(E))
        ES.append(E[s])
        EE.append(E[e])

        X.extend([float(f) for i in range(len(A))])
        x.extend([float(f) for i in range(len(E))])
        Y.extend(A)
        Z.extend(E)
    X = np.array(X)
    x = np.array(x)
    Y = np.array(Y)
    Z = np.array(Z)

    title = 'results'
    rang = np.array(rang)
    Am = np.array(Am)
    AE = np.array(AE)
    AS = np.array(AS)
    Em = np.array(Em)
    EE = np.array(EE)
    ES = np.array(ES)
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel("Percentage of Diet (%)", fontsize = 12)
    ax1.set_ylabel("Percentage (%)", fontsize = 12)
    colors = ["xkcd:"+x for x in ["windows blue", "amber", "coral", "blood red", "tealish green", "vermillion", "very dark blue", "light grey blue", "chestnut", "pale gold", "deep sea blue", "rose red"]][:3]
    ax1.plot(rang, Am, '-', color=colors[0],label='Attack sensitivity')
    ax1.fill_between(rang, AS, AE, alpha=0.2, color=colors[0])
    ax1.scatter(X, Y,c = colors[0], alpha=0.02)

    ax1.plot(rang, Em, '-', color=colors[1],label='Error sensitivity')
    ax1.fill_between(rang, ES, EE, alpha=0.2, color=colors[1])
    ax1.scatter(x, Z,c = colors[1], alpha=0.02)

    ax1.legend(loc = 'upper left')
    fig.tight_layout()
    fig.savefig('{}_{}_CI95.png'.format(filename,title),dpi = (300))
    plt.close(fig='all')


def compare_curve(filename1,filename2,maxi,variable="Edges",index=2):
    rang = list(range(1,maxi))
    X,x,xx,Y,Z,W = [[],[],[],[],[],[]]
    Am = []
    AE = []
    AS = []
    
    Em = []
    EE = []
    ES = []
    for i,f in [("{}_{}.txt".format(filename1,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        s = round((len(j)-1)*.025)
        e = round((len(j)-1)*.975)
        A = sorted([float(i[index]) for i in j])
        Am.append(sum(A)/len(A))
        AE.append(A[e])
        AS.append(A[s])

        X.extend([float(f) for i in range(len(A))])
        Y.extend(A)

    for i,f in [("{}_{}.txt".format(filename2,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        s = round((len(j)-1)*.025)
        e = round((len(j)-1)*.975)
        E = sorted([float(i[index]) for i in j])
        Em.append(sum(E)/len(E))
        EE.append(E[e])
        ES.append(E[s])
        
        x.extend([float(f) for i in range(len(E))])
        Z.extend(E)
    X = np.array(X)
    x = np.array(x)
    Y = np.array(Y)
    Z = np.array(Z)

    title = 'comparison'
    rang = np.array(rang)
    Am = np.array(Am)
    AE = np.array(AE)
    AS = np.array(AS)
    
    Em = np.array(Em)
    EE = np.array(EE)
    ES = np.array(ES)

    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel("Percentage of Diet (%)", fontsize = 12)
    ax1.set_ylabel(variable, fontsize = 12)
    colors = ["xkcd:"+x for x in ["windows blue", "amber", "coral", "blood red", "tealish green", "vermillion", "very dark blue", "light grey blue", "chestnut", "pale gold", "deep sea blue", "rose red"]]

    if index==2:
        cc = colors[2]
        ccc = colors[3]
    elif index==0:
        cc = colors[-5]
        ccc = colors[-2]
    elif index==1:
        cc = colors[-3]
        ccc = colors[1]
    
    ax1.plot(rang, Am, '-', color=cc,label='Allesina')
    ax1.fill_between(rang, AS, AE, alpha=0.2, color=cc)
    ax1.scatter(X, Y,c = cc, alpha=0.02)


    ax1.plot(rang, Em, '-', color=ccc,label='Guerreiro & Scotti')
    ax1.fill_between(rang, ES, EE, alpha=0.2, color=ccc)
    ax1.scatter(x, Z,c = ccc, alpha=0.02)
    
    ax1.legend(loc = 'upper left')
    fig.tight_layout()
    fig.savefig('{}_{}_{}CI95.png'.format(filename2,variable,title),dpi = (300))
    plt.close(fig='all')


if __name__=="__main__":
    import os


    home_path = os.getcwd()
    file_folder = os.path.join(home_path,"weighted_FOODWEBS")
    for folder in os.listdir(file_folder)[:]:
        #folder = "ESM2---Olmo_Gilabert_2019"
        print(folder)
        filename = os.path.join(os.path.join(file_folder,folder),"data")
        if folder=="ESM2---Olmo_Gilabert_2019":
            Type = "EL"

        else:
            Type = "AM"
        #Allesina et al. (2006) approach
        print("...old")
##        trophic_WEBS.Allesina(filename,Type)
        print(" graph")
        filenameg1 = os.path.join(os.path.join(file_folder,folder),"data_r")
        linear_curve(filenameg1,99)
        #new approach
        print("new...")
##        trophic_WEBS.GuerreiroScotti(filename,Type)#1000
        print(" graph")
        filenameg2 = os.path.join(os.path.join(file_folder,folder),"data_2")
        linear_curve(filenameg2,99)
        #Allesina vs new approach in terms of
        #edges removed
        compare_curve(filenameg1,filenameg2,99)
        #attack
        compare_curve(filenameg1,filenameg2,99,'Attack Sensitivity',0)
        #error
        compare_curve(filenameg1,filenameg2,99,'Error Sensitivity',1)
#    #multiple nutrients weighted network (multidimensional)
    r = input("press 'enter' to close")





























