#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                     #
#                     Miguel Guerreiro                                #
#                          2020                                       #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
from network_TW import *
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
    Em = []
    EE = []
    for i,f in [("{}_{}.txt".format(filename,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        A = [float(i[0]) for i in j]
        E = [float(i[1]) for i in j]
        Am.append(sum(A)/len(A))
        AE.append(confidence_interval_95(A))
        Em.append(sum(E)/len(E))
        EE.append(confidence_interval_95(E))

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
    Em = np.array(Em)
    EE = np.array(EE)
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel("Percentage of Diet (%)", fontsize = 12)
    ax1.set_ylabel("Percentage (%)", fontsize = 12)
    colors = ["xkcd:"+x for x in ["windows blue", "amber", "coral", "blood red", "tealish green", "vermillion", "very dark blue", "light grey blue", "chestnut", "pale gold", "deep sea blue", "rose red"]][:3]
    ax1.plot(rang, Am, '-', color=colors[0],label='Attack sensitivity')
    ax1.fill_between(rang, Am - AE, Am + AE, alpha=0.2, color=colors[0])
    ax1.scatter(X, Y,c = colors[0], alpha=0.02) # actual

    ax1.plot(rang, Em, '-', color=colors[1],label='Error sensitivity')
    ax1.fill_between(rang, Em - EE, Em + EE, alpha=0.2, color=colors[1])
    ax1.scatter(x, Z,c = colors[1], alpha=0.02) # actual

    ax1.legend(loc = 'upper left')
    fig.tight_layout()
    fig.savefig('{}_{}_CI95.png'.format(filename,title),dpi = (300))
    plt.close(fig='all')


def basic_stats(filename,Type):
    g = load_data(filename,Type)
    dds = g.degree(mode="in")
    d = {i:dds.count(i) for i in set(dds)}
    ww = [int(i*100) for i in g.es["weights"]]
    w = {i:ww.count(i) for i in set(ww)}

    rang = np.array(list(range(0,100)))
    weights = np.array([w[i] if i in w else 0 for i in list(range(0,100))])
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_ylabel("edges", fontsize = 12)
    ax1.set_xlabel("Percentage of diet (%)", fontsize = 12)
    ax1.bar(rang, weights, color = "xkcd:windows blue")
    fig.savefig('{}_basic_stat_weights.png'.format(filename),dpi = (300))
    plt.close(fig='all')

    rang = np.array(list(range(1,max(d.keys())+1)))
    degrees = np.array([d[i] if i in d else 0 for i in list(range(1,max(d.keys())+1))])
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_ylabel("edges", fontsize = 12)
    ax1.set_xlabel("in-degrees", fontsize = 12)
    ax1.bar(rang, degrees, color = "xkcd:amber")
    fig.savefig('{}_basic_stat_degrees.png'.format(filename),dpi = (300))
    plt.close(fig='all')


def compare_curve(filename1,filename2,maxi):
    rang = list(range(1,maxi))
    X,x,xx,Y,Z,W = [[],[],[],[],[],[]]
    Am = []
    AE = []
    Em = []
    EE = []
    for i,f in [("{}_{}.txt".format(filename1,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        A = [float(i[2]) for i in j]
        Am.append(sum(A)/len(A))
        AE.append(confidence_interval_95(A))
        X.extend([float(f) for i in range(len(A))])
        Y.extend(A)

    for i,f in [("{}_{}.txt".format(filename2,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        E = [float(i[2]) for i in j]
        Em.append(sum(E)/len(E))
        EE.append(confidence_interval_95(E))
        x.extend([float(f) for i in range(len(E))])
        Z.extend(E)
    X = np.array(X)
    x = np.array(x)
    Y = np.array(Y)
    Z = np.array(Z)

    title = 'edgesRemoved'
    rang = np.array(rang)
    Am = np.array(Am)
    AE = np.array(AE)
    Em = np.array(Em)
    EE = np.array(EE)
    
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
    ax1.set_xlabel("Percentage of Diet (%)", fontsize = 12)
    ax1.set_ylabel("edges", fontsize = 12)
    colors = ["xkcd:"+x for x in ["windows blue", "amber", "coral", "blood red", "tealish green", "vermillion", "very dark blue", "light grey blue", "chestnut", "pale gold", "deep sea blue", "rose red"]]
    ax1.plot(rang, Am, '-', color=colors[2],label='Allesina')
    #print(len(Am),len(AE),len(rang))
    ax1.fill_between(rang, Am - AE, Am + AE, alpha=0.2, color=colors[2])
    ax1.scatter(X, Y,c = colors[2], alpha=0.02) # actual


    ax1.plot(rang, Em, '-', color=colors[3],label='Guerreiro & Scotti')
    #print(len(Em),len(EE),len(rang))
    ax1.fill_between(rang, Em - EE, Em + EE, alpha=0.2, color=colors[3])
    ax1.scatter(x, Z,c = colors[3], alpha=0.02) # actual
    
    ax1.legend(loc = 'upper left')
    fig.tight_layout()
    fig.savefig('{}_{}_CI95.png'.format(filename2,title),dpi = (300))
    plt.close(fig='all')


if __name__=="__main__":
    import os
    import trophic_WEBS


    home_path = os.getcwd()
    file_folder = os.path.join(home_path,"weighted_FOODWEBS")
    for folder in os.listdir(file_folder)[:]:
        print(folder)
        filename = os.path.join(os.path.join(file_folder,folder),"data")
        if folder=="ESM2---Olmo_Gilabert_2019":
            nodes_names,edges_w,edges_coords = trophic_WEBS.open_EL(filename)
            matrix = trophic_WEBS.EL_2_matrix(edges_w=edges_w,edges_coords=edges_coords)
        else:
            matrix = trophic_WEBS.absolute_AM_2_relative_AM(trophic_WEBS.open_AM(filename))
##        #get basic stats
##        print("basic stats")
##        basic_stats(filename,Type)
        #new approach
        print("new...")
        A,E,D = trophic_WEBS.GuerreiroScotti(matrix)
        for i,_ in enumerate(A):
            with open('{}_2_{}.txt'.format(filename,i),'w') as handle:
                handle.write('')
            with open('{}_2_{}.txt'.format(filename,i),'a') as handle:
                for ii,_ in enumerate(A[i]):
                    handle.write('{}\t{}\t{}\n'.format(A[i][ii],E[i][ii],D[i][ii]))
        print(" graph")
        filenameg2 = os.path.join(os.path.join(file_folder,folder),"data_2")
        linear_curve(filenameg2,99)

        #Allesina et al. (2006) approach
        print("...old")
        A,E,D = trophic_WEBS.Allesina(matrix)
        for i,_ in enumerate(A):
            with open('{}_r_{}.txt'.format(filename,i),'w') as handle:
                handle.write('')
            with open('{}_r_{}.txt'.format(filename,i),'a') as handle:
                handle.write('{}\t{}\t{}\n'.format(A[i],E[i],D[i]))
        print(" graph")
        filenameg1 = os.path.join(os.path.join(file_folder,folder),"data_r")
        linear_curve(filenameg1,99)
        #Allesina vs new approach in terms of edges removed
        compare_curve(filenameg1,filenameg2,99)
#    #multiple extinctions approach

#    #multiple nutrients weighted network (multidimensional)
    r = input("press 'enter' to close")





























