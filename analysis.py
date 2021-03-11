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


#------------------------------------------------------
def Build_DO_n(**kwargs):
    """ builds a Dominator matrix secondary extinctions by node/vertex
    requires:
        DM - dominator matrix [list(list())]
    optional:
        N - total number of nodes/vertices present
        out_giant - total number of isolated nodes/vertices
    ensure: list(number of directly dependent nodes for each node)"""
    DM = kwargs["DM"]
    if "N" in kwargs:
        N = kwargs["N"]
    else:
        N = len(DM)-1# -1 -->root?
    if "out_giant" in kwargs:
        to_add = kwargs["out_giant"]
    else:
        to_add = 0
    length_matrix = len(DM)
    
    DO_n = [(to_add + sum([DM[f][c] for f in range(length_matrix)])-1)/N for c in range(1,length_matrix)]
    return DO_n


def ES(DO_n):
    """Error sensitivity:
    the average of that(the fraction of nodes
    dominated by each node) corresponds to the error sensitivity
requires: list() of fractions of nodes dominated by each node
ensures: Value of ES
"""
    return sum(DO_n)/len(DO_n)#-1)


def AS(DO_n):
    """Attack sensitivity:
    attack sensitivity is the largest fraction of secondary
    extinctions triggered by a primary extinction
requires: list() of fractions of nodes dominated by each node
ensures:Value of AS
"""
    return max(DO_n)


def sensitivity_test(DM):
    """
requires: dominator matrix[list(list())]
ensures: attack and error sensitivities
"""
    N = len(DM)
    DO_n = [(sum([DM[f][c] for f in range(N)])-1)/(N-1) for c in range(1,N)]#N>len(DM) ? fill with zer0s?
    #DO_n = [(sum([DM[f][c] for f in range(N)])-1)/(N_new-1) for c in range(1,N)]
    attack = AS(DO_n)
    error = ES(DO_n)
    return attack,error


#----------------------------------------------------------------------------------------
def picking_combos(length, amount,weights,threshold,threshold_step):
    """makes a binary list() with combinations that are below or equal to 
threshold"""
#sorted([i for i in binary_combos(len(for_combo_maker)) if sum(i)!=0 and sum(select_binary(for_combo_maker,i))<=threshold/100 and sum(select_binary(for_combo_maker,i))>=(threshold-threshold_step)/100],key=lambda x:sum(select_binary(for_combo_maker,x)),reverse=True)
    if length==0:
        d = [[]]*amount
        return d
    elif length==1:
        d = [[1]]*amount
        return d
    elif sum(weights)<=(threshold)/100:
        #print("sum all*!")
        d = [[1]*length]*amount
        return list(d)
    else:#error prone here
        d = [[0] for i in range(length)]
        n_success = 0
        out = []
        tries = 0
        while n_success<amount:
            if tries/(n_success+1)>10000:
                break
            iii = WithoutReposition(randint(0,length),list(range(length)))
            deletes = make_binary(iii,d)
            #following combo is not found
            if sum(select_binary(weights,deletes))<=threshold/100 and sum(select_binary(weights,deletes))>(threshold-threshold_step)/100:
                n_success += 1
                out.append(deletes)
                tries = 0
            tries +=1
        if tries/(n_success+1)>10000:
            #print("little success")
            d = [[0] for i in range(length)]
            n_success = 0
            out_i = []
            while n_success<amount:
                iii = WithReposition(1,list(range(length)))
                deletes = make_binary(iii,d)
                out_i.append(deletes)
                n_success += 1
            return out_i
        else:
            #print("found combos")
            return out


def percent_COMBO(g,
                  threshold_step=10,
                  threshold=10,
                  monte_carlo=20):
    """
ensures: files with attack and error sensitivity
"""
    root = 0
    deletions_t = []
    wk = g.es.select(weights_le=threshold/100)
    #picking up combinations of edges that are between previous threshold and current
    for predator_i in range(1,len(g.vs)):
        eating_e = [(e.index,e["weights"]) for e in wk.select(_to=predator_i)]
        print("{} % threshold| #{} component | {} edges to pick".format(threshold,predator_i,len(eating_e)))
        for_combo_maker = [i[1] for i in eating_e]
        eds = [i[0] for i in eating_e]
        deletions_b = picking_combos(len(for_combo_maker),monte_carlo,for_combo_maker,threshold,threshold_step)
        deletions_t.append([select_binary(eds, ds) for ds in deletions_b])

    #setting up file for output
    with open("{}_{}.txt".format(filename,threshold),"w") as handle:
        handle.write("")

    #running the monte-carlo sequence of trials
    for turn in range(monte_carlo):
        gg = deepcopy(g)
        delete_turn = [d[turn] for d in deletions_t]
        for_del = [item for sublist in delete_turn for item in sublist]#edges are unique!
        print("[Threshold:{}%]Running Monte-Carlo cycle #{} edges removed:{}".format(threshold,turn,len(for_del)))
        gg.delete_edges(for_del)
        nodes_names,edges_w,edges_coords = process_iGraph(gg)
        m = EL_2_matrix(edges_coords=edges_coords, edges_w=edges_w)#already rooted
        DM = DOM_MATRIX(m)
        #printing out the results
        A,E = sensitivity_test(DM)#<---store values for analysis
        with open("{}_{}.txt".format(filename,threshold),"a") as handle:
            handle.write("{}\t{}\n".format(A,E))

            

#----------------------------------------------------------------------------------------
def Diet_percent_edgeRemoval_analysis(filename='ESM2---Olmo_Gilabert_2019',
                                      threshold_step=10,
                                      threshold=10,
                                      threshold_max=20,
                                      Type='EL'):
    """calculates the Attack and Error sensitivity of dominator matrix, by removing edges that are below a threshold.


                        Allesina et al. (2006) approach
"""
    g = load_data(filename,Type)
    nodes_names,edges_w,edges_coords = process_iGraph(g)
    while threshold<=threshold_max:
        gg = deepcopy(g)
        deletions_t = []
        wk = gg.es.select(weights_le=threshold/100)
        for predator_i in range(len(nodes_names)):
            eating_e = [(e.index,e["weights"]) for e in wk.select(_to=predator_i)]
            eds = [i[0] for i in eating_e]
            deletions_t.extend(eds)
        delete_turn = deletions_t
        gg.delete_edges(delete_turn)
        nodes_names,edges_w,edges_coords = process_iGraph(gg)
        m = EL_2_matrix(edges_coords=edges_coords, edges_w=edges_w)
        matrix_root,r_nodes_names,r_coords,r_w = RootFW(m,nodes_names,edges_coords,edges_w)
        DM = DOM_MATRIX(matrix_root)
        A,E = sensitivity_test(DM)
        with open("{}_r_{}.txt".format(filename,threshold),"w") as handle:
            handle.write("{}\t{}\n".format(A,E))
        #print("Attack sensitivity:{}\nError sensitivity:{}".format(A,E))#"%0.2f" % 
        threshold += threshold_step


def Allesina(filename='ESM2---Olmo_Gilabert_2019',
             threshold_step=10,
             threshold=10,
             threshold_max=20,
             Type='EL'):
    gg = load_data(filename,Type)
    g = RootingGraph(gg)
    n_nodes = len(g.vs)#keep root -> basal nodes constant
    n_edges = len(g.es)
    basal_e = g.es.select(_source=0)#primary producers [_source==0(root)] edges
    basal = [edge.target for edge in basal_e]#primary producers vertice
    predators = [j for j in range(1,n_nodes) if j not in basal]
    while threshold<=threshold_max:
        with open("{}_r_{}.txt".format(filename,threshold),"w") as handle:
            handle.write("")
        gg = deepcopy(g)#<--
        deletions_t = []
        wk = gg.es.select(weights_le=threshold/100)
        for predator_i in predators:
            #with links bellow threshold
            eating_e = [(e.index,e["weights"]) for e in wk.select(_to=predator_i)]
            eds = [i[0] for i in eating_e]
            deletions_t.extend(eds)
        gg.delete_edges(deletions_t)

        length2delete = len(deletions_t)#new info for statistics
        #remove nodes that are disconnected by selecting only the giant component
        gg = gg.components().giant()

        #back to dominator matrix
        nodes_names,edges_w,edges_coords = process_iGraph(gg)
        if len(edges_coords)!=0 and len(edges_w)!=0:
            m = EL_2_matrix(edges_coords=edges_coords, edges_w=edges_w)#solve here
            DM = DOM_MATRIX(m)#DM = DOM_MATRIX(matrix_root)
            #keep original number of nodes for calculation of sensitivity
            isolated_nodes = n_nodes - len(DM)
            second_direct_extinctions2 = Build_DO_n(DM=DM, N=n_nodes, out_giant=isolated_nodes)
            A2 = AS(second_direct_extinctions2)
            E2 = ES(second_direct_extinctions2)
            with open("{}_r_{}.txt".format(filename,threshold),"w") as handle:
                handle.write("{}\t{}\t{}\n".format(A2,E2,length2delete))#/n_edges))
        else:
            with open("{}_r_{}.txt".format(filename,threshold),"w") as handle:
                handle.write("{}\t{}\t{}\n".format(1.0,1.0,length2delete))
        threshold += threshold_step


def Diet_percent_combination_edgeRemoval_analysis(filename='ESM2---Olmo_Gilabert_2019',
                                                  threshold_step=10,
                                                  threshold=10,
                                                  threshold_max=20,
                                                  monte_carlo=20,
                                                  Type='EL'):#add threading input number
    """calculates the Attack and Error sensitivity of dominator matrix, by removing combinations of in edges
that combined sum below a threshold.
"""
    gg = load_data(filename,Type)
    g = RootingGraph(gg)
    while threshold<=threshold_max:
        percent_COMBO(g,threshold_step,threshold,monte_carlo)#threshold_max
        threshold += threshold_step


def Diet_edgeRemoval_analysis_new(filename='ESM2---Olmo_Gilabert_2019',
                                                  threshold_step=1,
                                                  threshold=1,
                                                  threshold_max=99,
                                                  monte_carlo=33,
                                                  Type='EL'):#add threading input number
    """calculates the Attack and Error sensitivity of dominator matrix, by removing combinations of in edges
that combined sum equal or higher than a threshold.
"""
    gg = load_data(filename,Type)
    g = RootingGraph(gg)
    n_nodes = len(g.vs)#keep root -> basal nodes constant
    n_edges = len(g.es)
    basal_e = g.es.select(_source=0)#primary producers [_source==0(root)] edges
    basal = [edge.target for edge in basal_e]#primary producers vertice
    predators = [j for j in range(1,n_nodes) if j not in basal]
    while threshold<=threshold_max:
        with open("{}_2_{}.txt".format(filename,threshold),"w") as handle:
            handle.write("")
        delete_turn = []
        for predator_i in predators:#changed algorithm: ranked(sorted) links, cutted until threshold reached.
            #with links bellow threshold
            toEating = [(e.index,e["weights"]) for e in g.es.select(_to=predator_i) if e["weights"]<=(threshold/100)]
            to_predator = []
            for _ in range(monte_carlo):
                shuffle(toEating)#greedy approach, changed what happens when weight ties| --> shuffle
                eating_e = sorted(toEating,key=lambda x:x[1],reverse=True)
                place_holder = []
                deletions_t = []
                for i,w in eating_e:
                    deletions_t.append(w)
                    place_holder.append(i)
                    if sum(deletions_t)>=threshold/100:
                        break
                to_predator.append(place_holder)
            delete_turn.append(to_predator)
        #print(delete_turn)
        for turn in range(monte_carlo):
            #get delete edges for each monte-carlo turn
            gg = deepcopy(g)
            delete_row = ()
            if not isEmpty(delete_turn):
                delete_row = [j for sublist in [i[turn] for i in delete_turn if i!=[]] for j in sublist]#
                gg.delete_edges(delete_row)
                #remove nodes that are disconnected by selecting only the giant component
                gg = gg.components().giant()
            length2delete = len(delete_row)#new info for statistics
            #back to dominator matrix
            nodes_names,edges_w,edges_coords = process_iGraph(gg)
            if len(edges_coords)!=0 and len(edges_w)!=0:
                m = EL_2_matrix(edges_coords=edges_coords, edges_w=edges_w)#solve here
                DM = DOM_MATRIX(m)#DM = DOM_MATRIX(matrix_root)

                #keep original number of nodes for calculation of sensitivity
                isolated_nodes = n_nodes - len(DM)

                #1 calc
##                second_direct_extinctions1 = Build_DO_n(DM=DM, N=n_nodes)
##                A1 = AS(second_direct_extinctions1)
##                E1 = ES(second_direct_extinctions1)
##                with open("{}_1_{}.txt".format(filename,threshold),"a") as handle:
##                    handle.write("{}\t{}\n".format(A1,E1))
                #2 calc
                second_direct_extinctions2 = Build_DO_n(DM=DM, N=n_nodes, out_giant=isolated_nodes)
                A2 = AS(second_direct_extinctions2)
                E2 = ES(second_direct_extinctions2)
                with open("{}_2_{}.txt".format(filename,threshold),"a") as handle:
                    handle.write("{}\t{}\t{}\n".format(A2,E2,length2delete))#/n_edges))
                #3 calc
##                second_direct_extinctions3 = Build_DO_n(DM=DM)
##                A3 = AS(second_direct_extinctions3)
##                E3 = ES(second_direct_extinctions3)
##                with open("{}_3_{}.txt".format(filename,threshold),"a") as handle:
##                    handle.write("{}\t{}\n".format(A3,E3))
            else:
                #1 calc
##                with open("{}_1_{}.txt".format(filename,threshold),"a") as handle:
##                    handle.write("{}\t{}\n".format(1.0,1.0))
                #2 calc
                with open("{}_2_{}.txt".format(filename,threshold),"a") as handle:
                    handle.write("{}\t{}\t{}\n".format(1.0,1.0,length2delete))
                #3 calc
##                with open("{}_3_{}.txt".format(filename,threshold),"a") as handle:
##                    handle.write("{}\t{}\n".format(1.0,1.0))
        threshold += threshold_step


#--------------------------------------------------
def linear_curve(filename,maxi):
    #https://matplotlib.org/3.2.0/gallery/lines_bars_and_markers/fill_between_demo.html
    rang = list(range(1,maxi))
    X,x,xx,Y,Z,W = [[],[],[],[],[],[]]
    Am = []
    AE = []
    Em = []
    EE = []
##    Dm = []
##    DE = []
    for i,f in [("{}_{}.txt".format(filename,f),f) for f in rang]:#
        with open(i,"r") as handle:
            j = [i.split("\t") for i in handle.readlines()]
        A = [float(i[0]) for i in j]
        E = [float(i[1]) for i in j]
##        D = [float(i[2]) for i in j]
        Am.append(sum(A)/len(A))
        AE.append(confidence_interval_95(A))
        Em.append(sum(E)/len(E))
        EE.append(confidence_interval_95(E))
##        Dm.append(sum(D)/len(D))
##        DE.append(confidence_interval_95(D))

        X.extend([float(f) for i in range(len(A))])
        x.extend([float(f) for i in range(len(E))])
##        xx.extend([float(f) for i in range(len(D))])
        Y.extend(A)
        Z.extend(E)
##        W.extend(D)
    X = np.array(X)
    x = np.array(x)
##    xx = np.array(xx)
    Y = np.array(Y)
    Z = np.array(Z)
##    W = np.array(W)

    title = 'results'
    rang = np.array(rang)
    Am = np.array(Am)
    AE = np.array(AE)
    Em = np.array(Em)
    EE = np.array(EE)
##    Dm = np.array(Dm)
##    DE = np.array(DE)
    
    fig = plt.figure(figsize=(12, 12))
    ax1 = fig.add_subplot(1,1,1)
##    ax.set_title("[{}]".format(title), fontsize = 22)
    ax1.set_xlabel("Percentage of Diet (%)", fontsize = 12)
    ax1.set_ylabel("Percentage (%)", fontsize = 12)
    colors = ["xkcd:"+x for x in ["windows blue", "amber", "coral", "blood red", "tealish green", "vermillion", "very dark blue", "light grey blue", "chestnut", "pale gold", "deep sea blue", "rose red"]][:3]
    ax1.plot(rang, Am, '-', color=colors[0],label='Attack sensitivity')
    ax1.fill_between(rang, Am - AE, Am + AE, alpha=0.2, color=colors[0])
    ax1.scatter(X, Y,c = colors[0], alpha=0.02) # actual
    #ax1.legend(loc = 'lower right')
    #ax2 = ax1.twinx()
    #ax2.set_ylabel("Error sensitivity", fontsize = 12)
    ax1.plot(rang, Em, '-', color=colors[1],label='Error sensitivity')
    ax1.fill_between(rang, Em - EE, Em + EE, alpha=0.2, color=colors[1])
    ax1.scatter(x, Z,c = colors[1], alpha=0.02) # actual
##    print(len(rang),len(Dm))
####    ax1.plot(rang, Dm, '-', color=colors[2],label='Edges removed')
####    ax1.fill_between(rang, Dm - DE, Dm + DE, alpha=0.2, color=colors[2])
####    ax1.scatter(xx, W,c = colors[2], alpha=0.02) # actual
    
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
##    import os
##    home_path = os.getcwd()
##    folder = os.path.join(os.path.join(home_path,"weighted_FOODWEBS"),"chesa")
##    filename = os.path.join(folder,"data")
##    Type = 'AMW'
##    threshold_step = 1
##    threshold = 1
##    threshold_max = 20
##    monte_carlo = 10
##    Allesina(filename,threshold_step,threshold,threshold_max,Type)
    import os
    threshold_step = 1
    threshold = 1
    threshold_max = 99
    monte_carlo = 100


    home_path = os.getcwd()
    file_folder = os.path.join(home_path,"weighted_FOODWEBS")
    for folder in os.listdir(file_folder)[:]:
        #folder = os.listdir(file_folder)[3]#"ESM2---Olmo_Gilabert_2019"
        print(folder)#indent here on
        filename = os.path.join(os.path.join(file_folder,folder),"data")
        if folder=="ESM2---Olmo_Gilabert_2019":
            Type = 'EL'
        else:
            Type = 'AMW'
        #get basic stats
        print("basic stats")
        basic_stats(filename,Type)
        #new approach
        print("new...")
        Diet_edgeRemoval_analysis_new(filename,threshold_step,threshold,threshold_max,monte_carlo,Type)
        print(" graph")
        filenameg2 = os.path.join(os.path.join(file_folder,folder),"data_2")
        linear_curve(filenameg2,99)

        #Allesina et al. (2006) approach
        print("...old")
        Allesina(filename,threshold_step,threshold,threshold_max,Type)
        print(" graph")
        filenameg1 = os.path.join(os.path.join(file_folder,folder),"data_r")
        linear_curve(filenameg1,99)
        #Allesina vs new approach in terms of edges removed
        compare_curve(filenameg1,filenameg2,99)
#    #multiple extinctions approach

#    #multiple nutrients weighted network (multidimensional)
    r = input("press 'enter' to close")

































