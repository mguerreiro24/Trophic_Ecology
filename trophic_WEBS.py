#############################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++++++++++++++++++Built by Miguel Fernandes Guerreiro+++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++09/04/2021++++++++++++++++++++++++++++++++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#############################################################################
#imports
from copy import deepcopy
import numpy as np
from random import randint,shuffle

#functions
#reading
def process_EL(edges_fl):
    edges_en = [item for sublist in edges_fl for item in sublist[:-1]]
    nodes_names = list(dict.fromkeys(edges_en))
    edges_w = [float(w[-1]) for w in edges_fl]
    edges_coords = []
    for edges in edges_fl:
        count = [0,0]
        cc =[None,None]
        for n,node in enumerate(nodes_names):
            if edges[0]==node:
                cc[0] = n
                count[0] = 1
            if edges[1]==node:
                cc[1] = n
                count[1] = 1
            if count==[1,1]:
                edges_coords.append((cc[0],cc[1]))
                break
    return nodes_names,edges_w,edges_coords


def open_EL(filename, sep="\t"):
    """get edge list out of format and get values for graph
"""
    with open('{}.txt'.format(filename), "r") as handle:
        edges_fl = [i[:-1].split(sep) for i in handle.readlines()]
    nodes_names,edges_w,edges_coords = process_EL(edges_fl)
    return nodes_names,edges_w,edges_coords


def open_AM(filename,sep="\t"):
    with open('{}.txt'.format(filename),"r") as handle:
        matrix = [[float(ii) for ii in i[:-1].split(sep)] for i in handle.readlines()]
    return matrix


#processing
def transpose(l_l):
    return list(map(list, zip(*l_l)))#transpose


def matrix_2_EL(matrix):
    """receives matrix (no headers) and exports EL(edges_coords[i],edges_w[i])
"""
    EL = []
    for r in range(len(matrix)):
        for c in range(len(matrix[0])):
            if float(matrix[r][c])>0:
                EL.append((r,c,float(matrix[r][c])))#<--the format
    return EL


def EL_2_matrix(**kwargs):
    """receives EL(edges_coords[i],edges_w[i]) and exports matrix (no headers)
"""
    if "edge_list" in kwargs:#EL(edges_coords[i],edges_w[i])
        EL = kwargs["edge_list"]
        length = max(list(dict.fromkeys([item for sublist in EL for item in sublist[:-1]])))
        length += 1
        m = [[0]*length for i in range(length)]
        for i in EL:
            m[i[0]][i[1]] = i[2]
    elif "edges_coords" in kwargs and "edges_w" in kwargs:
        edges_coords = kwargs["edges_coords"]
        edges_w = kwargs["edges_w"]
        length = max(list(dict.fromkeys([item for sublist in edges_coords for item in sublist])))
        length += 1
        m = [[0]*length for i in range(length)]
        for i,w in enumerate(edges_w):
            m[edges_coords[i][0]][edges_coords[i][1]] = w
    else:
        raise Exception('Provide proper input')
    return m


def absolute_AM_2_relative_AM(absolute_AM):
    """
"""
    N = len(absolute_AM)
    totals = [sum([absolute_AM[c][r] for c in range(N)]) for r in range(N)]
    relative_M = [[0 for c in range(N)] for r in range(N)]
    for r,raw in enumerate(absolute_AM):
        for c,caw in enumerate(raw):
            if caw==0 and totals[c]==0:
                relative_M[r][c] = caw
            else:
                relative_M[r][c] = caw/totals[c]
    return relative_M



def RootFW(m):
    """receives list of lists(adjacency matrix), and adds ROOT node"""
    #new length of matrix
    n = len(m)
    n += 1
    #insert row and column with zer0s
    matrix_root = [[0 for i in range(n)]]
    for r in m:
        temp = [0]
        temp.extend(deepcopy(r))
        matrix_root.append(temp)
    #checking out predators(both top and intermediate)
    for i,item in enumerate(matrix_root[0]):
        for r in range(1,len(matrix_root)):
            matrix_root[0][i] += matrix_root[r][i]
    matrix_root[0][0] = 0
    #reversing previous step, so base becomes connected to ROOT
    for i in range(1, n):
        if matrix_root[0][i]>0:
           matrix_root[0][i] = 0
        else:
            matrix_root[0][i] = 1
    return matrix_root


def DOM_MATRIX(matrix_root):
    """receives rooted matrix and provides dominator matrix"""
    #create matrix full of 1 except first line with 0
    N = len(matrix_root)
    DM = [[1 for j in range(N)] for i in range(N)]
    for i in range(1,N):
        DM[0][i] = 0
    ##
    flag = True
    while flag:
        flag = False
        for i in range(1,N):
            tmp = [1]*N
            for j in range(N):
                if matrix_root[j][i]>0:
                    tmp = [item_one*item_two for item_one, item_two in zip(tmp, DM[j])]#tmp *= DM[j]
            tmp[i] = 1
            if sum(list(map(abs,[item_one-item_two for item_one, item_two in zip(tmp, DM[i])])))!=0:
                flag = True
            DM[i] = tmp
    return DM



def DOM_TREE(DM):
    """from dominator matrix provides dominator tree
"""
    N = len(DM)
    DT = [[0 for j in range(N)] for i in range(N)]

    for j in range(0,N):
        vett = [0 for i in range(0,N)]
        for k in range(0,N):
            if j!=k:
                vett = deepcopy(DM[j])
                vett[k] = 1
                if sum([abs(vett[p] - DM[k][p]) for p in range(0,N)])==0:
                    DT[j][k] = 1
    return DT


def get_TrophicLevels(matrix):#erro cadeias bidirecionais
    """Williams and Martinez (2004) Am Nat, Limits to Trophic Levels and Omnivory in Complex Food Webs: Theory and Data. 
# Algorithm proposed by Levine 1980 J Theor Biol
requires: square matrix
ensures: Trophic levels list"""
    pm = absolute_AM_2_relative_AM(matrix)
    #The algorithm needs consumers on rows
    W = transpose(pm)
    
    WM = np.matrix(W)

    I = np.matrix([[1 if j==i else 0 for j in range(len(pm))] for i in range(len(pm))])
    result = I-WM
    result = result.I
    return list(map(sum,result.tolist()))


def calc_sensitivity(d_m,total_nodes,len_matrix):
    """calculates attack and error sensitivity from dominator matrix
"""
    c_m = [i[1:] for i in d_m[1:]]
    #c_mT = list(map(list, zip(*c_m)))
    c_mT = transpose(c_m)
    assert(len(c_mT)==len(c_m)),'Error'
    #Attack sensitivity ---> largest fraction of taxa disconnected from the "root" (value based on initial food web size)
    #AS = max(map(sum,c_mT))/len(c_m)
    AS = (max(map(sum,c_mT))+(total_nodes - len_matrix))/total_nodes#working!
    
    #Error Sensitivity ---> average number of taxa disconnected from the giant component "root" due to random removal
    ES = sum([i/len(c_m) for i in map(sum,c_mT)])/len(c_m)#working
    return AS,ES


def Allesina_alg(matrix, threshold, PP):
    """remove edges and dependent chains below threshold
"""
    m = np.array([[1 if c>threshold else 0 for c in r] for r in matrix])
    while True:
        no_prey = np.array([sum(c) for c in transpose(m)])
        disconnected = no_prey - PP#-1 ->PP; 0 -> disconected
        if np.count_nonzero(disconnected==0)==0:
            break
        for ri,r in reversed(list(enumerate(disconnected.tolist()))):
            if r==0:
                m = np.delete(m,ri,0)
                m = np.delete(m,ri,1)
                PP = np.delete(PP,ri,0)
    return m


def Allesina(filename,Type):
    """calculates the Attack and Error sensitivity of dominator matrix, by removing edges that are below a threshold.


                        Allesina et al. (2006) approach
"""
    if Type=="AM":
        matrix = absolute_AM_2_relative_AM(open_AM(filename))
    elif Type=="EL":
        nodes_names,edges_w,edges_coords = open_EL(filename)
        matrix = EL_2_matrix(edges_w=edges_w,edges_coords=edges_coords)
    #get primary producers list to never remove
    ##by rooting and keeping this list
    nodes = len(matrix)
    r_m = RootFW(matrix)
    PP = np.array(r_m[0][1:])#primary producers are indexes with 1
    
    for percent in range(0,101):
        threshold = percent/100
        print(threshold)
        m = Allesina_alg(matrix, threshold, PP)
        #construction of the rooted food web
        r_m = RootFW(m)
        #construction of the dominator matrix
        d_m = DOM_MATRIX(r_m)
        #number of nodes in giant component minus the root
        giant_node = len(d_m) - 1
        assert(giant_node==len(m)),'Error'
        #calculating the sensitivities
        A,E = calc_sensitivity(d_m,nodes,len(m))
        
        with open('{}_r_{}.txt'.format(filename,percent),'w') as handle:
            handle.write('{}\t{}\t{}\n'.format(A,E,nodes-giant_node))


def Guerreiro_alg(matrix, threshold, PP):
    """remove edges and dependent chains weight-SUM below threshold 
"""
    #m = np.array([[1 if c>threshold else 0 for c in r] for r in matrix])
    W = transpose(matrix)#list of columns
    m = np.array([[1 if c>0 else 0 for c in r] for r in matrix])
    for ci,c in enumerate(W):
        w = list(enumerate(c))
        shuffle(w)#for randomizing ties
##        w = sorted(w,key=lambda x:x[1],reverse=True)
        w = sorted(w,key=lambda x:x[1],reverse=False)
        w = [i for i in w if i[1]<=threshold]
        outs = 0
        for e in w:
            if outs>=threshold:
                break
            outs += e[1]
            m[e[0],ci] = 0
    del W,w,outs
    
    while True:
        no_prey = np.array([sum(c) for c in transpose(m)])
        disconnected = no_prey - PP#-1 ->PP; 0 -> disconected
        if np.count_nonzero(disconnected==0)==0:
            break
        for ri,r in reversed(list(enumerate(disconnected.tolist()))):
            if r==0:
                m = np.delete(m,ri,0)
                m = np.delete(m,ri,1)
                PP = np.delete(PP,ri,0)
    return m


def GuerreiroScotti(filename, Type, monte_carlo=30):
    #get primary producers list to never remove
    ##by rooting and keeping this list
    if Type=="AM":
        matrix = absolute_AM_2_relative_AM(open_AM(filename))
    elif Type=="EL":
        nodes_names,edges_w,edges_coords = open_EL(filename)
        matrix = EL_2_matrix(edges_w=edges_w,edges_coords=edges_coords)
    nodes = len(matrix)
    r_m = RootFW(matrix)
    PP = np.array(r_m[0][1:])#primary producers are indexes with 1
    
    for percent in range(0,101):
        with open('{}_2_{}.txt'.format(filename,percent),'w') as handle:
            handle.write('')
        threshold = percent/100
        print(threshold)
        for _ in range(monte_carlo):
            m = Guerreiro_alg(matrix, threshold, PP)
            #construction of the rooted food web
            r_m = RootFW(m)
            #construction of the dominator matrix
            d_m = DOM_MATRIX(r_m)
            #number of nodes in giant component minus the root
            giant_node = len(d_m) - 1
            assert(giant_node==len(m)),'Error'
            #calculating the sensitivities
            A,E = calc_sensitivity(d_m,nodes,len(m))
            with open('{}_2_{}.txt'.format(filename,percent),'a') as handle:
                handle.write('{}\t{}\t{}\n'.format(A,E,nodes-giant_node))


def max_prey_only(filename, Type):
    #get primary producers list to never remove
    ##by rooting and keeping this list
    def clean_edges_out(matrix,PP):
        W = transpose(matrix)#list of columns
        m = np.array([[1 if c>0 else 0 for c in r] for r in matrix])
        for ci,c in enumerate(W):
            w = list(enumerate(c))
            w = sorted(w,key=lambda x:x[1],reverse=True)
            for ei,e in enumerate(w):
                if ei==0:
                    m[e[0],ci] = 1.0
                else:
                    m[e[0],ci] = 0.0
        del W,w
##        print(m)
##        while True:
##            no_prey = np.array([sum(c) for c in transpose(m)])
##            disconnected = no_prey - PP#-1 ->PP; 0 -> disconected
##            if np.count_nonzero(disconnected==0)==0:
##                break
##            for ri,r in reversed(list(enumerate(disconnected.tolist()))):
##                if r==0:
##                    m = np.delete(m,ri,0)
##                    m = np.delete(m,ri,1)
##                    PP = np.delete(PP,ri,0)
        return m
    if Type=="AM":
        matrix = absolute_AM_2_relative_AM(open_AM(filename))
    elif Type=="EL":
        nodes_names,edges_w,edges_coords = open_EL(filename)
        matrix = EL_2_matrix(edges_w=edges_w,edges_coords=edges_coords)
    nodes = len(matrix)
    r_m = RootFW(matrix)
    PP = np.array(r_m[0][1:])#primary producers are indexes with 1
    m = clean_edges_out(matrix,PP)

    #construction of the rooted food web
    r_m = RootFW(m)
    #construction of the dominator matrix
    d_m = DOM_MATRIX(r_m)
    #number of nodes in giant component minus the root
    giant_node = len(d_m) - 1
    assert(giant_node==len(m)),'Error'
    #calculating the sensitivities
    A,E = calc_sensitivity(d_m,nodes,len(m))
    with open('{}_strongest_edges.txt'.format(filename),'w') as handle:
        handle.write('{}\t{}\t{}\n'.format(A,E,nodes-giant_node))



if __name__=="__main__":
    print("start")
##    m = np.array([[0,0,1,1,0,0,0],
##         [0,0,0,1,1,0,0],
##         [0,0,0,0,0,1,1],
##         [0,0,0,0,0,0,1],
##         [0,0,0,0,0,0,0],
##         [0,0,0,0,0,0,0],
##         [0,0,0,0,0,0,0]])
##    r_m = RootFW(m)
##    d_m = DOM_MATRIX(r_m)
##    As,Es = calc_sensitivity(m,len(m),len(m))
##    print(As,Es)
    nodes_names,edges_w,edges_coords = open_EL("ESM2---Olmo_Gilabert_2019")
    m = EL_2_matrix(edges_w=edges_w,edges_coords=edges_coords)
    r_m = RootFW(m)
    d_m = DOM_MATRIX(r_m)
    As,Es = calc_sensitivity(d_m,len(d_m),len(d_m))
    print(As,Es)
