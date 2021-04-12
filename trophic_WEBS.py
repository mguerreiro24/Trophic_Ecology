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


def calc_sensitivity(d_m):
    """calculates attack and error sensitivity from dominator matrix
"""
    c_m = [i[1:] for i in d_m[1:]]
    #c_mT = list(map(list, zip(*c_m)))
    c_mT = transpose(c_m)
    AS = max(map(sum,c_mT))/len(c_m)
    #ES = sum([i/len(c_m) for i in map(sum,c_mT)])/len(c_mT)
    ES = sum([i/len(c_m) for i in map(sum,c_mT)])/(len(c_mT)-1)#dominator matrix
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


def Allesina(matrix):
    """calculates the Attack and Error sensitivity of dominator matrix, by removing edges that are below a threshold.


                        Allesina et al. (2006) approach
"""
    #get primary producers list to never remove
    ##by rooting and keeping this list
    r_m = RootFW(matrix)
    PP = np.array(r_m[0][1:])#primary producers are indexes with 1
    nodes = len(matrix)
    AS = []
    ES = []
    for percent in range(0,101):
        threshold = percent/100
        m = Allesina_alg(matrix, threshold, PP)
        #construction of the rooted food web
        r_m = RootFW(m)
        #construction of the dominator matrix
        d_m = DOM_MATRIX(r_m)
        #number of noves in giant component minus the root
        giant_node = len(d_m) - 1
        #calculating the sensitivities
        A,E = calc_sensitivity(d_m)
        AS.append(A)
        ES.append(E)
    return AS,ES


def Guerreiro_alg(matrix, threshold, PP):
    """remove edges and dependent chains weight-SUM below threshold 
"""
    #m = np.array([[1 if c>threshold else 0 for c in r] for r in matrix])
    W = transpose(matrix)#list of columns
    m = np.array([[1 if c>0 else 0 for c in r] for r in matrix])
    for ci,c in enumerate(W):
        w = list(enumerate(c))
        shuffle(w)
        w = sorted(w,key=lambda x:x[1],reverse=True)
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

def GuerreiroScotti(matrix, threshold, PP, monte_carlo):
    print("man at work")
    #get primary producers list to never remove
    ##by rooting and keeping this list
    r_m = RootFW(matrix)
    PP = np.array(r_m[0][1:])#primary producers are indexes with 1
    nodes = len(matrix)
    AS = []
    ES = []
    for percent in range(0,threshold):
        threshold = percent/100
        for i in range(monte_carlo):
            m = Guerreiro_alg(matrix, threshold, PP)
            #construction of the rooted food web
            r_m = RootFW(m)
            #construction of the dominator matrix
            d_m = DOM_MATRIX(r_m)
            #number of noves in giant component minus the root
            giant_node = len(d_m) - 1
            #calculating the sensitivities
            A,E = calc_sensitivity(d_m)
            AS.append(A)
            ES.append(E)



if __name__=="__main__":
    print("start")
    m = np.array([[0,0,1,1,0,0,0],
         [0,0,0,1,1,0,0],
         [0,0,0,0,0,1,1],
         [0,0,0,0,0,0,1],
         [0,0,0,0,0,0,0],
         [0,0,0,0,0,0,0],
         [0,0,0,0,0,0,0]])
    r_m = RootFW(m)
    d_m = DOM_MATRIX(r_m)
    As,Es = Allesina(m)
    
