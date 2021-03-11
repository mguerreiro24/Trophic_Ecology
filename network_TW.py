"""https://igraph.org/python/doc/tutorial/tutorial.html"""

from igraph import *
from copy import deepcopy


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


def process_AM(edges_fl):
    edges_w = []
    edges_coords = []
    for r in range(0,len(edges_fl)):
        for c in range(0,len(edges_fl)):
            if float(edges_fl[r][c])>0:
                edges_w.append(float(edges_fl[r][c]))
                edges_coords.append((r,c))
    return edges_w,edges_coords


def process_iGraph(g):
    nodes_names = g.vs["name"]
    edges_w = g.es["weights"]
    edges_coords = g.get_edgelist()
    return nodes_names,edges_w,edges_coords
#-----------------------------------------------------------------------------------------#
##def load_EL(filename, sep="\t"):returns input of process_EL(edges_fl)
##    
def load_AM(filename,h=True,sep="\t"):
    """open adjacency matrixs with header
"""
    if h:
        with open('{}.txt'.format(filename), "r") as handle:
            edges_fl = [list(map(float,i[:-1].split(sep))) for i in handle.readlines()[1:]]
        with open('{}.txt'.format(filename), "r") as handle:
            nodes_names = [i for i in handle.readline()[:-1].split(sep)]
    else:
        with open('{}.txt'.format(filename), "r") as handle:
            edges_fl = [list(map(float,i[:-1].split(sep))) for i in handle.readlines()]
        nodes_names = [i for i in range(len(edges_fl))]
    return nodes_names,edges_fl


#-----------------------------------------------------------------------------------------#
def open_EL(filename, sep="\t"):
    """get edge list out of format and get values for graph
"""
    with open('{}.txt'.format(filename), "r") as handle:
        edges_fl = [i[:-1].split(sep) for i in handle.readlines()]
    nodes_names,edges_w,edges_coords = process_EL(edges_fl)
    return nodes_names,edges_w,edges_coords


def open_AM_H(filename,h=True,sep="\t"):
    """open adjacency matrixs with header and get values for graph
"""
    if h:
        with open('{}.txt'.format(filename), "r") as handle:
            edges_fl = [i[:-1].split(sep) for i in handle.readlines()[1:]]
        with open('{}.txt'.format(filename), "r") as handle:
            nodes_names = [i for i in handle.readline()[:-1].split(sep)]
    else:
        with open('{}.txt'.format(filename), "r") as handle:
            edges_fl = [i[:-1].split(sep) for i in handle.readlines()]
        nodes_names = [i for i in range(len(edges_fl))]
    edges_w,edges_coords = process_AM(edges_fl)
    return nodes_names,edges_w,edges_coords


def open_Matrix(filename,h=True,sep="\t"):
    if not h:
        with open(filename,"r") as handle:
            h = [[float(j) for j in i[:-1].split(sep)] for i in handle.readlines()]
    else:
        with open(filename,"r") as handle:
            h = [[float(j) for j in i[:-1].split(sep)] for i in handle.readlines()[1:]]
    return h


def from_out_O_fmt_files(nodes_names,edges_w,edges_coords):
    """receives objects processed from out of format files 
and outputs a directed graph """
    g = Graph(directed=True)
    g.add_vertices(len(nodes_names))
    g.add_edges(edges_coords)
    g.vs["name"] = nodes_names
    g.es["weights"] = edges_w
    return g


#----------------------------------------------------------------------------------------
def RootingGraph(g):#<--------
    """roots iGraph
"""
    nodes_names = ["root"]
    nodes_names.extend(g.vs["name"])

    edges_coords = [(i[0],i[1]) for i in g.get_edgelist()]
    #getting base edges
    IN = [i[0] for i in edges_coords]
    OUT = [i[1] for i in edges_coords]
    unrooted_connected = set([i for i in IN if i not in OUT])#not including isolated!
    #with isolated nodes now...?
    unrooted = [i for i,_ in enumerate(g.vs["name"]) if i in unrooted_connected or (i not in IN and i not in OUT)]
    new_edges = [(0,i) for i in unrooted]
    new_edges.extend([(i[0]+1,i[1]+1) for i in edges_coords])

    weig = [1.0 for i in unrooted]
    weig.extend(g.es["weights"])
    return from_out_O_fmt_files(nodes_names,weig,new_edges)


def RootFW(m,nodes_names,edges_coords,edges_w):
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
    r_coords = []#plus, retrieving new edges from root to bases
    r_w = []#plus, retrieving new edges weights from root to bases
    for i in range(1, n):
        if matrix_root[0][i]>0:
           matrix_root[0][i] = 0
        else:
            matrix_root[0][i] = 1
            r_coords.append((0,i))
            r_w.append(1)
    #processing labels and edges(opposi)
    r_nodes_names = ["Root"]
    r_nodes_names.extend(nodes_names)
    for r in edges_coords:
        r_coords.append((r[0]+1,r[1]+1))
    for r in edges_w:
        r_w.append(r)
    return matrix_root,r_nodes_names,r_coords,r_w


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
                #if matrix_root[j][i]==1:
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


def Dom_tree(matrix_root):
    """receives rooted matrix and provides dominator tree"""
    DM = DOM_MATRIX(matrix_root)

    return DOM_TREE(DM)#DM here is the Dominator matrix


def get_TrophicLevels(matrix):#erro cadeias bidirecionais
    """Williams and Martinez (2004) Am Nat, Limits to Trophic Levels and Omnivory in Complex Food Webs: Theory and Data. 
# Algorithm proposed by Levine 1980 J Theor Biol
requires: square matrix
ensures: Trophic levels list"""
    pm = absolute_AM_2_relative_AM(matrix)
    #The algorithm needs consumers on rows
    W = list(map(list, zip(*pm)))#transpose

    from numpy import matrix
    WM = matrix(W)

    I = matrix([[1 if j==i else 0 for j in range(len(pm))] for i in range(len(pm))])
    result = I-WM
    result = result.I
    return list(map(sum,result.tolist()))

if __name__=="__main__":
    print("hello")
    filename="dummy5"
    absolute_AM = open_Matrix(filename+'.txt',h=False,sep=" ")






















