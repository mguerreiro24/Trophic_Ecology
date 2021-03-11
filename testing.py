from analysis import *

import os
threshold_step = 1
threshold = 1
threshold_max = 99
monte_carlo = 100


home_path = os.getcwd()
file_folder = os.path.join(home_path,"weighted_FOODWEBS")
folder = "ESM2---Olmo_Gilabert_2019"
filename = os.path.join(os.path.join(file_folder,folder),"data")
Type = 'EL'
g = load_data(filename,Type)
nodes_names,edges_w,edges_coords = process_iGraph(g)
gg = deepcopy(g)
deletions_t = []

wk = gg.es.select(weights_le=.26)
print(len(wk))
wk1 = gg.es.select(weights_le=.9)
print(len(wk1))



