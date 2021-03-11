#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#                                                                     #
#                     Miguel Guerreiro                                #
#                          2021                                       #
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
import os
from time import time
#taxonL = ('class','class','class','class','class','class','class','class','phylum','class','family','order')
#taxonN = ('Cephalopoda','Actinopterygii','Cephalaspidomorphi','Myxini','Chondrichthyes','Malacostraca','Ostracoda','Maxillopoda','Ctenophora','Scyphozoa','Porpitidae','Siphonophorae')

def buffering_readlines(filename,sep="\n",buffersize=102400):
    """buffering data from file
requires: filename[str()]    filepath to file to buffer
          sep[str()]         separator of entries
          buffersize[int()]  amount of characters to load per each iteration of buffereing process
ensures: generator object
"""
    with open(filename,"r", encoding="utf8") as input_file:
        buffer_data = input_file.read(buffersize)
        while len(buffer_data)>0:
            Muninn = buffer_data
            processing = [i for i in buffer_data.split(sep)[:-1]]
            for i in processing:
                yield i
            buffer_data = buffer_data.split(sep)[-1] + input_file.read(buffersize)
            if Muninn==buffer_data:
                yield buffer_data
                break


def collect_data(filename, target_i, source_i, interaction_i, taxon_target, degrees="all"):
    """
"""
    if degrees.lower()=="all":
        I,O = True,True
    elif degrees.lower()=="in":
        I,O = True,False
    elif degrees.lower()=="out":
        I,O = False,True
    out = []
    gen = buffering_readlines(filename)
    for ii,i in enumerate(gen):
        if ii==0:
            continue
        data = [j for j in i.split("\t")]
        if data[interaction_i]=="get eaten by" or data[interaction_i]=="eats":
            if taxon_target in data[target_i] or taxon_target in data[source_i]:
                flag = False
                if data[interaction_i]=="get eaten by":
                    if (I and taxon_target in data[target_i]) or (O and taxon_target in data[source_i]):
                        prey = {k:v for v,k in zip(data[source_i].split(" | "),data[source_i+2].split(" | "))}
                        preyTaxonIds = data[1].split(" | ")
                        predator = {k:v for v,k in zip(data[target_i].split(" | "),data[target_i+2].split(" | "))}
                        predatorTaxonIds = data[37].split(" | ")
                        flag = True
                elif data[interaction_i]=="eats":
                    if (I and taxon_target in data[source_i]) or (O and taxon_target in data[target_i]):
                        prey = {k:v for v,k in zip(data[target_i].split(" | "),data[target_i+2].split(" | "))}
                        preyTaxonIds = data[37].split(" | ")
                        predator = {k:v for v,k in zip(data[source_i].split(" | "),data[source_i+2].split(" | "))}
                        predatorTaxonIds = data[1].split(" | ")
                        flag = True
                if flag:
                    out.append({"prey":prey,
                                "preyTaxonIds":preyTaxonIds,
                                "predator":predator,
                                "predatorTaxonIds":predatorTaxonIds,
                                "sourceLastSeenAtUnixEpoch":data[-1],
                                "sourceDOI":data[-2],
                                "sourceArchiveURI":data[-3],
                                "sourceCitation":data[-5],
                                "referenceUrl":data[-6],
                                "referenceDoi":data[-7],
                                "referenceCitation":data[-8],
                                "eventDateUnixEpoch":data[-10],
                                "decimalLongitude":data[-13],
                                "decimalLatitude":data[-14]})

    return out


def main(filename, target={'Cephalopoda': 'class'}):
    """
"""
    #separated by " | "
    tg = 40#index of target taxonomical path
    so = 4#index of source taxonomical path 
    #+2 level of taxonomy of each of the taxon names in the path
    it = 34#index of interaction type | "get eaten by" or "eats"
    home_path = os.getcwd()
    file_folder = os.path.join(home_path,"GLOBI")
    filename = os.path.join(file_folder,"interactions.tsv")
    c_time = time()
    data = collect_data(filename, tg, so, it, list(target.keys())[0])
    print("cephs {} s".format(time()-c_time))
##    predators = []
##    for i in data:
##        try:
##            if i["predator"]["class"]!=list(target.keys())[0]:
##                predators.append(i["predator"]["species"])
##        except KeyError:
##            pass
##    predators = set(predators)
##    dataa = []
##    print("preds {}".format(len(predators)))
##    for i in predators:
##        c_time = time()
##        print(i)
##        dataa.extend(collect_data(filename, tg, so, it, i, "IN"))
##        print(time()-c_time)
    return data#,dataa


if __name__=="__main__":
    home_path = os.getcwd()
    file_folder = os.path.join(home_path,"GLOBI")
    filename = os.path.join(file_folder,"interactions.tsv")
    data= main(filename)#,dataa = main(filename)
    with open("interaction_cephs.csv","w") as handle:
        handle.write("\n".join(map(",".join,[[" | ".join(i["prey"].values())," | ".join(i["predator"].values())] for i in data]))+"\n")
    
        
