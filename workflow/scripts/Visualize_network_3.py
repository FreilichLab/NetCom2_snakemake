#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 13:05:23 2020

@author: ofir
"""
import matplotlib._color_data as mcd
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import os, sys
#from scipy.spatial import ConvexHull
import matplotlib.patches as mpatches
import random

#For reproduceability
seed = 123
random.seed(seed)
np.random.seed(seed)

"""
Run as:
python3 Network_Figures.py <This is base folder> <drop fragments with size smaller then> \
                          <integer of hubness to filter> <prefix> <treatment_compounds txt file> \
                          <treatment_seeds txt file> <treatment_ECs txt file> \
                          <treatment_compounds_knockout> <EC_w_taxaon_dominant_knockout> <treatment_compounds_pathway>
                          

"""

########################NB################################

#Sort rows treatment_compounds_pathway_ by number of compounds

#Sort rows treatment_compounds_pathway_ by number of compounds
# this is to make sure the the first colured pathways and micobres 
#are colored and refrelct the ones with most compounds
#Currently the order is random
###########################################
#try:
netcom_module_relative_path = sys.argv[1]
min_fragments_size = int(sys.argv[2])
filter_hubness = int(sys.argv[3])

treatment_compounds = sys.argv[4]
f = open(treatment_compounds, "r")
treatment_compounds = f.readline()
treatment_compounds = treatment_compounds.split("; ")
f.close()

print("###################treatment_compounds######")
print(type(treatment_compounds))
print(treatment_compounds)
print("###################treatment_compounds######")

treatment_seeds = sys.argv[5]
f = open(treatment_seeds, "r")
treatment_seeds = f.readline()
#original
treatment_seeds = treatment_seeds.split("; ")
f.close()

print("#############treatment_seeds###########")
print(type(treatment_seeds))
print(treatment_seeds)
print("#############treatment_seeds###########")

treatment_ECs = sys.argv[6] 

#for gene_fuc_ec (from removal)

df_full = pd.read_csv(treatment_ECs, header = 0, sep = '\t')

treatment_ECs = df_full['Unnamed: 0'].to_list()

#treatment_ECs = df_full['EC'].to_list()

#for treatmet_ECs from netcom1
#f = open(treatment_ECs, "r")
#treatment_ECs = f.readline()
#treatment_ECs = treatment_ECs.split("; ")
#f.close()

#for i in range(0,len(treatment_ECs)):
        #fix_ecs to full 4 numbers, if missing add -
#        full_ec = ['-','-','-','-']
#        i_arr = treatment_ECs[i].split(sep = ".")
#        for j in range(0,len(i_arr)):
#            full_ec[j] = i_arr[j]
#        treatment_ECs[i] =  ".".join(full_ec)


print("#########treatment_ECs########")
print(type(treatment_ECs))
print(treatment_ECs)
print("#########treatment_ECs########")

treatment_compounds_knockout_ = sys.argv[7]
f = open(treatment_compounds_knockout_, "r")
treatment_compounds_knockout_ = f.readlines()
treatment_compounds_knockout_ = [i.strip("\n") for i in treatment_compounds_knockout_]
treatment_compounds_knockout_ = [i.split(" ") for i in treatment_compounds_knockout_]
f.close()
########################NB################
#Sort rows treatment_compounds_knockout_ by number of compounds
#######################################
df = pd.DataFrame(treatment_compounds_knockout_)
#count of many compounds in a row and add as a column
df["compounds_count"] = df.count(axis =1)
#sort by compounds count in desecing order
df.sort_values(by = 'compounds_count',ascending=False,inplace = True)
df.drop(["compounds_count"], axis = 1, inplace = True)
#This command coverts back list of list with only non None values
treatment_compounds_knockout_ = [[y for y in x if pd.notna(y)] for x in df.values.tolist()]
print("#########treatment_compounds_knockout_#########")
print(type(treatment_compounds_knockout_))
print(treatment_compounds_knockout_)
print("#########treatment_compounds_knockout_#########")

EC_w_taxaon_dominant_knockout_ = sys.argv[8]
f = open(EC_w_taxaon_dominant_knockout_, "r")
EC_w_taxaon_dominant_knockout_ = f.readlines()
EC_w_taxaon_dominant_knockout_ = [i.strip("\n") for i in EC_w_taxaon_dominant_knockout_]
EC_w_taxaon_dominant_knockout_ = [i.split(" ") for i in EC_w_taxaon_dominant_knockout_]
f.close()

#remove first element which is "all"

EC_w_taxaon_dominant_knockout_.reverse()

EC_w_taxaon_dominant_knockout_.pop()

EC_w_taxaon_dominant_knockout_.reverse()

print("###########EC_w_taxaon_dominant_knockout_########")
print(type(EC_w_taxaon_dominant_knockout_))
print(EC_w_taxaon_dominant_knockout_)
print("###########EC_w_taxaon_dominant_knockout_########")

#for i in range(0,len(EC_w_taxaon_dominant_knockout_)):
        #fix_ecs to full 4 numbers, if missing add -
#        full_ec = ['-','-','-','-']
#        i_arr = EC_w_taxaon_dominant_knockout_[i][1].split(sep = ".")
#        for j in range(0,len(i_arr)):
#            full_ec[j] = i_arr[j]
#        EC_w_taxaon_dominant_knockout_[i][1] =  ".".join(full_ec)

treatment_compounds_pathway_ = sys.argv[9]
f = open(treatment_compounds_pathway_, "r")
treatment_compounds_pathway_ = f.readlines()
f.close()
treatment_compounds_pathway_ = [i.strip("\n").split(";") for i in treatment_compounds_pathway_]

########################NB################
#Sort rows treatment_compounds_pathway_ by number of compounds
#######################################

df = pd.DataFrame(treatment_compounds_pathway_)
#count of many compounds in a row and add as a column
df["compounds_count"] = df.count(axis =1)
#sort by compounds count in desecing order
df.sort_values(by = 'compounds_count',ascending=False,inplace = True)
df.drop(["compounds_count"], axis = 1, inplace = True)
#This command coverts back list of list with only non None values
treatment_compounds_pathway_ = [[y for y in x if pd.notna(y)] for x in df.values.tolist()]

treatment_compounds_network_notFiltered_ = sys.argv[10]
treatment_removal_network_knockout_png_ = sys.argv[11]
treatment_removal_network_knockout_pdf_ = sys.argv[12]
treatment_compounds_count_table = sys.argv[13]
treatment_compounds_table = sys.argv[14]

####compute count and complementry compound tables for (taxa,pathway) Start

count_table = {}

compound_table = {}

for i in range(0, len(treatment_compounds_pathway_)):

    #go over pathways

    pathway = treatment_compounds_pathway_[i][0]

    #remove empty as denoted as '', since '' counts as 1 in count table

    if treatment_compounds_pathway_[i][1] == "":

        treatment_compounds_pathway_[i].pop()

    set_pathways = set(treatment_compounds_pathway_[i])

    count_table[pathway] = {}

    compound_table[pathway] = {}

    for j in range(0, len(treatment_compounds_knockout_)):
        
        taxa = treatment_compounds_knockout_[j][0]

        set_compounds = set(treatment_compounds_knockout_[j])

        count = len(set_pathways.intersection(set(set_compounds)))

        compound_list = list(set_pathways.intersection(set(set_compounds)))

        compound_str = ";".join(compound_list)

        #This is nested dictionries

        count_table[pathway][taxa] =  count

        compound_table[pathway][taxa] = compound_str


df_count_table = pd.DataFrame.from_dict(count_table)

df_compound_table = pd.DataFrame.from_dict(compound_table)

df_count_table.to_csv(treatment_compounds_count_table, sep = "\t")

df_compound_table.to_csv(treatment_compounds_table, sep = "\t")

####compute count and complementry compound tables for (taxa,pathway) End 

for i, vali in enumerate(treatment_compounds_pathway_):
    for j, valj in enumerate(vali):
        treatment_compounds_pathway_[i][j] = valj.strip()



print("##########treatment_compounds_pathway_###########")
print(type(treatment_compounds_pathway_))
print(treatment_compounds_pathway_)
print("##########treatment_compounds_pathway_###########")



#netcom_module_relative_path = "./"
#min_fragments_size = 1
#filter_hubness = 100
#
#
#treatment_compounds = "allCompounds_BjSa.txt"
#f = open(netcom_module_relative_path+treatment_compounds, "r")
#treatment_compounds = f.readline()
#treatment_compounds = treatment_compounds.split(";")
#f.close()
#
#treatment_seeds = "Seeds_BjSa.txt"
#f = open(netcom_module_relative_path+treatment_seeds, "r")
#treatment_seeds = f.readline()
#treatment_seeds = treatment_seeds.split(";")
#f.close()
#
#treatment_ECs = "EC_ALL.txt"
#f = open(netcom_module_relative_path+treatment_ECs, "r")
#treatment_ECs = f.readline()
#treatment_ECs = treatment_ECs.split(";")
#f.close()
#
##treatment_compounds_knockout_ = "Compounds_BjSa_order.txt"
#treatment_compounds_knockout_ = "Compounds_BjSa_order.txt"
#f = open(netcom_module_relative_path+treatment_compounds_knockout_, "r")
#treatment_compounds_knockout_ = f.readlines()
#treatment_compounds_knockout_ = [i.strip("\n") for i in treatment_compounds_knockout_]
#treatment_compounds_knockout_ = [i.split(" ") for i in treatment_compounds_knockout_]
#f.close()
#
##EC_w_taxaon_dominant_knockout_ = "Enzymes_BjSa_order.txt"
#EC_w_taxaon_dominant_knockout_ = "Enzymes_BjSa_Order_1.txt"
#f = open(netcom_module_relative_path+EC_w_taxaon_dominant_knockout_, "r")
#EC_w_taxaon_dominant_knockout_ = f.readlines()
#EC_w_taxaon_dominant_knockout_ = [i.strip("\n") for i in EC_w_taxaon_dominant_knockout_]
#EC_w_taxaon_dominant_knockout_ = [i.split(" ") for i in EC_w_taxaon_dominant_knockout_]
#f.close()
#
#treatment_compounds_pathway_ = "Pathways_BjSa_order.txt"
#f = open(netcom_module_relative_path+treatment_compounds_pathway_, "r")
#treatment_compounds_pathway_ = f.readlines()
#treatment_compounds_pathway_ = [i.strip("\n").split(";") for i in treatment_compounds_pathway_]
#for i, vali in enumerate(treatment_compounds_pathway_):
#    for j, valj in enumerate(vali):
#        treatment_compounds_pathway_[i][j] = valj.strip()
#f.close()
#print(treatment_compounds_pathway_)


#High_alpha_pathways = ['Glucosinolate biosynthesis', 'Linoleic acid metabolism', 'Tyrosine metabolism', 'Geraniol degradation', 'Porphyrin and chlorophyll metabolism']

High_alpha_pathways = ['Nicotinate and nicotinamide metabolism',
'Lysine degradation',
'Vitamin B6 metabolism',
'Pantothenate and CoA biosynthesis',
'Pentose phosphate pathway']


#Load database files
#########################
df_el_ = pd.read_csv(netcom_module_relative_path+"data/DB/full_enzymes_labels_jun.txt", sep="|")
df_ecMapping_ = pd.read_csv(netcom_module_relative_path+"data/DB/ec_reac_mapping_jun.txt", sep="|")
df_reactions_ = pd.read_csv(netcom_module_relative_path+"data/DB/reactions_3_balanced.txt", sep="|")
df_ec_to_compoundIndex_ = pd.read_csv(netcom_module_relative_path+"data/DB/compounds_lables_jun_1.txt", sep="|")

def CreateCompoundsNetwork(All_compounds_B=treatment_compounds,
                           Seeds_B=treatment_seeds,
                           ECs_All=treatment_ECs,
                           df_el = df_el_,
                           df_ecMapping = df_ecMapping_,
                           df_reactions = df_reactions_,
                           df_ec_to_compoundIndex = df_ec_to_compoundIndex_,
                           min_fragments_size = 1,#size of fragments to drop
                           filter_hubness = filter_hubness,#filter by max number of node connections
                           treatment_compounds_network_notFiltered = treatment_compounds_network_notFiltered_,
                           treatment_removal_network_knockout_png = treatment_removal_network_knockout_png_,
                           treatment_removal_network_knockout_pdf = treatment_removal_network_knockout_pdf_,
                           treatment_compounds_knockout = treatment_compounds_knockout_,
                           EC_w_taxaon_dominant_knockout = EC_w_taxaon_dominant_knockout_,
                           treatment_compounds_pathway = treatment_compounds_pathway_):


    #df_el = df_el_
    #df_ecMapping = df_ecMapping_
    #df_reactions = df_reactions_
    #df_ec_to_compoundIndex = df_ec_to_compoundIndex_
    
    
    All_compounds_B = All_compounds_B+Seeds_B    
    
    IndexEnzymeList = df_el[["Index", "Enzyme_Code"]].values.tolist()
    DictEL = {}
    
    for i in IndexEnzymeList:
        DictEL[i[1]] = i[0]
        DictEL[i[0]] = i[1]
        
    #EClist = [DictEL[i] for i in EClist]
    EClist_tmp = []
    Errors = []
    for i in ECs_All:
        try:
            EClist_tmp.append(DictEL[i])
        except:
            Errors.append(i)
            print("Cant find EC "+i)
    #Write ECs couldnt be found
    with open(netcom_module_relative_path+'unfound_ECs.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % x for x in Errors)
    
    df_ecMapping = df_ecMapping[df_ecMapping["Enzyme_index"].isin(EClist_tmp)]
    ListOfRelevantReactions = df_ecMapping["Reactions"].values.tolist()
    
    df_ecMapping["Reactions_list"] = [ i.split(":")[1].strip().split(" ") for i in ListOfRelevantReactions]
    
    #This to be used as Reaction to Enzyme_index 
    df_ecMapping_expolded = df_ecMapping[["Enzyme_index", "Reactions_list"]].explode("Reactions_list")
    
    flat_list = []
            
    for i in ListOfRelevantReactions:
        #print(i)
        for j in i.split(":")[1].strip().split(" "):
            #print(j)
            try:
                flat_list.append(int(j.strip()))
            except:
                print("invalid for int")
    
    #Save enzyme information of each reaction.
    
    df_reactions = df_reactions[df_reactions["Reactions"].isin(flat_list)]
    
    #Fix reaction directions - flip =
    
    tempList = []
    
    for i in df_reactions.values.tolist():
        if i[1] == "=":
            x = i
            x[3], x[2] = x[2], x[3]
            tempList.append(x)
    
    df_reactions = pd.concat([df_reactions, pd.DataFrame(tempList, columns=list(df_reactions.columns))])
    
    #Save the Enzymes linked to reaction. 
    #reactions = df_reactions["Reactions"].values.tolist()
    
    
    l = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Left"].values.tolist()]
    r = [i.split(":")[1].lstrip().rstrip().split(",") for i in df_reactions["Right"].values.tolist()]
    
      
    DictComp = {}
    
    for i in df_ec_to_compoundIndex[["Index", "Compound_Code"]].values.tolist():
        DictComp[i[1]] = i[0]
        DictComp[i[0]] = i[1]
        
      
    df = pd.DataFrame([l, r, df_reactions["Reactions"].values.tolist()]).T           
    df = df.explode(0)         
    df = df.explode(1)
    
    #map reaction to enzyme (col 4)
    ##Might be a problem here, when multiple enzymes are related to a single reaction
    
    ReactionToEnzymeDict = { i[1]:i[0] for i in df_ecMapping_expolded.values.tolist() }
    cols = ["Left", "Right", "Reaction"]
    df.columns = cols
    
    def map_enz(x):
        return ReactionToEnzymeDict[str(x)]
    
    df["Enzyme"] = df["Reaction"].apply(map_enz)
    
    df_grouped = df.groupby(["Left","Right"])["Enzyme"].apply(list)
    df_grouped = df_grouped.to_frame().reset_index()
    
    def enzInd_to_enzEC(l):
        return [DictEL[i] for i in l] 
    
    df_grouped["Enzyme"] = df_grouped["Enzyme"].apply(enzInd_to_enzEC)
    
    
    #Fix Compound IDs
    
    def CompInd_to_CompC(x):
        try:
            return DictComp[int(x)]
        except:
            return np.nan
    
    df_grouped["Left"] = df_grouped["Left"].apply(CompInd_to_CompC)
    df_grouped["Right"] = df_grouped["Right"].apply(CompInd_to_CompC)
    
    #Save only rows contain compounds from the all_compounds list
    #df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_B+All_compounds_A)) &
    #                        (df_grouped["Right"].isin(All_compounds_B+All_compounds_A))]    
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(All_compounds_B)) &
                            (df_grouped["Right"].isin(All_compounds_B))]    
    
  
    #Colors must match to the patches
    #all_patch_colors =  ["brown", "coral", "cyan", "darksalmon", "darkviolet", "gold", "green", "lavenderblush",
    #                     "lime", "magenta", ]
    
    all_patch_colors = ['red','blue','green','yellow', 'violet','peachpuff' ]


    # strong_colors = {
    #         'blue': '#0000FF',
    #         'darkorange': '#FF8C00',
    #         'teal': '#008080',
    #         'limegreen': '#32CD32',
    #         'darkred': '#8B0000',
    #         'deeppink': '#FF1493'
    #                        }

    # strong_colors = list(strong_colors.keys())
    #Set node color
    node_colors_Dict = {}
    
    all_nodes = list(set(df_grouped["Left"].values.tolist()+df_grouped["Right"].values.tolist()))

    for i in all_nodes:
        node_colors_Dict[i] = "gray"

    for i, vali in enumerate(treatment_compounds_knockout):
        for j in vali:
            if j in all_nodes:
                if i >= len(all_patch_colors):
                    node_colors_Dict[j] = all_patch_colors[len(all_patch_colors) - 1]
                else:
                    node_colors_Dict[j] = all_patch_colors[i]
    
   
    

    
    def ColorEnzyme(l):
        c = "gray"
        l = [str(i) for i in l]
        for i, vali in enumerate(EC_w_taxaon_dominant_knockout):
            if [x for x in l if x in vali]:
                if i >= len(all_patch_colors):
                    c = all_patch_colors[len(all_patch_colors) - 1]
                else:
                    c = all_patch_colors[i]
                
        return c
        
    df_grouped["Color"] = df_grouped["Enzyme"].apply(ColorEnzyme)
    
    
    def setEdgeWidthByColor(x):
        if x != "gray":
            return 1.5
        else:
            return 0.75
        
    df_grouped["Edge_width"] = df_grouped["Color"].apply(setEdgeWidthByColor)
    
    df_grouped.to_csv(treatment_compounds_network_notFiltered)
    
    #df_grouped["count"] = df.groupby(["Left","Right"]).count().reset_index()["Reaction"]
    
    
    #Filter df by nodes appearance
    temp_df = df_grouped["Left"].value_counts().reset_index()
    keeplist = temp_df[temp_df["Left"] <= filter_hubness]["index"].values.tolist()
    temp_df = df_grouped["Right"].value_counts().reset_index()
    keeplist += temp_df[temp_df["Right"] <= filter_hubness]["index"].values.tolist()
    keeplist = list(set(keeplist))
    
    df_grouped = df_grouped[(df_grouped["Left"].isin(keeplist)) & (df_grouped["Right"].isin(keeplist))]
    
    
    
    
    
    #Network
    G = nx.Graph()
    
    for i in df_grouped[["Left", "Right", "Color","Edge_width"]].values.tolist():
        G.add_edge(i[0], i[1], color=i[2], width=i[3])
        
    #Filter out isolates
    #G.remove_nodes_from(list(nx.isolates(G)))
    
    #Filter out fragments
    for component in list(nx.connected_components(G)):
        if len(component)<min_fragments_size:
            for node in component:
                G.remove_node(node)
    
    
    
    edges,colors = zip(*nx.get_edge_attributes(G,'color').items())
    edges,widths = zip(*nx.get_edge_attributes(G,'width').items())

    node_color = []
    
    for node in list(G.nodes()):
        node_color.append(node_colors_Dict[node])
          
    nodes_of_largest_component  = max(nx.connected_components(G), key = len)
    largest_component = G.subgraph(nodes_of_largest_component)
    
    pos1 = nx.spring_layout(G,k=0.075,iterations=200, seed=seed)
    
    pos2 = nx.spring_layout(G, pos=pos1,fixed=nodes_of_largest_component,
                           k=0.001,iterations=0,seed=seed)
    #Label colored nodes
    
    flatten_compounds = [item for sublist in treatment_compounds_knockout for item in sublist]
    
    nodes_labeldict = {}
    for i in G.nodes():
        if i in flatten_compounds:
            nodes_labeldict[i] = i
        else: nodes_labeldict[i] = ""
  

    fig = plt.figure(figsize=(50,50))
    nx.draw_networkx_edges(G, pos2, alpha=0.25, width=widths, edge_color=colors, arrows=False)
    nx.draw_networkx_nodes(G, pos2, node_color=node_color, node_size=30,
                                   alpha=0.5)
    nx.draw_networkx_labels(G,pos2,nodes_labeldict,font_size=4,font_color='black')
    
# =============================================================================
    PathwayDict = {}
    for i, vali in enumerate(treatment_compounds_pathway):
        for j in vali[1:]:
            if j in PathwayDict.keys():
                PathwayDict[j].append(vali[0])
            else:        
                PathwayDict[j] = [vali[0]]
    
    
    ax=plt.gca()

    #Create node-position and colors
    patchNodePositions = []
    PatchNodeNames = []
    handles = []
    handles2 = []

    #Legend 1
    for i, vali in enumerate(treatment_compounds_pathway):
        #node to position
        for j in vali:
            if j in list(nodes_of_largest_component):
                # this was the place i've made a change. the code didn't work, so i turned it back.
                # probably this is the block that is responsible for the mix-up in pathways order.
                # worst case - i can do it manually - add each color in the corresponding index to specific pathway.
                patchNodePositions.append(pos2[j])
                PatchNodeNames.append(j)
             
        if len(PatchNodeNames) > 0:
            #if i+len(treatment_compounds_knockout) >= len(all_patch_colors):
            if i >= len(all_patch_colors):
                handles.append(mpatches.Patch(color=all_patch_colors[len(all_patch_colors) -1], label=vali[0]))
                #handles.append(mpatches.Patch(color=all_patch_colors[i+len(treatment_compounds_knockout)], label=vali[0]))
            else:
                #handles.append(mpatches.Patch(color=all_patch_colors[i+len(treatment_compounds_knockout)], label=vali[0]))
                handles.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0]))
        
        for cl in PatchNodeNames:
            # plot circles using the RGBA colors
            if len(set(PathwayDict[cl]).intersection(set(High_alpha_pathways))) > 0:
                #if i+len(treatment_compounds_knockout) >= len(all_patch_colors):
                if i >= len(all_patch_colors):
                    circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[len(all_patch_colors) -1], fill=True, alpha=0.8)
                    #circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(treatment_compounds_knockout)], fill=True, alpha=0.8)
                else:
                    #circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(treatment_compounds_knockout)], fill=True, alpha=0.8)
                    circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i], fill=True, alpha=0.8)
                ax.add_artist(circle)
                print(PathwayDict[cl])

            else:
                #if i+len(treatment_compounds_knockout) >= len(all_patch_colors):
                if i >= len(all_patch_colors):
                    circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[len(all_patch_colors) -1], fill=True, alpha=0.2)
                    #circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(treatment_compounds_knockout)], fill=True, alpha=0.2)
                else:
                    #circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i+len(treatment_compounds_knockout)], fill=True, alpha=0.2)
                    circle = plt.Circle(pos2[cl], 0.01, color=all_patch_colors[i], fill=True, alpha=0.2)

                ax.add_artist(circle)
        
        PatchNodeNames = []



    #legend 2
    for i, vali in enumerate(treatment_compounds_knockout):  
        if i >= len(all_patch_colors):
            handles2.append(mpatches.Patch(color=all_patch_colors[len(all_patch_colors) -1] , label=vali[0]))
            #handles2.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0]))
        else:
            handles2.append(mpatches.Patch(color=all_patch_colors[i], label=vali[0]))

    

    # plot the legend


    leg1 = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc=('upper left'), fontsize=20, title="Pathways")
    
    #lower_leg_anch = leg1.get_window_extent().transformed(ax.transAxes.inverted()).get_points()    
    
    leg1.get_title().set_fontsize('20') #legend 'Title' fontsize
    ax.add_artist(leg1)    
    leg2 = plt.legend(handles=handles2, bbox_to_anchor=(1.05, 0.3), loc=('upper left'), fontsize=20, title="Order")
    leg2.get_title().set_fontsize('20') #legend 'Title' fontsize
    ax.add_artist(leg2)    


    #Title = str("Nodes:"+str(len(node_color))+" Edges:"+str(len(colors))+" Seeds:"+str(Count_seeds)+" Unique:"+str(Count_unique)+" Total colored nodes:"+str(Count_unique+Count_seeds) +" Colored edges:"+str(count_colored_edges))
    Title = str("Nodes:"+str(len(node_color))+" Edges:"+str(len(colors)))
    plt.title(Title, fontdict=None, loc='center')
    plt.axis('equal')    
    fig.savefig(treatment_removal_network_knockout_png, bbox_inches='tight')
    fig.savefig(treatment_removal_network_knockout_pdf, bbox_inches='tight')
    plt.close()
  
CreateCompoundsNetwork(All_compounds_B = treatment_compounds,
                           Seeds_B = treatment_seeds,
                           ECs_All = treatment_ECs,
                           df_el = df_el_,
                           df_ecMapping = df_ecMapping_,
                           df_reactions = df_reactions_,
                           df_ec_to_compoundIndex = df_ec_to_compoundIndex_,
                           min_fragments_size = 1,#size of fragments to drop
                           filter_hubness = filter_hubness,#filter by max number of node connections                          
                           treatment_compounds_network_notFiltered = treatment_compounds_network_notFiltered_,
                           treatment_removal_network_knockout_png = treatment_removal_network_knockout_png_,
                           treatment_removal_network_knockout_pdf = treatment_removal_network_knockout_pdf_,
                           treatment_compounds_knockout = treatment_compounds_knockout_,
                           EC_w_taxaon_dominant_knockout = EC_w_taxaon_dominant_knockout_,
                           treatment_compounds_pathway = treatment_compounds_pathway_)



