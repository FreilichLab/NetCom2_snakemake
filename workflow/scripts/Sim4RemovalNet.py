#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os

netcom_module_relative_path = sys.argv[1]

FULL_PATH = os.path.abspath(netcom_module_relative_path)+"/"

sys.path.append(FULL_PATH)

from netcom import simulation

EC_w_taxaon_dominant_knockout = sys.argv[2]

treatment_seeds = sys.argv[3]

sim_log_dir = sys.argv[4]

treatment = sys.argv[5]

treatment_compounds = sys.argv[6]

treatment_compounds_complementry_table =  sys.argv[7]

treatment_compounds_complementry_list =  sys.argv[8]

f = open(EC_w_taxaon_dominant_knockout, "r")

EC_w_taxaon_dominant_knockout_list = f.readlines()

EC_w_taxaon_dominant_knockout_list = [i.strip("\n").strip().split(" ") for i in EC_w_taxaon_dominant_knockout_list]

f.close()

f = open(treatment_seeds, "r")

treatment_seeds_list = f.readlines()

treatment_seeds_list = [i.strip("\n").strip().split("; ") for i in treatment_seeds_list]

f.close()

Results = []
Results_complementry = []
all_compounds = []
all_compounds_complementry = []
#Simulation taxaon_dominant_knockout_list

for i in EC_w_taxaon_dominant_knockout_list:

    sim_df_T1, steps_df_T1 = simulation(input1 = i.copy(), input2=[treatment_seeds_list[0][1::]], BaseFolder=FULL_PATH)

    sim_res = sim_log_dir + treatment + "_" + i[0] + "_" + "simulation_res.tsv"

    sim_res_steps = sim_log_dir + treatment + "_" + i[0] + "_" + "simulation_steps.tsv"

    sim_df_T1.to_csv(sim_res,sep = "\t")    

    steps_df_T1.to_csv(sim_res_steps,sep = "\t")

    All_compounds_A_sim = sim_df_T1.values.tolist()[0][0]

    Results.append(i[0]+" "+" ".join(All_compounds_A_sim))

    print(i[0])
    if i[0] == 'all':
        all_compounds = All_compounds_A_sim
        set_all_compounds = set(all_compounds)
    else:
        set_All_compounds_A_sim = set(All_compounds_A_sim)
        All_compounds_A_sim_complementry = list(set_all_compounds - set_All_compounds_A_sim)
        all_compounds_complementry.extend(All_compounds_A_sim_complementry)
        Results_complementry.append(i[0]+" "+" ".join(All_compounds_A_sim_complementry))

all_compounds_complementry = list(set(all_compounds_complementry))

#save list of strings as text file

with open(treatment_compounds, 'w') as f:

    for item in Results:

        f.write("%s\n" % item)

f.close()

with open(treatment_compounds_complementry_table, 'w') as f:

    for item in Results_complementry:

        f.write("%s\n" % item)

f.close()

with open(treatment_compounds_complementry_list, 'w') as f:

        f.write(",".join(all_compounds_complementry))

f.close()
