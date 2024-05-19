import pandas as pd
import sys

#exec(open("vis_pathway.py").read())

#pathway_table_file =  ./results/netcom/cf_Compounds_pathway.csv
pathway_table_file = sys.argv[1]

treatment_compounds_complementry_list = sys.argv[2]

treatment_pathway_table_file = sys.argv[3]

df_full = pd.read_csv(pathway_table_file, header = 0, sep = ',')

file1 = open(treatment_compounds_complementry_list, 'r')

file1_Lines = file1.readlines()

file1.close()

compounds_complementry_arr = file1_Lines[0].split(sep = ",")

compounds_complementry_set = set(compounds_complementry_arr)

df_full['entities'] = df_full['entities'].apply(lambda x: x.replace("[","").replace("]","").replace("'","").replace(" ","").split(sep = ","))

df_full['entities'] = df_full['entities'].apply(lambda x: list(compounds_complementry_set.intersection(set(x))))

df_full['entities'] = df_full['entities'].apply(lambda x: ";".join(x))

df_full['Pathway_vis'] = df_full['Pathway'] + ";" + df_full['entities']

df_full['Pathway_vis'].to_csv(treatment_pathway_table_file, sep = '\t', index = False,header = False)

