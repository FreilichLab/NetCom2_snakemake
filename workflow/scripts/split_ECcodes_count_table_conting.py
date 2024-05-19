'''
splits ECcodes_count_table_conting.txt into two count table based on 
treatments.
Input:
1. ECcodes_count_table_conting.txt 
2. a list with treatments 

assumtions:
1. Removal_Pipeline.py (1-4)
2. Step 5 with spcific taxonimic rank (-l)
-l : specifying taxonomic rank according to ncbi taxonomy etc
    species,genus,order,phylum
'''
#exec(open("split_ECcodes_count_table_conting.py").read())

import sys
import pandas as pd

EC_contig_taxa_ctable =  sys.argv[1]

treatment = sys.argv[2]

treatment_EC_contig_taxa_ctable = sys.argv[3]

df_full = pd.read_csv(EC_contig_taxa_ctable, header = 0, sep = '\t')

#select coloumns with treatment as subtring (i.e., treatment columns)

t = df_full.columns[df_full.columns.str.contains(treatment)]

#extract corresponding aboundace data 

df_tretment_counts = df_full[t]

#extract EC and taxa data

df_taxa = df_full[['EC','taxa']]

#add corresponding aboundace data to data

df_treatment = pd.concat([df_taxa,df_tretment_counts],axis = 1)

df_treatment.to_csv(treatment_EC_contig_taxa_ctable, sep = '\t', index = False)
