#######################test_remove.py####################
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 08:49:55 2022

@author: ofir
"""


import sys, os
import pandas as pd
import re
import numpy as np
from multiprocessing import Pool
from ete3 import NCBITaxa
import skbio
from skbio.diversity.alpha import simpson



def create_temp_file(input_, tempfile):
    with open(input_) as f:
        with open(tempfile, 'w') as t:
            data = f.read()
            data = re.sub(r'([0-9][0-9]*_[0-9][0-9]*)_([0-9][0-9]*\t)', r'\1\t\2', data)
            data = re.sub(' ', '_', data)
            t.write(data)

def most_common(lst):
    current = max(set(lst), key=lst.count)
    maxnum = lst.count(current)
    listmax = []
    for el in set(lst):
        if lst.count(el) == maxnum:
            listmax.append(el)
    if len(listmax) == 1 :
        return listmax[0]
    else:
        return 'even'

def manipulate(cont):
    results_list = []
    #global number_current_contig
    #number_current_contig += 1
    #print(str(number_current_contig) + " / " + str(contnum) )
    #print(("Starting to work on contig {} of {}").format(str(), str(contnum)))
    contdata = data.where(data['contig'] == cont).dropna()
    smalldata = contdata['annotation'].str.split(';', expand=True)
    smalldata.replace('', np.nan, inplace=True)
    smalldata.replace('Not_assigned', np.nan, inplace=True)
    most_freq_tax = ''
    common_tax = ''
    levelnum = len(smalldata.columns)
     #I set a parameter (diff) that switch to 1 from 0 if the intra-contig taxonomy differs, so the program stop adding things to the common taxonomy column.
    diff = 0
    for num in range(levelnum):
        anlist = smalldata[num].dropna().tolist()
        if len(set(anlist)) == 1:                        
            most_freq_tax += (';' + anlist[0])
            if diff == 0:
                common_tax += (';' + anlist[0])
        elif len(set(anlist)) > 1:
            toadd = most_common(anlist)
            if toadd != 'even':
                most_freq_tax += (';' + toadd)
                smalldata = smalldata[smalldata[num] == toadd]
                diff = 1
            else:
                break
    if most_freq_tax == common_tax:
        #o.write(cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'Yes' + '\n')
        results_list.append([cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'Yes' + '\n'])
    else:
        #o.write(cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'No' + '\n')
        results_list.append([cont + '\t' + most_freq_tax + '\t' + common_tax + '\t' + 'No' + '\n'])
    
    return results_list

def get_last_taxonomy(string):
    x = string.split(";")
    return x[-1]

def get_tax_rank(string):
    x = string.replace('_', ' ').split(", ")
    for i in range(len(x)):
        try:
            taxid = ncbi.get_name_translator([x[i]])[x[i]]
            rank = ncbi.get_rank(taxid)[taxid[0]]
            if rank != 'no rank' :
                break
        except KeyError:
            continue
    try:
        rank
    except NameError:
        rank = "Unknown_rank"
    return rank.replace(' ', '_')

def get_level_taxonomy(string, level):
    x = string.split(";")
    for i in x:
        if get_tax_rank(i) == level :
            result = i
            break
        if level != 'order':
            if get_tax_rank(i) == 'order':
                result = 'Unidentified_' + str(i)
    try:
        result
    except NameError:
        result = 'Unidentified'
    return result

gene_taxonomy_annotation = sys.argv[1]
gene_functional_annotation = sys.argv[2]
contig_count_table = sys.argv[3]
level = sys.argv[4]

#gene_taxonomy_annotation = 'GeneTaxonomicAnnotation.txt'
#gene_functional_annotation = 'GeneECAnnotation.txt'
#contig_count_table =  'count_table_contigs.txt'
#BaseFolder = './'
#gene_level = 'genus'


#########Script #1 START######################

file = gene_taxonomy_annotation

tempfile = sys.argv[5] #temp_file

get_contig_taxonomy_ranked = sys.argv[6]

output = get_contig_taxonomy_ranked

#Check input file
try:
   with open(file) as f:
       print('I will now start to work. Feel free to do something relaxing while you wait.')
except FileNotFoundError:
    print('It seems that the input file you specified does not exist. I am very sorry, but I cannot proceed :(.')
    sys.exit()

#At first I create the temporary file with the data formatted as I need
create_temp_file(file, tempfile)

data = pd.read_csv(tempfile, names=['contig', 'protein', 'annotation'], sep = '\t')
contlist = set(data['contig'].tolist())
contnum = len(contlist)

results_list_final = Pool().map(manipulate, list(contlist))

results_list_final_flatten = [ent for sublist in results_list_final for ent in sublist]
results_list_final_flatten = [['Contig' + '\t' + 'most_frequent_taxonomy' + '\t' + 'common_taxonomy' + '\t' + 'Same_taxonomy?' + '\n']]+ results_list_final_flatten

o = open(output, "w")
for element in results_list_final_flatten:    
    o.write(element[0])
o.close()



ncbi = NCBITaxa()

data = pd.read_csv(output, sep='\t')

data['taxonomy_rank'] = data['most_frequent_taxonomy'].apply(lambda x: get_tax_rank(get_last_taxonomy(x)))

data.to_csv(path_or_buf = output, sep = '\t', index = False)   

print('Everything done. The final output is a file named \'' + str(output) + '\' with five columns.')
print('1: contig names')
print('2: most frequent taxonomy in the contig')
print('3: taxonomy common to all proteins on the contig')
print('4: the proteins on this contig share the same taxonomy? If you find a \"Yes\" in this column, then column 2 and 3 will be the same.')
print('5: the taxonomy level of the most frequent taxonomy on the contig. Es: if the most frequent taxonomy is "Penicillium", in this column you will find "genus"')

os.remove(tempfile)


#########Script 1 END######################

#########Script 2 START######################
"""
Step 2:
python prepare_edgeR_input_taxonomy_v2.py -o edgeR_input_taxonomy_1.contig.order.txt -c NTC_G210.count_tab.matrix -t get_contig_splitted_2_1.txt -l order
(-l you can choose a specific rank which you're interested in )
"""

count_table = contig_count_table

taxonomy_file = get_contig_taxonomy_ranked

contig_taxa_ctable = sys.argv[7] #conting_taxa_ctable

output = contig_taxa_ctable

# Prepare dataframes and dictionary
#Prepare taxonomy dictionary
ncbi = NCBITaxa()

data = pd.read_csv(taxonomy_file, sep='\t')
data.set_index(data[data.columns[0]], inplace = True)


#contigs_list = data['Contig'].values.tolist()
data["get_level"] = data.apply(lambda x: get_level_taxonomy(x["most_frequent_taxonomy"], level), axis = 1)

dic = {i[0] : i[1] for i in data[["Contig", "get_level"]].values.tolist()}

# Prepare count table
ctable = pd.read_csv(count_table, sep = '\t')
ctable.set_index(ctable['ID'], inplace = True)
ctable.drop(['ID'], axis = 1, inplace = True)
try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')
 
    
print("Proceed to analysis")
    
#prepare empy output dataframe
outdf = pd.DataFrame(index = list(set(list(dic.values()))), columns = ctable.columns, data = 0.0)

# Fill output dataframe
for cont in dic.keys():
    for col in outdf.columns:
        outdf[col][dic[cont]] += ctable[col][cont]
  
outdf.to_csv(path_or_buf = output, sep = '\t')

#########Script 2 END######################

#########Script 3 START######################
#Select only genes with ECs annotation
#Here I assume the gene funtional input file has only ECs annotation

#gene_functional_annotation = sys.argv[2] #gene functional

gene_functional_EC = sys.argv[8] #gene_functional_EC

cmd = 'cp' + ' ' + gene_functional_annotation + ' ' + gene_functional_EC

os.system(cmd)

#fix ECs - converted perl script
file1 = open(gene_functional_EC, 'r')
file1_Lines = file1.readlines()
file1.close()


#Wed 06 Mar 2024 09:59:54 IST
# Gon addition to include to better parse EC numbers
#recived smae output, i.e. inputUltimate_fixed  as original.
# ; spit gives more than one EC per conting, hope this works
#if not choose one, to make it work
fixed_lines = []
for line in file1_Lines:
    if 'EC:' in line:
        line_arr = line.strip('\n').split("\t")
        EC_arr = line_arr[1].split(';')
        for ec in EC_arr:
            ec_num = ec.split('EC:')[1]
            ec_num = ec_num.replace(']"','')
            ec_str = line_arr[0] + '\t'+ ec_num
            fixed_lines.append(ec_str)

inputUltimate_fixed = fixed_lines

#########Script 3 END######################

#########Script 4 START######################

#dictionary is fix ECs from script 3
#That is fixed gene_functional_EC.
dictionary = inputUltimate_fixed #perl output

gene_func_ec_contig_ctable = sys.argv[9] #gene_functional_EC_contig_ctable

output = gene_func_ec_contig_ctable

# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
#with open(dictionary) as t:
for line in dictionary:
    line = line.rstrip().split('\t')
    line[0] = re.sub('_\d+$','',line[0])
    for el in line[1].rstrip().split(' '):
        if line[0] in dic:
            dic[line[0]].append(el)
        else:
            dic[line[0]] = [el]
for id in list(dic.keys()):
    dic[id] = list(set(dic[id]))

# Prepare count table
ctable = pd.read_csv(count_table, sep = '\t')
ctable.set_index(ctable['ID'], inplace = True)
ctable.drop(['ID'], axis = 1, inplace = True)
try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')

print("Proceed to genewise analysis")
newindex = []
for el in dic.values():
    for code in el:
        newindex.append(code)
       
newindex = list(set(newindex))
outdf = pd.DataFrame(index = newindex, columns = ctable.columns, data = 0.0)
for gene in dic.keys():
    for col in outdf.columns:
        try:
            for code in dic[gene]:
                outdf[col][code] += ctable[col][gene]
        except KeyError:
            continue

outdf.to_csv(path_or_buf = output, sep = '\t')

#########Script 4 END######################

#########Script 5 START######################

EC_contig_taxa_ctable = sys.argv[10] #EC_contig_taxa_ctable

output= EC_contig_taxa_ctable

proteins = gene_taxonomy_annotation 

# Prepare dataframes and dictionary
#Prepare dictionary from protein to EC code

dic = {}
for line in dictionary:
    line = line.rstrip().split('\t')
    for el in line[1].rstrip().split(' '):
        if el in dic:
            dic[el].append(line[0])
        else:
            dic[el] = [line[0]]
for id in list(dic.keys()):
    dic[id] = list(set(dic[id]))

# Prepare count table

try:
    ctable.drop(['Len'], axis = 1, inplace = True)
except KeyError:
    print('')

print("Prepare temporary files")      

## prepare protein data
with open(proteins) as f:
    with open(tempfile, 'w') as t:
        data = f.read()
        data = re.sub(r'([0-9][0-9]*_[0-9][0-9]*)_([0-9][0-9]*\t)', r'\1\t\2', data)
        data = re.sub(' ', '_', data)
        t.write(data)

annotation = pd.read_csv(tempfile, names=['contig', 'protein', 'annotation'], sep = '\t')
annotation['protein'] = annotation['protein'].astype(str)
annotation['protein'] = annotation['contig'] + '_' + annotation['protein']
   


#Prepare taxonomy of desired level
ncbi = NCBITaxa()

#read taxonomy file and prepare taxonomy dictionary of desired level
taxa = pd.read_csv(taxonomy_file, sep='\t')
taxa.set_index(taxa['Contig'], inplace = True)

dictax = {}
for cont in taxa['Contig'].tolist():
    dictax[cont] = get_level_taxonomy(taxa['most_frequent_taxonomy'][cont], level)
    
print("Proceed to contigwise analysis")

codes = []
for code in dic.keys():
        codes.append(code)
           
outdf = pd.DataFrame(index = [], columns = ctable.columns, data = 0.0)
samples = ctable.columns
outdf['EC'] = ''
outdf['taxa'] = ''

# Fill output dataframe
for code in codes:
    templist = dic[code]
    tempdf = annotation.where(annotation['protein'].isin(templist)).dropna()
    contlist = tempdf['contig'].drop_duplicates().tolist()
    for cont in contlist:
        try:
            contdf = tempdf.where(tempdf['contig'] == cont).dropna()
            protnum = len(contdf)
            tax = dictax[cont]
            subdf = outdf.where((outdf['EC'] == code) & (outdf['taxa'] == tax)).dropna()
            if len(subdf) > 1:
                print('Something unexpected happened while filling the database. Exiting now.')
                sys.exit()
            elif len(subdf) == 1:
                ind = subdf.index[0]
                for sample in samples:
                    outdf[sample][ind] += round((protnum * ctable[sample][cont]), 2)
            elif len(subdf) == 0 :
                ind = code + '___' + tax
                outdf = outdf.append(pd.DataFrame(index = [ind], columns = outdf.columns, data = 0.0))
                for sample in samples:
                    outdf[sample][ind] += round((protnum * ctable[sample][cont]), 2)
                outdf['EC'][ind] = code
                outdf['taxa'][ind] = tax
        except:
           print("error in parsing raw")

os.remove(tempfile)

cols = outdf.columns.tolist()
cols = cols[-2:] + cols[:-2]
outdf = outdf[cols]

outdf.to_csv(path_or_buf = output, sep = '\t', index = False)

#########Script 5 END######################
