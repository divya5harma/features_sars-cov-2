# -*- coding: utf-8 -*-
"""
Created on Thu May 20 21:16:52 2021

@author: Divya
"""
# code for finding epitope identity from residue numbers
import pandas as pd
import itertools
import csv
df = pd.read_excel('table_for_pdbs.xlsx',sheet_name='final')
# match  = 0
# print(len(df))
x = []
for i in range(len(df)):
    x.append([])
    match = 0
    epi = df['Epitope Region'][i].split(';')
    epi = "".join(epi)
    epi = epi.split(' ')
    # print(epi)
    for j in range(len(df)):
        epi2 = df['Epitope Region'][j].split(';')
        epi2 = "".join(epi2)
        epi2 = epi2.split(' ')
        # print('..............')
        # print(epi2)
        # comparisons = [a == b for (a, b) in itertools.product(epi, epi2)]
        comparisons = [match+1 for (a, b) in itertools.product(epi, epi2) if a==b]
        comparisons = sum(comparisons)
        ident = comparisons/min(len(epi),len(epi2))*100
        ident = "{:.2f}".format(ident)
        x[i].append(ident)
        # print(match)
    
    # print('outer loop end')
print(x)  
file = open('max_epitope_identity.csv', 'w+', newline ='') 
# writing the data into the file 
with file:
  write = csv.writer(file) 

  write.writerows(x)        
    # print(epi)
# for j in range(len(df)):
#     epi2 = df['Epitope Region'][j].split(';')
#     epi2 = "".join(epi2)
#     epi2 = epi2.split(' ')
#     print(epi2)
    
