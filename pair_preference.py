# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 12:09:02 2021

@author: Divya
"""

# -*- coding: utf-8 -*-
"""
Created on Sat May  1 11:53:37 2021

@author: HP
"""

import os
import glob
import pandas as pd
from pandas import DataFrame
import copy
import re
import csv

# combine data from multipe files in csv format
def combine_files(file_location):
    filenames = glob.glob(file_location)
    print(len(filenames))
    # concat all files
    combined_csv = pd.concat([pd.read_csv(f) for f in filenames ])
    # print(combined_csv)
    #export to csv
    # combined_csv.to_csv(r'C:\Users\HP\Documents\more_abs\28_abs_data\epitope\pair_prefrence\file_pair\combine_pairs.csv', index=False, encoding='utf-8-sig')
    return combined_csv

# Get one letter code from three letter code
def three_2_one(combined_csv):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    r1 = []
    r2 = []
    df1 = combined_csv
    for idx,row in df1.iterrows():
        for key in d.keys():
            if row['0'] == key:
                r1.append(d[key])
            if row['2'] == key:
                r2.append(d[key])

    df1['r1'] = r1
    df1['r2'] = r2
    df1['pair'] = df1['r1'] + df1['r2']
    # print(df1)
    return df1

# Get number of pairs from df1
def num_pairs(df1):
    mains = {'AA': 0, 'AC': 0, 'AD': 0, 'AE': 0, 'AF': 0, 'AG': 0, 'AH': 0, 'AI': 0, 'AK': 0, 'AL': 0, 'AM': 0, 'AN': 0, 'AP': 0, 'AQ': 0, 'AR': 0, 'AS': 0, 'AT': 0, 'AV': 0, 'AW': 0, 'AY': 0, 'CA': 0, 'CC': 0, 'CD': 0, 'CE': 0, 'CF': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CK': 0, 'CL': 0, 'CM': 0, 'CN': 0, 'CP': 0, 'CQ': 0, 'CR': 0, 'CS': 0, 'CT': 0, 'CV': 0, 'CW': 0, 'CY': 0, 'DA': 0, 'DC': 0, 'DD': 0, 'DE': 0, 'DF': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DK': 0, 'DL': 0, 'DM': 0, 'DN': 0, 'DP': 0, 'DQ': 0, 'DR': 0, 'DS': 0, 'DT': 0, 'DV': 0, 'DW': 0, 'DY': 0, 'EA': 0, 'EC': 0, 'ED': 0, 'EE': 0, 'EF': 0, 'EG': 0, 'EH': 0, 'EI': 0, 'EK': 0, 'EL': 0, 'EM': 0, 'EN': 0, 'EP': 0, 'EQ': 0, 'ER': 0, 'ES': 0, 'ET': 0, 'EV': 0, 'EW': 0, 'EY': 0, 'FA': 0, 'FC': 0, 'FD': 0, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'GA': 0, 'GC': 0, 'GD': 0, 'GE': 0, 'GF': 0, 'GG': 0, 'GH': 0, 'GI': 0, 'GK': 0, 'GL': 0, 'GM': 0, 'GN': 0, 'GP': 0, 'GQ': 0, 'GR': 0, 'GS': 0, 'GT': 0, 'GV': 0, 'GW': 0, 'GY': 0, 'HA': 0, 'HC': 0, 'HD': 0, 'HE': 0, 'HF': 0, 'HG': 0, 'HH': 0, 'HI': 0, 'HK': 0, 'HL': 0, 'HM': 0, 'HN': 0, 'HP': 0, 'HQ': 0, 'HR': 0, 'HS': 0, 'HT': 0, 'HV': 0, 'HW': 0, 'HY': 0, 'IA': 0, 'IC': 0, 'ID': 0, 'IE': 0, 'IF': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IK': 0, 'IL': 0, 'IM': 0, 'IN': 0, 'IP': 0, 'IQ': 0, 'IR': 0, 'IS': 0, 'IT': 0, 'IV': 0, 'IW': 0, 'IY': 0, 'KA': 0, 'KC': 0, 'KD': 0, 'KE': 0, 'KF': 0, 'KG': 0, 'KH': 0, 'KI': 0, 'KK': 0, 'KL': 0, 'KM': 0, 'KN': 0, 'KP': 0, 'KQ': 0, 'KR': 0, 'KS': 0, 'KT': 0, 'KV': 0, 'KW': 0, 'KY': 0, 'LA': 0, 'LC': 0, 'LD': 0, 'LE': 0, 'LF': 0, 'LG': 0, 'LH': 0, 'LI': 0, 'LK': 0, 'LL': 0, 'LM': 0, 'LN': 0, 'LP': 0, 'LQ': 0, 'LR': 0, 'LS': 0, 'LT': 0, 'LV': 0, 'LW': 0, 'LY': 0, 'MA': 0, 'MC': 0, 'MD': 0, 'ME': 0, 'MF': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'MK': 0, 'ML': 0, 'MM': 0, 'MN': 0, 'MP': 0, 'MQ': 0, 'MR': 0, 'MS': 0, 'MT': 0, 'MV': 0, 'MW': 0, 'MY': 0, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0, 'NG': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0, 'NN': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'NY': 0, 'PA': 0, 'PC': 0, 'PD': 0, 'PE': 0, 'PF': 0, 'PG': 0, 'PH': 0, 'PI': 0, 'PK': 0, 'PL': 0, 'PM': 0, 'PN': 0, 'PP': 0, 'PQ': 0, 'PR': 0, 'PS': 0, 'PT': 0, 'PV': 0, 'PW': 0, 'PY': 0, 'QA': 0, 'QC': 0, 'QD': 0, 'QE': 0, 'QF': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QK': 0, 'QL': 0, 'QM': 0, 'QN': 0, 'QP': 0, 'QQ': 0, 'QR': 0, 'QS': 0, 'QT': 0, 'QV': 0, 'QW': 0, 'QY': 0, 'RA': 0, 'RC': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0, 'RT': 0, 'RV': 0, 'RW': 0, 'RY': 0, 'SA': 0, 'SC': 0, 'SD': 0, 'SE': 0, 'SF': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SK': 0, 'SL': 0, 'SM': 0, 'SN': 0, 'SP': 0, 'SQ': 0, 'SR': 0, 'SS': 0, 'ST': 0, 'SV': 0, 'SW': 0, 'SY': 0, 'TA': 0, 'TC': 0, 'TD': 0, 'TE': 0, 'TF': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TK': 0, 'TL': 0, 'TM': 0, 'TN': 0, 'TP': 0, 'TQ': 0, 'TR': 0, 'TS': 0, 'TT': 0, 'TV': 0, 'TW': 0, 'TY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0, 'VF': 0, 'VG': 0, 'VH': 0, 'VI': 0, 'VK': 0, 'VL': 0, 'VM': 0, 'VN': 0, 'VP': 0, 'VQ': 0, 'VR': 0, 'VS': 0, 'VT': 0, 'VV': 0, 'VW': 0, 'VY': 0, 'WA': 0, 'WC': 0, 'WD': 0, 'WE': 0, 'WF': 0, 'WG': 0, 'WH': 0, 'WI': 0, 'WK': 0, 'WL': 0, 'WM': 0, 'WN': 0, 'WP': 0, 'WQ': 0, 'WR': 0, 'WS': 0, 'WT': 0, 'WV': 0, 'WW': 0, 'WY': 0, 'YA': 0, 'YC': 0, 'YD': 0, 'YE': 0, 'YF': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YK': 0, 'YL': 0, 'YM': 0, 'YN': 0, 'YP': 0, 'YQ': 0, 'YR': 0, 'YS': 0, 'YT': 0, 'YV': 0, 'YW': 0, 'YY': 0}
   
    for idx,row in df1.iterrows():
        for key in mains.keys():
            if row['pair'] == key:
                mains[key] += 1
    # print(mains)
    return mains

# total no. of amino acids in spike and antibody protein
def num_aa(df2_ab):
   amino_acids = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0} 
   for i in df2_ab:
       for j in range(len(i)):
           for key in amino_acids:
               if i[j] == key:
                   amino_acids[key] += 1                   
   return amino_acids

# total no. of aa in ab and spike at interface
def num_aa_if(df3_spike):
    amino_acids_if = {"A":0, "R":0, "N":0, "D":0, "C":0, "Q":0, "E":0, "G":0, "H":0, "I":0, "L":0, "K":0, "M":0, "F":0, "P":0, "S":0, "T":0, "W":0, "Y":0, "V":0}
    for i in df3_spike:
        for key in amino_acids_if:
            if i == key:
                amino_acids_if[key] += 1
    return amino_acids_if

#Get pair preference
def pair_pref(mains,aa_ab,aa_spike):
    pair_pref = copy.deepcopy(mains)
    for keys in pair_pref:                                  # get pair preference  
        n1 = 0
        n2 = 0
        if pair_pref[keys] != 0:
            if keys[0] in aa_spike.keys():
                n1 += aa_spike[keys[0]]
            if keys[1] in aa_ab.keys():
                n2 += aa_ab[keys[1]]
            pair_pref[keys] = (pair_pref[keys] / (n1*n2))*100
    return pair_pref

file_location = os.path.join('C:/Users/HP/Documents/more_abs/reanalysis/pair_preference', 'contacts', '*.csv')
file = combine_files(file_location)

# convert 3 letter code to 1 letter
df_up = three_2_one(file)
# df_up.to_csv(r'C:\Users\HP\Documents\more_abs\28_abs_data\epitope\pair_prefrence\file_pair\num_pairs.csv', index=False)

# total no. of pairs
pair_nos = num_pairs(df_up)
dat1 = pd.DataFrame()
dat1 = dat1.append(pair_nos,ignore_index=True)
dat1 = dat1.transpose()
dat1.to_csv('C:/Users/HP/Documents/more_abs/reanalysis/pair_preference/num_pair.csv')
print('succesful')

# total no. of amino acids in spike and antibody protein
df2 = pd.read_excel(r'C:\Users\HP\Documents\more_abs\reanalysis\pair_preference\proteinseq.xlsx')
df2_ab = df2['ab']
df2_spike = df2['spike']
aa_ab = num_aa(df2_ab)
aa_spike = num_aa(df2_spike)
# print(aa_ab)
# print(aa_spike)

# total no. of aa in ab and spike at interface
df3_spike = df_up['r1']
df3_ab = df_up['r2']
aa_ab_if = num_aa_if(df3_ab)
aa_spike_if = num_aa_if(df3_spike)
# print(aa_ab_if)
# print(aa_spike_if)

# Get pair preference for all proteins
pp_prot = pair_pref(pair_nos,aa_ab,aa_spike)
# print(pp_prot)
dataset_prot = pd.DataFrame()
dataset_prot = dataset_prot.append(pp_prot,ignore_index=True)
dataset_prot = dataset_prot.transpose()
# Get pair preference considering total residues at interface
pp_if = pair_pref(pair_nos,aa_ab_if,aa_spike_if)
# print(pp_if)
dataset_if = pd.DataFrame()
dataset_if = dataset_if.append(pp_if,ignore_index=True)
dataset_if = dataset_if.transpose()

# write pair preference to an excel file

writer = pd.ExcelWriter(r'C:\Users\HP\Documents\more_abs\reanalysis\pair_preference\pair_preference_final.xlsx')

# write dataframe to excel sheets
dataset_prot.to_excel(writer, 'pair_protein')
dataset_if.to_excel(writer,'pair_interface')

# dataset_aacount2.to_excel(writer,'aa_count_ab')
# dataset_pref.to_excel(writer,'pair_preference')
# dataset_pref_if.to_excel(writer,'pair_preference_if')
# final.to_excel(writer,'final')
# save the excel file

writer.save()

print('DataFrame is written successfully to Excel Sheet.')
