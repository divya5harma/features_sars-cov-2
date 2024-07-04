# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 15:34:17 2021

@author: Divya
"""

import os
import glob
import pandas as pd
from pandas import DataFrame
import copy
import re
import csv
import numpy as np
from statistics import mean

filenames = glob.glob(os.path.join('C:/Users/HP/Documents/more_abs/reanalysis/num_contacts', 'reswise_contacts_4a', '*.csv'))
total = []
for f in filenames:
    df = pd.read_csv(f)
    tot = len(df)
    total.append(tot)
    
tot_e = pd.DataFrame(
    {'pdb': filenames,
     'total_energy': total,
    })   
tot_e.to_csv('number_contacs4a.csv',index=False) 
print(tot_e)
print('succesful')
