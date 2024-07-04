# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 12:49:58 2021

@author: Divya
"""

import pandas as pd
from tqdm import tqdm
#from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from scipy.stats import pearsonr
from itertools import combinations

dataset = pd.read_csv('input_regres_train.csv')
testset = pd.read_csv('input_test.csv')

trainingset=dataset

def rSubset(arr, r):
    return list(combinations(arr, r))
arr =[i for i in range(2,37)]

"""if __name__ == "__main__":
    arr =[i for i in range(5,100)]
    #arr=[1,2,3,4]
    print (arr, len(arr))
    r = 5
    comb4=rSubset(arr, r)
print (len(comb4))
"""
data=[]
for i in tqdm(rSubset(arr,6)):
    lis=list(i)
    names=[dataset.columns[j] for j in i]
    X = trainingset.iloc[:,lis]
    Y = trainingset.iloc[:,1]
    X_t = testset.iloc[:,lis]
    Y_t = testset.iloc[:,1]

    regressor = LinearRegression()  
    regressor.fit(X,Y)
    Y_train_pred = regressor. predict(X)
    m = pearsonr(Y, Y_train_pred)[0]
    n = metrics.mean_absolute_error(Y, Y_train_pred)

    Y_test_pred = regressor.predict(X_t)
    o = pearsonr(Y_t, Y_test_pred)[0]
    p = metrics.mean_absolute_error(Y_t, Y_test_pred)
   
    data.append((names,m,n,o,p))

results=pd.DataFrame(data,columns=["A","B","C","D","E"])
# results.to_csv(r'4_check.txt', header=None, index=None, sep='\t', mode='w')
results.to_csv('6_comb.csv', header=None, index=None, mode='w')
