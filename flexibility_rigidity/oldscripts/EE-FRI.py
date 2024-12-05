#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#-------------Flexibility Regidity Index---------------------------#
#-------------Multiple kernel(EE) with different hyperparameters---#
# Navya Ramesh


# import necessary modules
from biopandas.mol2 import PandasMol2
import numpy as np
from biopandas.pdb import PandasPdb
ppdb = PandasPdb()
import math
import re 
import os
import sys

# read pdbs and mol files from V2007
file_list=os.listdir(r"/mnt/home/tippanan/prj/v2007")
#print (file_list)
#print((file_list))
import itertools
import pandas as pd
import numpy as np

#read the train and test csv files
train_data = pd.read_csv(r"/mnt/home/tippanan/prj/v2007_refine_listt.csv") 
test_d = pd.read_csv(r"/mnt/home/tippanan/prj/v2007_core_listt.csv")
train_data = np.asarray(train_data)
test_d = np.asarray(test_d)
train_list = train_data[:,0]
test_list = test_d[:,0]
print(test_list)
mol_lt = []
prt_lt = [ ]
dist = [ ]
eta = [ ]
eta_lst = [ ]

#array to store data
exp = np.zeros((1106,36*3))
res2 = [ ]
res3 = [ ]
lor = [ ]
prt_lt = []
ml_lt = []
dist = 0

#specifing atoms and Vanderwal's radius
lig_ele_rad = [1.70, 1.55, 1.52, 1.80, 1.80, 1.47, 1.75, 1.85, 1.98]
lig_ele_list = ['C','N','O','S','P','F','Cl','Br','I']
r_i = 0
t1 = 1
k1 = 5
k2 = 5
t2 = 12.5
v = 40
counter = 0
car = {"C": 1.70,  "N": 1.55,  "O": 1.52,  "S": 1.80, "P":1.80 ,"F":1.47 ,"Cl":1.75,"Br":1.85,"I":1.98}

ftr_df = {"CC":0,"CN":1,"CO":2,"CS":3,"CP":4,"CF":5,"CCl":6,"CB":7,"CI":8,
          "NC":9,"NN":10,"NO":11,"NS":12,"NP":13,"NF":14,"NCl":15,"NB":16,"NI":17,
          "OC":18,"ON":19,"OO":20,"OS":21,"OP":22,"OF":23,"OCl":24,"OB":25,"OI":26,
          "SC":27,"SN":28,"SO":29,"SS":30,"SP":31,"SF":32,"SCl":33,"SB":34,"SI":35 }

# read the pdb and mol files 
for i in range(len(train_list)):
    prt_lt = []
    ml_lt = []
    counter = counter + 1
    pmol = PandasMol2().read_mol2('/mnt/ufs18/home-023/tippanan/prj/v2007/'+train_list[i]+'/'+train_list[i]+'_ligand.mol2')
   # print(pmol.df.head(10))
    mol_lst = pmol.df[pmol.df['atom_type'] != 'H' ][['atom_name','x', 'y', 'z']]
    ppdb = PandasPdb()
    pbd_re = ppdb.read_pdb('/mnt/ufs18/home-023/tippanan/prj/v2007/'+train_list[i]+'/'+train_list[i]+'_protein.pdb')
    prot_lst =ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] != 'H'][['atom_name','x_coord', 'y_coord', 'z_coord']]
    prot_lst  = np.asarray(prot_lst)
    mol_lst = np.asarray(mol_lst)
    for x, y in itertools.product(prot_lst,mol_lst):
        
        # computer euclidean distance between atoms
        dist = np.linalg.norm(x[1:4]-y[1:4], axis = 0)
        if dist< 39:
            if (x[0][0] in lig_ele_list ) and (y[0][0] in lig_ele_list):
                
                # get the radius based on the atoms and get exponential kernel values
                r_i = car[x[0][0]]
                r_j = car[y[0][0]]
                p_char = x[0][0]
                m_char = y[0][0]
                print(r_i,x[0][0],x)
                print(r_j,y[0][0],y)
                ftr_cc = p_char+m_char
                df_num = ftr_df[ftr_cc]
                print(ftr_cc)
                print(df_num)
                eta = (t1*(r_i+r_j))
                #print("dist is",dist)
                lor_v = 1//(1+(dist/eta)**v1)
                #exp[i,36+df_num] = exp[i,36+df_num] +lor_v
                print("row is", i)
                #exp_r = math.exp(-(dist/eta)**k1)
                exp[i,df_num] = exp[i,df_num]+lor_v
    for x, y in itertools.product(prot_lst,mol_lst):
        dist = np.linalg.norm(x[1:4]-y[1:4], axis = 0)
        if dist< 47:
            if (x[0][0] in lig_ele_list ) and (y[0][0] in lig_ele_list):
                r_i = car[x[0][0]]
                r_j = car[y[0][0]]
                p_char = x[0][0]
                m_char = y[0][0]
                #print(r_i,x[0][0],x)
                #print(r_j,y[0][0],y)
                ftr_cc = p_char+m_char
                df_num = ftr_df[ftr_cc]
                #print(ftr_cc)
                #print(df_num)
                eta = (t2*(r_i+r_j))
                lor_v = 1//(1+(dist/eta)**v2)
                #print("dist is",dist)
                #exp_r = math.exp(-(dist/eta)**k2)
                exp[i,36+df_num] = exp[i,36+df_num]+lor_v
                #exp[i,36+df_num] = exp[i,36+df_num] +lor_v
                print("row is", i)
                print("column is",df_num+36)
print(exp)
np.savetxt("Lor_2kerneltrain_2.5-1-12_2.5-5.5-21.csv", exp, delimiter=",")

