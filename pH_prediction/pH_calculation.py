# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 20:21:55 2019

@author: jeromlu2
"""
#%%python packages

#third party packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#my packages
from solvers import bisekcija
from pH_calculation_functions import acid_conc

Kw = 10**(-14)

fname_data = './data.xlsx'

data_mol = pd.read_excel(fname_data, 
                         sheet_name = 'components_init_mol', index_col = 1)

data_comp = pd.read_excel(fname_data, sheet_name = 'comp_dis_sub', index_col = 1)

pKa = {}
charge = {}
z_a = {}

for acid_ID in np.unique(data_comp.acid_ID):
    mask = acid_ID == data_comp.acid_ID
    pKa[acid_ID] = data_comp.loc[mask,'pKa'].values
    charge[acid_ID] = data_comp.loc[mask,'charge'].values
    z_a[acid_ID] = data_comp.loc[mask,'conjug_charge'].values



comp_conc = np.empty(len(data_comp))

pH = 5
phos_type = 2.

min_Na_conc = data_mol.loc['NaCl', 'molarity [mM]'] + \
                phos_type * data_mol.loc['Phosphate', 'molarity [mM]'] + \
                data_mol.loc['EDTA', 'molarity [mM]']
min_Cl_conc = data_mol.loc['NaCl', 'molarity [mM]'] + \
                data_mol.loc['Arginine', 'molarity [mM]'] + \
                data_mol.loc['Cysteine', 'molarity [mM]']

I_prev = -1.
I_current = (min_Na_conc + min_Cl_conc) * 0.5

def calculate_charge_IS(data_mol, pH, I_prev, c_Cl, c_Na):
    c_H = 10.**-pH
    c_OH = Kw / c_H
    
    I_current = (c_Cl + c_Na) * 0.5
    total_charge = c_H - c_OH - c_Cl + c_Na
    
    for acid_ID in data_mol.index[:-1]:
        Ca =  data_mol.loc[acid_ID, 'molarity [mM]']
        if Ca > 0:
            conc = acid_conc(Ca, pKa[acid_ID], z_a[acid_ID], c_H, I_prev)
            sum_term_charge = charge[acid_ID] * conc
            #sum_term = sum_term_charge * charge[acid_ID]
            #I_current = I_current + 0.5 * sum_term.sum()
            
            total_charge = total_charge + sum_term_charge.sum()
    print('I_current ', I_current)
    print('Charge ', total_charge)
    return total_charge, I_current


c_Cl0 = min_Cl_conc
c_Na0 = min_Na_conc
I_current = (min_Na_conc + min_Cl_conc) * 0.5
total_charge = 1

#calulate required c_Cl or c_Na
total_charge0, I_current0 = calculate_charge_IS(data_mol, pH, I_current, c_Cl0, c_Na0)    

if total_charge < 0:
    c_Na = c_Na + abs(total_charge)
elif total_charge > 0:
    c_Cl = c_Cl + abs(total_charge)
else:
    print('Exact match with zero!')

points = np.empty((20,2))
c_Na = 50. * c_Na0 + 0.03
c_Nas = np.linspace(5, 28, 20)

for i, c_Na in enumerate(c_Nas):
     total_charge, I_current= calculate_charge_IS(data_mol, pH, I_current, c_Cl, c_Na)
     points[i,:] = total_charge, I_current
    
plt.plot(c_Nas, points[:, 0])
plt.plot(c_Nas, points[:, 1])
    
#%%
c_Na = 2. * c_Na0 + 0.03
total_charge, I_current = calculate_charge_IS(data_mol, pH, I_current0, c_Cl, c_Na0)


slope = (total_charge - total_charge0) / (c_Na - c_Na0)
c_Na0 = c_Na
c_Na = slope * total_charge + c_Na
if c_Na > 0:
    c_Na = 0
total_charge0, I_current0 = total_charge, I_current
total_charge, I_current = calculate_charge_IS(data_mol, pH, I_current, c_Cl, c_Na0)

slope = (total_charge - total_charge0) / (c_Na - c_Na0)
c_Na0 = c_Na
c_Na = slope * total_charge + c_Na
if c_Na > 0:
    c_Na = 0
total_charge0, I_current0 = total_charge, I_current
total_charge, I_current = calculate_charge_IS(data_mol, pH, I_current, c_Cl, c_Na0)


I_current = I_current + 0.5 * abs(total_charge)


    
n = 0




#%%
while abs(I_prev - I_current) > 0.00005 and total_charge > 0.002 and n < 200:
    I_prev = I_current
    total_charge, I_current, c_Cl, c_Na = calculate_charge_IS(data_mol, pH, I_current, c_Cl, c_Na)
    delta = abs(I_prev - I_current)
    #print(n, ' prev, current ', I_prev, I_current)
    n = n+1
    
    
    
#%%***************** Ploting concentrations**************************************    
acid_ID = 'Phosphate'

pHs = np.linspace(1, 13, 50)

concentrations = np.empty((len(pHs),len(pKa[acid_ID])))

for i, pH in enumerate(pHs):
    conc = acid_conc(50, pKa[acid_ID], z_a[acid_ID], 10**-pH, 10)
    concentrations[i,:]=conc


plt.plot(pHs, concentrations)
    

for i, pH in enumerate(pHs):
    conc = acid_conc(50, pKa[acid_ID], z_a[acid_ID], 10**-pH, 18)
    concentrations[i,:]=conc


plt.plot(pHs, concentrations)



























