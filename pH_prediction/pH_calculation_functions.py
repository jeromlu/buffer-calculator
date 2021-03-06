# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 20:49:12 2019

@author: jeromlu2
"""
import numpy as np

Kw = 10**-14

def pKa_effective(pKa, z_a, I, A = 0.5114, b =0.1):
    '''
    Calculates pKa effective from pKa. 
    Source: Pabst and Carta 2007.
    
    Parameters:
    -----------
    pKa : float
        pKa of acid, from literature
    z_a : int
        charge of conjugate acid species
    I : float
        ionic strength of solution
    A : float
        A is a constant that incorporates temperature effects (A = 0.5114 at 298 K)
    b : float
        b is a buffer-dependent parameter
        
    Returns:
    --------
    pKa_eff : float
        corrected pKa from extended Debye–Huckel theory
        
    '''
    pKa_eff = pKa +(z_a - 1) * ( (A * np.sqrt(I)) / (1 + np.sqrt(I)) - b * I)
    return pKa_eff

def ionic_strength(charges, concentrations):
    '''
    From chagres and concentrations array calculates ionic strength.
    
    Parameters:
    -----------
    charges : np.array of floats
        charges of individual species
    concentrations : np.array of floats
        concentrations of individual speicies [mM]
    Returns:
    --------
    float ionic strength [mM]
    '''
    
    I = 0.5 * (charges**2 * concentrations).sum()
    return I


def acid_conc_slow(Ca, pKa, z_a, pH, I):
    
    #number of chemical species in solution
    n = len(pKa)
    
    #calculation of H+ and OH- concentration
    c_H = 10**-pH
    c_OH = Kw / c_H
    
    temp_H = np.power(c_H, -np.arange(n))    
    
    Ka_eff = 10**pKa_effective(pKa, z_a, I)
    
    denum = Ka_eff * np.array([1, Ka_eff[0], Ka_eff[1]]) * np.array([1, 1, Ka_eff[1]]) * temp_H 
    
    conc = Ca /(1 + denum.sum())
    
    return conc

def acid_conc(Ca, pKa, z_a, c_H, I):
    
    n = len(pKa)
    conc = np.ones(n)
    
    #number of chemical species in solution
    #n = len(pKa)
    
    #calculation of H+ and OH- concentration
    
    
    #Ka_eff_0 = 10**pKa_effective(pKa[0], z_a[0], I)
    #Ka_eff_1 = 10**pKa_effective(pKa[1], z_a[1], I)
    #Ka_eff_2 = 10**pKa_effective(pKa[2], z_a[2], I)
    
    Ka_eff = 10.**-pKa_effective(pKa, z_a, I)
    
    if n == 3:
        denum = 1 + Ka_eff[0] / c_H + Ka_eff[0] * Ka_eff[1] / c_H**2 + Ka_eff[0] * Ka_eff[1] * Ka_eff[2] / c_H**3
    elif n == 2:
        denum = 1 + Ka_eff[0] / c_H + Ka_eff[0] * Ka_eff[1] / c_H**2
    else:
        denum = 1 + Ka_eff[0] / c_H
    
    conc[-1] = Ca / denum
    for i, Ka in enumerate(Ka_eff):
     conc[i] = conc[i-1] * Ka / c_H
    
    return conc



def solution_charge(charges, concentrations, pH, c_Cl, c_Na):
    
    c_H = 10.**-pH
    c_OH = Kw / c_H     
    sol_charge = charges * concentrations + c_H - c_OH - c_Cl + c_Na
    
    return sol_charge
    
    





