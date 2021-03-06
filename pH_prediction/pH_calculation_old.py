# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:07:06 2019

@author: jeromlu2
"""
import numpy as np
import pandas as pd

from solvers import solve_quadrat_eq

def Q1(x, Ca, z, Ka):
    return (- Ca * z[0] * Ka[0]) / (x + Ka[0])

def dQ1(x, Ca, z, Ka):
    pass
    
    
def Q2(x, Ca, z, Ka):
    
    nom = -Ca * (z[0] * Ka[0] * x + z[1] * Ka[1])
    denom = x**2 + Ka[0] * x + Ka[1] 
    return nom / denom

def Q3(x, Ca, z, Ka):
    
    nom = -Ca * (z[0] * Ka[0] * x**2 + z[1] * Ka[1] * x + z[2] * Ka[2])
    denom = x**3 + Ka[0] * x**2 + Ka[1] * x + Ka[2] 
    return nom / denom


def acid_pH(Ca, Ka, Kw = 10**(-14), z = -1):
    """
    Calculates the pH and pOH of the acid.
    
    Parameters
    ----------
    Ca : float
        concentration of the main component (acid) in [mM]
    Ka : float
        equilibrium constant of the acid component in [mM]
        """
    p = np.zeros(4)

    #coefficient of the 
    p[0] = 1
    p[1] = Ka
    p[2] = (z * Ka * Ca / 1000 + Kw)
    p[3] = - Ka * Kw

    sol = np.roots(p)
    print('\n','sol ', sol)
    mask = (sol > 0) & (np.imag(sol) == 0)
    if mask.sum() >= 1:
        H_conc = sol[mask][0]
    else:
        H_conc = 1
        print('no solution')
    pH = - np.log10(H_conc)

    return pH
acid_pH = np.vectorize(acid_pH)

def strong_acid_pH(Ca, Ka, Kw = 10**(-14), z = -1):
    """
    Calculates the pH and pOH of the acid.
    
    Parameters
    ----------
    Ca : float
        concentration of the main component (acid) in [mM]
    Ka : float
        equilibrium constant of the acid component in [mM]
        """
        
    #coefficient of the 
    a = 1
    b = z * Ca / 1000
    c = Kw
    
    sol1, sol2 = solve_quadrat_eq(a, b, c)

    sol = np.array([sol1, sol2])
    print(sol)
    mask = (sol > 0)
    if mask.sum() >= 1:
        H_conc = sol[mask][0]
    else:
        H_conc = 1
        print('no solution ')
    pH = - np.log10(H_conc)
    return pH

strong_acid_pH = np.vectorize(strong_acid_pH)

def strong_base_pH(Cb, Kb, Kw = 10**(-14)):
    
    #coefficient of the 
    a = 1
    b = - Cb / 1000
    c = - Kw
    
    sol1, sol2 = solve_quadrat_eq(a, b, c)

    sol = np.array([sol1, sol2])
    print(sol)
    mask = (sol > 0)
    if mask.sum() >= 1:
        H_conc = sol[mask][0]
    else:
        H_conc = 1
        print('no solution ')
    pH = 14 + np.log10(H_conc)
    return pH
strong_base_pH = np.vectorize(strong_base_pH)

def super_strong_acid(Ca, Ka, Kw = 10**(-14)):
    return - np.log10(Ca/1000)
super_strong_acid = np.vectorize(super_strong_acid)

if __name__ == "__main__":
    
    Ka = 10**(-7.20)
    Ca = 1000 #mM
    
    pH = acid_pH(Ca, Ka)
    #print('pH = ', pH)
    
    #hydrochloric acid pKa = -4
    Ka = 10**(-13.5)
    Kb = 10**(-0.5)
    Ca = 1000 #M
    
    
    Ca_array = np.array([1000., 500., 100., 10., 0.1, 0.01, 0.001, 0.0001, 0.00005, 0.00001, 0.000001]) #mM
    df = pd.DataFrame([], index = Ca_array)
    #df['exact'] = acid_pH(Ca_array, Ka)
    df['strong'] = strong_acid_pH(Ca_array, Ka)
    df['super'] = super_strong_acid(Ca_array, Ka)
    df['strong_base'] = strong_base_pH(Ca_array, Kb)
    
    
"""    
    #Nitric acid pKa = -4
    Ka = 10**(1.)
    Ca = 1 #mM
    
    pH, pOH = acid_pH(Ca, Ka)
    print('pH = ', pH)
    
    pH, pOH = strong_acid_pH(Ca, Ka)
    print('pH strong = ', pH)
    
    #Acetic acid pKa = -4
    Ka = 10**(-4.75)
    Ca = 1 #mM
    
    pH, pOH = acid_pH(Ca, Ka)
    print('pH = ', pH)
    
    pH, pOH = strong_acid_pH(Ca, Ka)
    print('pH strong = ', pH)
    
    #Hydrocyanic acid pKa = -4
    Ka = 10**(-9.31)
    Ca = 10 #mM
    
    pH, pOH = acid_pH(Ca, Ka)
    print('pH = ', pH)
"""