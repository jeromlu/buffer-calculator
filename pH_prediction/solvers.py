# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 20:06:37 2019

@author: JEROMLU2
"""
import numpy as np


def solve_quadrat_eq(a0, b0, c0, eps = 10**(-12), verbose = True):
    
    if verbose:
        print('a ', a0,' b ', b0, ' c ', c0)
    
    max_abs = 0
    for val in [a0, b0, c0]:
        if abs(val) > max_abs:
            max_abs = abs(val)
        
        
    a = a0 / max_abs 
    b = b0 / max_abs
    c = c0 / max_abs
    
    D = b**2 - 4*a*c
    
    if abs(b) < eps:
        D = 0
    
    if D < 0:
        return 1, 1
    
    
    if D > 0:
        
        sol1 = (- b - np.sqrt(D) ) / (2 * a)    
        sol2 = c / (a * sol1)
    else:
    
        sol2 = (- b + np.sqrt(D)) / (2 * a)
        sol1 = c / (a * sol2)
    
    if verbose:
        test =  ( sol1 * a + b) + c / sol1
        print(test)
        test = sol2**2 * a + sol2 * b + c 
        print(test)
        


    
    return sol1, sol2

def bisekcija(func, a, b, eps = 10**(-3)):
    
    if np.sign(func(a)) == np.sign(func(b)):
        print('Function has same sign in both starting points.')
        return
    
    s = b - a
    n = 0
    while (abs(s) > eps) or (n > 200):
        
        s = s/2
        
        c = a + s
        if np.sign(func(a)) == np.sign(func(c)):
            a = c
        else:
            b = c
        
        n = n + 1
    print('Koraki: ', n)    
    return c

if __name__ == '__main__':
    sol1, sol2 = solve_quadrat_eq(1.2345678, 76543210.5, 0.1122334455)
    
    print(sol1, sol2)
    