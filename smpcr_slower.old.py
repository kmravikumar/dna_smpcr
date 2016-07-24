#!/usr/bin/python

import numpy as np
from numpy.random import randint
from pcr_reactions import *
#from plot import plot_template
from sys import exit


def set_primers(Nmut):
    """
    Set up initial mixture elements
    """
    primers = []
    for i in range(0,Nmut):
        ele = [i]
        primers.append(ele) 
    primers = np.array(primers).astype('int8')
    Ntotal = Nmut*1.0
    total = np.zeros(Nmut)
    total[0] = 1
    print Ntotal, total
    return primers


def do_round(mixture, primers, Nmut):
    """
    Add initial_mixture elements one-by-one
    to the current mixture and get new products
    """

    total = np.zeros(Nmut).astype('int32')
    next_mixture = []
    for template in mixture:
        for primer in primers:
            product1 = reaction1(template,primer)
            for template2 in mixture:
                product2 = reaction2(template2,product1)
                next_mixture.append(product2)
                total[len(product2)-1] += 1
    
    Ntotal = np.sum(total)*1.0
    total = total/Ntotal
    print Ntotal, total
    return total, Ntotal, next_mixture


Nmut = 3
Nrounds = 4

primers = set_primers(Nmut)
initial_mixture = np.copy(primers)
print primers
print initial_mixture

for i in np.arange(0,Nrounds-1):
    total, Ntotal, next_mixture = do_round(initial_mixture, primers, Nmut)
    initial_mixture = list(next_mixture)





