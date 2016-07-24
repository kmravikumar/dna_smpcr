#!/usr/bin/python

## DNA PCR - modified SMPCR routine by D Chung, M Toney (Davis)
## - Krishna Ravikumar (07/23/2016)
## NOTE: Add Paper reference here .....
##
## USAGE: python smpcr.py Nmut Nrounds
## Nmut --  number of mutations and 
## Nrounds -- number of rounds
##


# import modules/scripts
import numpy as np
from numpy.random import randint
from pcr_reactions import *
from sys import argv, exit

# import progress bar module
BAR_FLAG=False
try:
    from progress.bar import Bar
    BAR_FLAG=True
except:
    pass

# set precision for printing output
np.set_printoptions(precision=2)


def set_primers(Nmut):
    """
    Set up initial mixture elements
    """
    primers = []
    for i in range(0,Nmut):
        ele = [i]
        primers.append(ele) 
    primers = np.array(primers).astype('int8')
    return primers


def are_lists_equal(a,b):
    """
    check if all elements of list a,b are same
    """
    if len(a) != len(b):
        return False
    for i,j in zip(a,b):
        if i!=j:
            return False
    #print a,b
    return True


def check_in_mix(product, mixture):
    """
    check if product is in mixture
    """
    idx = 0.1
    FLAG = False
    for i,ele in enumerate(mixture):
        if are_lists_equal(ele, product):
            idx = i
            FLAG = True
            break
    return idx, FLAG 


def do_round(mixture, mix_prop, primers, Nmut):
    """
    Add initial_mixture elements one-by-one
    to the current mixture and get new products
    """

    total = np.zeros(Nmut)
    next_mixture = []
    next_mix_prop = []
    if BAR_FLAG:
        bar = Bar('Processing', max=len(mixture))

    # for each template
    for i,template in enumerate(mixture):
        if BAR_FLAG: bar.next()
        # for each primer do first round PCR
        for primer in primers:
            product1 = reaction1(template,primer)
            product1_prop = mix_prop[i]
            # for each product from first round do secound round PCR
            for i2,template2 in enumerate(mixture):
                product2 = reaction2(template2,product1)
                product2_prop = product1_prop * mix_prop[i2] 

                # if product already exists then increment proportions
                idx, FLAG = check_in_mix(product2,next_mixture)
                if FLAG:
                    next_mix_prop[idx] += product2_prop
                else:
                    next_mixture.append(product2)
                    next_mix_prop.append(product2_prop)

                # depending on the size of the mutations
                # calculate the ratio of number of mutations
                total[len(product2)-1] += product2_prop

    if BAR_FLAG: bar.finish()
    Ntotal = np.sum(total)
    total = total/Ntotal
    # rescale proportions as probabilities
    next_mix_prop = np.array(next_mix_prop)/np.sum(next_mix_prop)
    #print Ntotal, total, len(next_mixture)

    return total, Ntotal, next_mixture, next_mix_prop


def print_output(pcr_round, ratio):
    """
    Print percentage of mutants 
    from each PCR round
    """
    
    print "\n"
    print "**************** ROUND {} *********************".format(pcr_round)
    print "Mutation ratios (%)"
    print ratio*100.0
    print "***********************************************\n"

    return 0



def do_PCR(Nmut,Nrounds):
    """
    For a given Nmut, Nrounds
    do the PCR rounds to find the 
    mutant distributions
    """
    Ntotal = Nmut*1.0
    total = np.zeros(Nmut)
    total[0] = 1.0
    print_output(1,total)
    
    primers = set_primers(Nmut)
    initial_mixture = np.copy(primers)
    initial_mix_prop = np.ones(len(initial_mixture))
    initial_mix_prop = initial_mix_prop/np.sum(initial_mix_prop)
    
    # repeat for rounds 2 - Nrounds
    for i in np.arange(2,Nrounds+1):
        total, Ntotal, next_mixture, next_mix_prop = \
            do_round(initial_mixture, initial_mix_prop, primers, Nmut)

        initial_mixture = list(next_mixture)
        initial_mix_prop = np.copy(next_mix_prop)
        print_output(i,total)

def print_help():
    """
    print help and exit
    """

    print "ERROR"
    print "USAGE: python smpcr.py Nmut Nrounds"
    print "Nmut --  number of mutations and "
    print "Nrounds -- number of rounds"

    exit(0)
    

if __name__=="__main__":
    """
    main function
    """
    
    if len(argv) < 2:
        print_help()

    # number of mutations
    Nmut = int(argv[1])
    # number of rounds
    Nrounds = int(argv[2])
    # call PCR rounds
    do_PCR(Nmut,Nrounds)
    


