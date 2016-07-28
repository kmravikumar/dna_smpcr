#!/usr/bin/python

## SMPCR library method by Dong hee Chung and Michael D Toney (UC Davis)
## This program calculates the distribution of mutations for each round of SMPCR
## Code written by: Krishna Ravikumar (07/27/2016)
## For errors/comments/suggestions contact: mdtoney@ucdavis.edu
##
## NOTE: Add Paper reference here .....
##
## USAGE: python smpcr.py Nmut Nrounds [output_filename]
## Nmut --  number of mutations and 
## Nrounds -- number of rounds
## output_filename -- store output [default: smpcr.out]
##
## Output to outfile:    Round   %1mutation   %2mutations  . . .  %Nmutations   Unique_mutations
##
## Unique_combinations DOES NOT include the wild type gene, only mutations. 

# import modules/scripts
import numpy as np
from numpy.random import randint
from sys import argv, exit
import time
from itertools import combinations

# import progress bar module
BAR_FLAG=False
try:
    from progress.bar import Bar
    BAR_FLAG=True
except:
    pass

# Current date and time
now = time.strftime("%c")
print(now)

# set precision for printing output
np.set_printoptions(precision=1)



def print_output(pcr_round, ratio, mix_length, fp, HEADER_FLAG=0):
    """
    Output percentage of mutants from each PCR round
    """
    
    print("ROUND {} - {:5d} unique mutation combinations in mixture".format(pcr_round, mix_length))
    print("Mutation ratios (%)")
    print(ratio*100.0)
    print("\n")
    print("\n")

    if HEADER_FLAG:
        header = "Percentage of single-mutations, double-mutations, ... etc in the sample: \n"
        header += "--------------------------------------------------------------------\n"
        header += "Round "
        for i in range(0,len(ratio)):
            header += str(i+1)+"-Mut(%) "
        header += "Unique Mutations(#)\n"
        fp.write(header)

    wdata = "{:5d} ".format(pcr_round)
    for ele in ratio:
        wdata += "{:8.1f} ".format(ele*100) 
    wdata += "{:7d}\n".format(mix_length)
    fp.write(wdata)

    return 0


def print_input_parameters(Nmut, Nrounds, out_file, fp):
    """
    Print input parameters used 
    in the python script
    """
    summary = "********** "+ str(now) + " **********\n"
    summary += "OUTPUT from smpcr.py  \n"
    summary += "SMPCR library method by Dong hee Chung and Michael D Toney (UC Davis)\n"
    summary += "This program calculates the distribution of mutations for each round of SMPCR \n"
    summary += "Code written by: Krishna Ravikumar (07/23/2016) \n"
    summary += "For Errors/Comments contact : mdtoney@ucdavis.edu \n\n"
    summary += "Ref: Add Paper reference here ..... \n\n"
    summary += "USAGE: python smpcr.py Nmut Nrounds [output_filename]\n"
    summary += "Number of Mutation Sites:  "+str(Nmut)+"\n"
    summary += "Number of PCR Rounds:  "+str(Nrounds)+"\n"
    summary += "Output File Name:   "+out_file+"\n"
    summary += "**********************************************\n\n"
    fp.write(summary)
    return 0


def print_help():
    """
    Print help and exit
    """

    print("ERROR")
    print("USAGE: python smpcr.py Nmut Nrounds output_filename")
    print("Nmut --  number of mutations and ")
    print("Nrounds -- number of rounds")
    print("output_filename -- output file name [default: smpcr.out]")

    exit(0)
    

def print_dict(mut_dict, fp):
    """
    Print all mutations and 
    their proportions in each PCR round
    """
    
    header = "--------------------------------------------------------------------\n"
    header += "\n"
    header += "Mutation, " + " Percentage in PCR Round 1, 2, 3 ... etc ---->\n"
    header += "--------------------------------------------------------------------\n"
    fp.write(header)

    for mut in sorted(mut_dict, key=lambda mut_dict: len(mut_dict)):
        wdata = "{},".format(mut)
        for prop in mut_dict[mut]:
            wdata += "{:6.3f},".format(prop*100.0)
        wdata += "\n"
        fp.write(wdata)


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
    Check if all elements of list a,b are same
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
    Check if product is in mixture
    """
    idx = 0.1
    FLAG = False
    for i,ele in enumerate(mixture):
        if are_lists_equal(ele, product):
            idx = i
            FLAG = True
            break
    return idx, FLAG 


def reaction1(template, primer):
    """
    Megaprimer synthesis
    """
    temp = template[template > primer[0]]
    return np.hstack((primer,temp))


def reaction2(template, primer):
    """
    Mutagenesis PCR
    """
    temp = template[template<primer[0]]
    return np.hstack((temp,primer))


def create_all_combinations(Nmut, Nrounds):
    """
    Create all possible single, double, ...etc
    mutations for Nmutations
    """
    mut_dict = {}
    for i in range(1,Nmut+1):
        comb =  combinations(range(Nmut),i)
        for c in comb:
            mut_dict[c] = np.zeros(Nrounds)
    #print (len(mut_dict))
    return mut_dict

def normalize(mut_dict, i):
    """
    Normalize the population proportions in each round
    and find unique mutations
    """

    S = 0.0
    unique = 0
    for ele in mut_dict:
        val = mut_dict[ele][i]
        if val !=0: unique += 1
        S += mut_dict[ele][i]
    for ele in mut_dict:
        mut_dict[ele][i] /= S
    return unique


def do_round(primers, Nmut, mut_dict, iround):
    """
    Add initial_mixture elements one-by-one to the current mixture and get new products
    """

    total = np.zeros(Nmut)
    if BAR_FLAG:
        bar = Bar('  Calculating next round:   ', max=len(mut_dict))

    # for each template
    for templatet in mut_dict:
        if BAR_FLAG: 
            bar.next()
        template = np.array(templatet)
        product1_prop = mut_dict[templatet][iround-2] 
        if product1_prop ==0: continue

        # for each primer do first round PCR
        for primer in primers:
            product1 = reaction1(template,primer)

            # for each product from first round do second round PCR
            for template2t in mut_dict:
                template2 = np.array(template2t)
                product2 = reaction2(template2,product1)
                product2_prop = product1_prop * mut_dict[template2t][iround-2] 
                if product2_prop == 0: continue
                product2t = tuple(product2)
                mut_dict[product2t][iround-1] += product2_prop

                # depending on the size of the mutations
                # calculate the ratio of number of mutations
                total[len(product2)-1] += product2_prop

    if BAR_FLAG: 
        bar.finish()
    Ntotal = np.sum(total)
    total = total/Ntotal
    # rescale proportions as probabilities
    unique =  normalize(mut_dict,iround-1)

    #print Ntotal, total, len(next_mixture)
    return total, unique


def do_PCR(Nmut,Nrounds,out_file):
    """
    Do the PCR rounds to find the mutant distributions
    """
    fp = open(out_file,"w")
    # Print summary of input parameters
    print_input_parameters(Nmut, Nrounds, out_file, fp)
    Ntotal = Nmut*1.0
    total = np.zeros(Nmut)
    total[0] = 1.0
    print_output(1,total,Nmut,fp,1)

    mut_dict = create_all_combinations(Nmut, Nrounds)   
    for i in range(0,Nmut):
        mut_dict[(i,)][0] = 1.0
    unique = normalize(mut_dict,0)

    primers = set_primers(Nmut)
    
    # repeat for rounds 2 to Nrounds
    for i in np.arange(2,Nrounds+1):
        total, unique = \
            do_round(primers, Nmut, mut_dict, i)
        # print single, double..etc mutation summary
        print_output(i,total,unique,fp)

    # print proportions in each round
    print_dict(mut_dict,fp)

    fp.close()
    return 0



if __name__=="__main__":
    """
    Main function
    """

    if len(argv) < 3:
        print_help()
    out_file = "smpcr.out"
    if len(argv) == 4:
        out_file = argv[3].strip()

    # number of mutations
    Nmut = int(argv[1])
    # number of rounds
    Nrounds = int(argv[2])
    # call PCR rounds
    do_PCR(Nmut,Nrounds,out_file)
    

