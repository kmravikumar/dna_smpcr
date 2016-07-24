#!/usr/bin/python

from pcr_reactions import * 
import numpy as np

def test_reaction1():
    """
    function to test reaction1
    with random primers and templates
    """
    Nmut = randint(4,10,1)
    a = np.arange(0,Nmut)
    
    t = randint(0,2,Nmut)
    template = a[t==1]

    primer = randint(0,Nmut,1)
    if len(primer) == 0 or len(template)==0:
        return 0
    
    print template, "  ", primer, "  ", reaction1(template,primer)
    

def test_reaction2():
    """
    function to test reaction2
    with random primers and templates
    """
    Nmut = randint(4,10,1)

    a = np.arange(0,Nmut)
    t = randint(0,2,Nmut)
    template = a[t==1]

    a = np.arange(0,Nmut)
    p = randint(0,2,Nmut)
    primer = a[p==1]
    if len(primer) == 0 or len(template)==0:
        return 0

    print template, "  ", primer, "  ", reaction2(template,primer)
    

def do_tests():
    """
    Do reaction1, reaction2 tests
    """
    print "Reaction1"
    print "Template  ", "Primer  ", "Product"    
    for i in range(0,30):
        test_reaction1()

    print "Reaction2"
    print "Template  ", "Primer  ", "Product"    
    for i in range(0,30):
        test_reaction2()


if __name__=="__main__":
    do_tests()
