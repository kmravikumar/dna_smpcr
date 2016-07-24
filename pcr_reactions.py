#!/usr/bin/python

# PCR reactions for each round
import numpy as np

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
