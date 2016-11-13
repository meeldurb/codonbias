#!usr/bin/env python

"""
Author: Melanie van den Bosch
Script for calculating oligo properties
"""

# import statements
from __future__ import division
import math

# functions are below

# probably the input sequence is read as a string  
## class OligoProperties(self):
## i can probably make a class of this one
def count_aminoacids(sequence):
    A = sequence.count("A")
    C = sequence.count("C")
    G = sequence.count("G")
    T = sequence.count("C")
    return A, C, G, T

def molecular_weight(sequence):
    """ Returns the molecular weight in Da of a sequence
    """
    # weight constants for aminoacids
    A_WEIGHT = 313.21
    C_WEIGHT = 289.18
    G_WEIGHT = 329.21
    T_WEIGHT = 304.2
    # counting the aminoacids
    A, C, G, T = count_aminoacids(sequence)
    # calculate total weight of sequence
    # This is the Anhydros Mol Weight as retrieved from 
    # http://biotools.nubic.northwestern.edu/OligoCalc.html
    seq_weight = A * A_WEIGHT + C * C_WEIGHT + G * G_WEIGHT + \
                 T * T_WEIGHT - 61.95
    return seq_weight

def GC_content(sequence):
    """ Returns the GC content in % of a sequence
    """
    # counting the aminoacids
    A, C, G, T = count_aminoacids(sequence)
    # dividing frequency GC by total amino acids    
    GC_percentage = G+C/length(sequence)
    return GC_percentage
    
def length(sequence):
    """ Returns the length in numbers of a sequence
    """
    result = len(sequence)
    return result

def melting_temp(sequence, Na_conc=None):
    # needs some work still
    """ calculates the melting temperature of a sequence
    by taking the other parameters in regard
    sequence -- the sequence that you want to calculate the Tm of
    Na_conc -- the concentration of Na+(mM) present
    """
    # counting the aminoacids
    A, C, G, T = count_aminoacids(sequence)
    # separate calculations for each length
    if Na_conc == None:
        if len(sequence) <= 13:
            Tm = (A+T)*2 + (G+C)*4
        else:
            Tm = 64.0 + 41*((G+C - 16.4)/(A+T+G+C))
    if Na_conc != None:
        if len(sequence) <= 13:
            Tm = (A+T)*2 + (G+C)*4 - (16.6*math.log(0.050) + \
                 16.6*math.log(Na_conc, 10))
        if len(sequence) >= 14:
            Tm = 100.5 + (41*(G+C)/(A+T+G+C) - (820/(A+T+G+C) + \
                      16.6*math.log(Na_conc, 10)))
    return Tm



if __name__== "__main__":
    sequence = "ACTGCCGTAGGCTACCCAGT"
    
