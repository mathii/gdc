#Given a vcf and a panel file, export a frequency file with N+5 columns and M+1 rows
# N populations and M SNPs - one header row and the first 5 columns are 
#SNPID, CHR, POS, REF, alter_code1

from __future__ import division, print_function
from collections import defaultdict
import argparse, sys
import pdb


################################################################################

def parse_options():
    """
    Try using argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default="-")
    parser.add_argument('-p', '--panel', type=str, action="store", default=None, help=
                        "Two column file mapping IDs to populations")

    return parser.parse_args()

################################################################################

def read_panel(panel_file):
    """
    Open the panel file and return a dictionary that maps sample to popuation
    """
    map={}
    pf=open(panel_file, "r")
    for line in pf:
        bits=line.split()
        map[bits[0]]=bits[1]
        
    pops=list(set(map.values())) 
    pops.sort()   
    
    return map, pops

################################################################################

def main(options):
    """
    run through the file and output to stdout. 
    """
    map, pops=read_panel(options.panel)
    print("\t".join(["SNPID", "CHR", "POS", "REF", "ALT"]+pops))
    
    samples=[]
    for line in options.input:
        if line.startswith("##"):
            pass
        elif line.startswith("#"): #Header line
            samples=line.split()[9:]
        else: #This is a data line. 
            counts=defaultdict(int)
            totals=defaultdict(int)
            bits=line.split()
            first_cols=[bits[2], bits[0], bits[1], bits[3], bits[4]]
            for i in range(9, len(bits)):
                if samples[i-9] not in map:
                    continue
                for j in [0,1]:
                    if bits[i][j]=="0":
                        totals[map[samples[i-9]]]+=1
                    if bits[i][j]=="1":
                        counts[map[samples[i-9]]]+=1
                        totals[map[samples[i-9]]]+=1
            
            freqs=[0.0]*len(pops)
            try:
                for i,pop in enumerate(pops):
                    freqs[i]=counts[pop]/totals[pop]
                print("\t".join(first_cols+[format(x, "1.4f") for x in freqs]))
            except ZeroDivisionError:
                continue
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)
