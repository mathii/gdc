#Apply a mask to a fasta file
#The mask contains [0-9] and possibly other characters (counts as -1)
#Replaces everything in the orginal fasta with Ns if it is below the integer level specified

from __future__ import division, print_function
import argparse, gdc, pdb, sys
from pyfaidx import Fasta
 
################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, default="", help=
                        "input fasta file")
    parser.add_argument('-m', '--mask', type=str, default="", help=
                        "Mask file (fasta)")
    parser.add_argument('-c', '--level', type=int, default=1, help=
                        "Filter level")

    return parser.parse_args()

################################################################################

def main(options):
    """
    print masked fasta
    """
    fa=Fasta(options.fasta)
    mask=Fasta(options.mask)

    if not set(fa.keys()).issubset(set(mask.keys())): 
        raise Exception("Mask does not include all chromosomes in fasta")
    chroms=[x for x in fa.keys()]
    
    for chrom in chroms:
        if len(fa[chrom][:]) != len(mask[chrom][:]):
            raise Exception("Chromosome "+chrom+" is different length in mask than fasta")
    
    for chrom in chroms: 
        print("Masking chromosome "+chrom, file=sys.stderr)
        faseq=fa[chrom][:].seq
        maskseq=mask[chrom][:].seq
        pdb.set_trace()
        new_seq="".join([x if y.isdigit() and int(y)>=options.level else "N" for x,y in zip(faseq,maskseq)])
        print("Printing chromosome "+chrom, file=sys.stderr)
        print(">"+chrom)
        print(textwrap.fill(new_seq,50))
        
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

