#Apply a mask to a fasta file
#The mask contains [0-9] and possibly other characters (counts as -1)
#Replaces everything in the orginal fasta with Ns if it is below the integer level specified
#If no mask file is specificed, then it just produces an output that reads 0/1/2/3/4 according to 
#whether the input is uncalled (0) or called

from __future__ import division, print_function
import argparse, sys
from pyfaidx import Fasta
 
WRAP_LENGTH=50

IUPAC_ACGT=["A","C","G","T"]
IUPAC_HETS=["R", "Y", "S", "W", "K", "M"]
IUPAC_CALLS=IUPAC_ACGT+IUPAC_HETS
 
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

    if options.mask:
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
            new_seq="".join([x if y.isdigit() and int(y)>=options.level else "N" for x,y in zip(faseq,maskseq)])
            wrap_seq="\n".join(new_seq[i:i+WRAP_LENGTH] for i in range(0,len(new_seq), WRAP_LENGTH))
            print(">"+chrom)
            print(wrap_seq)
    else:
        for chrom in chroms: 
            print("Masking chromosome "+chrom, file=sys.stderr)
            faseq=fa[chrom][:].seq
            new_seq="".join(["1" if x in IUPAC_CALLS else "0" for x in faseq])
            wrap_seq="\n".join(new_seq[i:i+WRAP_LENGTH] for i in range(0,len(new_seq), WRAP_LENGTH))
            print(">"+chrom)
            print(wrap_seq)
        
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

