#eigenstrat (packed or unpacked) to vcf
#Writes to stdout, unlike many of these scripts. 
#Usage: python eigenstrat2vcf.py -i root 
#Data files are root.snp, root.ind and root.geno

from __future__ import division, print_function
import argparse, gdc, pyEigenstrat

#Remember, in eigenstrat, 2 means "2 ref copies"
GT_DICT={2:"0/0", 1:"0/1", 0:"1/1", 9:"./."}

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, default="", help=
                        "Root for eigenstrat files - i.e {root.snp, root.geno, root.ind}")

    return parser.parse_args()

################################################################################

def main(options):
    """
    Convert
    """

    data=pyEigenstrat.load(options.input)

    #Write header. 
    print("##fileformat=VCFv4.0")
    print("##source=eigenstrat2vcf.R")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(data.ind["IND"]))

    #Now line by line write data
    for i,d in enumerate(data):
        this_snp=data.snp[i]
        line="\t".join([this_snp["CHR"], str(this_snp["POS"]), this_snp["ID"], 
                         this_snp["REF"], this_snp["ALT"], "100", "PASS", ".", "GT" ])
        line=line+"\t"+"\t".join([GT_DICT[x] for x in d])
        print(line)
        
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

