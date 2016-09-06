#This is for converting Shop Mallick's polysite format to vcf
#Probably only useful to you if you are working on the SGDP 
#Very specific to this particular format.

from __future__ import division, print_function
import argparse, sys, pdb

#Remember, in eigenstrat, 2 means "2 ref copies"

CODES={
       "A":"AA",
       "C":"CC",
       "G":"GG",
       "T":"TT",
       "R":"AG",
       "Y":"CT",
       "S":"GC",
       "W":"AT",
       "K":"GT",
       "M":"AC",
       "-":"..",
       "N":"..",
    }

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), default="-")
    parser.add_argument('-c', '--chrom', type=str, default="")

    return parser.parse_args()

################################################################################

def main(options):
    """
    Convert
    """
    samples=[]    
    include_ancients=False
    include_refs=False
    
    reading_header=True
    for line in options.input:
        if len(line)==1:
            continue
        elif line[:2]=="##" and reading_header:
            bits=line.split()
            if len(bits)<4:
                continue
            elif bits[1]!="..":
                continue
            elif bits[2][0]=="3" and include_refs: #Refs
                samples.append(bits[3])
            elif bits[2][0]=="4": #C team
                samples.append(bits[7])
            elif bits[2][0]=="5": #B team with cteam processing
                samples.append(bits[7])
            elif bits[2][0]=="7" and include_ancients: #Ancients
                samples.append(bits[4].split(":")[0])
            elif bits[2][0]=="8": #A/B team, original
                samples.append(bits[4].split(":")[0])
        elif line[0]=="#" and reading_header:
            reading_header=False
            print("##fileformat=VCFv4.2")
            print("##source=polysites2vcf.py")
            print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
            print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(samples))
        elif not reading_header:
            bits=line.split()
            chrom=bits[0]

            if options.chrom and options.chrom!=chrom:
                continue

            poss=bits[1]
            idd=chrom+"_"+poss
            ref=bits[2][0]
            
            alleles=(bits[3]+bits[4]+bits[7]).upper()
            if include_refs and include_ancients:
                alleles=(bits[2]+bits[3]+bits[4]+bits[6]+bits[7]).upper()
            elif include_refs:
                alleles=(bits[2]+bits[3]+bits[4]+bits[7]).upper()
            elif include_ancients:
                alleles=(bits[3]+bits[4]+bits[6]+bits[7]).upper()
            
            gts = [CODES[x] for x in alleles]
            alt_alleles=list(set([x for x in "".join(gts) if (x!=ref and x!=".")]))
            
            if not len(alt_alleles):
                continue
            
            alt=",".join(alt_alleles)
            
            allele_map={ref:"0", ".":"."}
            for i,a in enumerate(alt_alleles):
                allele_map[a]=str(i+1)
            
            gt_strings=[allele_map[x[0]]+"/"+allele_map[x[1]] for x in gts]
            print("\t".join([chrom, poss, idd, ref, alt, "100", ".", ".", "GT"]+gt_strings))
            
        else: 
            print(line, file=sys.stderr)
            raise Exception("Header line in unexpected place")
            


        
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

