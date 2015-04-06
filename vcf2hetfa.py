# Convert a vcf file to hetfa format

from __future__ import division
import sys, getopt, gdc, gzip
from pyfaidx import Fasta
import pdb
################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "out":None, "ref":None, "haplotypes":False, "sample":None, "chrom":None, "refcheck" :True }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:r:s:c:h", ["vcf", "out", "ref", "sample", "chrom", "haplotypes"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-v","--vcf"]:         options["vcf"] = a
        if o in ["-r","--ref"]:         options["ref"] = a
        if o in ["-s","--sample"]:      options["sample"] = a
        if o in ["-c","--chrom"]:       options["chrom"] = a
        if o in ["-h","--haplotypes"]:  options["haplotypes"] = True
        elif o in ["-o","--out"]:       options["out"] = a
        
        if not options["haplotypes"]:
            raise Exception("Only implemented for haplotype option (-h)")
        
    print "found options:"
    print options
    return options

################################################################################

def main(options):
    """
    Iterate over the vcf and output one fasta file for each chromosome. 
    """

    ref_fa=Fasta(options["ref"])

    out0=gzip.open(options["out"]+".0.fa.gz", "w")
    out1=gzip.open(options["out"]+".1.fa.gz", "w")
    out0.write(">"+options["chrom"]+"\n")
    out1.write(">"+options["chrom"]+"\n")

    vcf=gdc.open2(options["vcf"])
    sample_idx=None
    last_pos=0
    for line in vcf:
        if line.startswith("##"):
            continue
        elif line.startswith("#"):
            bits=line.split()
            sample_idx=bits.index(options["sample"])
        else: #data line
            bits=line.split()
            gt=bits[sample_idx]
            pos=int(bits[1])
            if pos==last_pos:
                continue
            ref=bits[3]
            alt=bits[4]
            
            if len(ref)==1 and len(alt)==1 and gt in ["0|0", "1|0", "0|1", "1|1"]: #This is a phased biallelic site
                #This is the sequence from the last position to the base before the current position (note that pos is 1-based)
                ref_seq=ref_fa[options["chrom"]][last_pos:(pos-1)].seq
                if options["refcheck"] and ref_fa[options["chrom"]][pos-1].seq!=ref:
                    raise Exception("Reference mismatcah at pos "+str(pos))
                
                if gt[0]=="0":
                    out0.write(ref_seq+ref)
                elif gt[0]=="1":
                    out0.write(ref_seq+alt)
                else:
                    raise Exception("Untrapped bad genotype in haplotype 0 at pos"+str(pos))
                
                if gt[2]=="0":
                    out1.write(ref_seq+ref)
                elif gt[2]=="1":
                    out1.write(ref_seq+alt)
                else:
                    raise Exception("Untrapped bad genotype in haplotype 1 at pos"+str(pos))
                
            else:   #This is either unphased or missing or multiallelic
                out0.write("N"*(pos-last_pos))
                out1.write("N"*(pos-last_pos))
            
            last_pos=pos

    #Fill in the reference at the end and terminate with newline. 
    tail_seq=ref_fa[options["chrom"]][last_pos:].seq
    out0.write(tail_seq+"\n")
    out1.write(tail_seq+"\n")
    
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
