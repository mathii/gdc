# Convert a vcf file to hetfa format
# To make psmc output, use the fq2psmcfa program from the utils folder of the psmc directory: 
# vcf2hetfa.py -v vcf -r ref -s sample -c chrom | fold | fq2psmcfa - > sample.psmcfa

from __future__ import division, print_function
import sys, getopt, gdc, gzip
from pyfaidx import Fasta
import pdb

HETFA_MAP={("A","A"):"A", ("C","C"):"C", ("G", "G"):"G", ("T","T"):"T", ("A", "C"):"M", ("A","G"):"R", ("A","T"):"W", ("C","G"):"S", ("C", "T"):"Y", ("G","T"):"K"}


################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "out":None, "ref":None, "haplotypes":False, 
              "sample":None, "chrom":None, "refcheck" :True, "mask":None,
              "mask_value":None, "haploid":False }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:r:s:c:m:a:h", 
                                   ["vcf", "out", "ref", "sample", "chrom", "mask", "mask_value", "haplotypes"])
    except Exception as err:
        print(str(err), file=sys.stderr)
        sys.exit()

    for o, a in opts:
        if o in ["-v","--vcf"]:         options["vcf"] = a
        elif o in ["-r","--ref"]:         options["ref"] = a
        elif o in ["-s","--sample"]:      options["sample"] = a
        elif o in ["-c","--chrom"]:       options["chrom"] = a
        elif o in ["-m","--mask"]:        options["mask"] = a
        elif o in ["-a","--mask_value"]:  options["mask_value"] = int(a)
        elif o in ["-h","--haplotypes"]:  options["haplotypes"] = True
        elif o in ["-o","--out"]:       options["out"] = a
        
    if bool(options["mask"])!=bool(options["mask_value"]):
        raise Exception("Must specify both mask and mask_value") 
        
    if options["haplotypes"] and not options["out"]:
        raise Exception("Must specify output file if using haplotypes")
        
    print("found options:", file=sys.stderr)
    print(options, file=sys.stderr) 
    return options

################################################################################

def get_ref_seq(options, ref_fa, mask, last_pos, pos):
    """
    Extract and mask the reference sequence. Replace masked regions with N
    """
    ref_seq=ref_fa[options["chrom"]][last_pos:(pos-1)].seq
    if mask:
        mask_seq=mask[options["chrom"]][last_pos:(pos-1)].seq
        ref_seq="".join([x if y!="N" and int(y) >= options["mask_value"] else "N" for x,y in zip(ref_seq, mask_seq)])
    return ref_seq

################################################################################                                                                             

def check_mask(mask, options, pos):
    """
    Check the mask for a single postion. Return True if masked and False if ok. 
    """
    masked=False
    
    if mask:
        mask_chr=mask[options["chrom"]][pos-1].seq
        if mask_chr!="N" and int(mask_chr) < options["mask_value"]:
            masked=True
    return masked
        
################################################################################

def output_fastas(options):
    """
    output two .fa files, one for each chromosome. 
    """
    ref_fa=Fasta(options["ref"])

    out0=gzip.open(options["out"]+".0.fa.gz", "w")
    out1=gzip.open(options["out"]+".1.fa.gz", "w")
    out0.write(">"+options["chrom"]+"\n")
    out1.write(">"+options["chrom"]+"\n")

    mask=None
    if options["mask"]:
        mask=Fasta(options["mask"])

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

            masked=check_mask(mask, options, pos)
            ref_seq=get_ref_seq(options, ref_fa, mask, last_pos, pos)
            
            if len(ref)==1 and len(alt)==1 and gt in ["0|0", "1|0", "0|1", "1|1"] and not masked: #This is a phased biallelic site
                #This is the sequence from the last position to the base before the current position (note that pos is 1-based)
                if options["refcheck"] and ref_fa[options["chrom"]][pos-1].seq!=ref:
                    raise Exception("Reference mismatch at pos "+str(pos))
                
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
                out0.write(ref_seq+"N")
                out1.write(ref_seq+"N")
            
            last_pos=pos

    #Fill in the reference at the end and terminate with newline. 
    tail_seq=ref_fa[options["chrom"]][last_pos:].seq
    out0.write(tail_seq+"\n")
    out1.write(tail_seq+"\n")
    out0.close()
    out1.close()

################################################################################

def output_hetfa(options):
    """
    output a single hetfa 
    """
    ref_fa=Fasta(options["ref"])
    mask=None
    if options["mask"]:
        mask=Fasta(options["mask"])
               
    out=None
    if options["out"]:    
        out=gzip.open(options["out"]+".hetfa.fa.gz", "w")
    else: 
        out=sys.stdout
        
    out.write(">"+options["chrom"]+"\n")
    
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

            masked=check_mask(mask, options, pos)
            ref_seq=get_ref_seq(options, ref_fa, mask, last_pos, pos)
            
            if len(ref)==1 and len(alt)==1 and gt[0] in ["0","1"] and gt[2] in ["0","1"] and not masked: #This is a biallelic site
                #This is the sequence from the last position to the base before the current position (note that pos is 1-based)
                if options["refcheck"] and ref_fa[options["chrom"]][pos-1].seq!=ref:
                    raise Exception("Reference mismatch at pos "+str(pos))
                
                genotype=int(gt[0])+int(gt[2])
                if genotype==0:
                    out.write(ref_seq+ref)
                elif genotype==1:
                    hetfa_code=HETFA_MAP[tuple(sorted([ref,alt]))]
                    out.write(ref_seq+hetfa_code)
                elif genotype==2:
                    out.write(ref_seq+alt)
                else:
                    raise Exception("Untrapped bad genotype in haplotype 0 at pos"+str(pos))
                                
            else:   #This is either unphased or missing or multiallelic
                out.write(ref_seq+"N")
                
            last_pos=pos

    #Fill in the reference at the end and terminate with newline. 
    tail_seq=ref_fa[options["chrom"]][last_pos:].seq
    out.write(tail_seq+"\n")
    
################################################################################


def main(options):
    """
    Iterate over the vcf and output one fasta file for each chromosome. 
    """
    if options["haploid"]:
        output_fastas(options)
    else:
        output_hetfa(options)
        
    
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
