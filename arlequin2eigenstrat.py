# Convert an arelequin .arp file to eigenstrat format
# Generate .arp file using fastsimcoal2 -s0 option
# Only DNA sequence suppoerted (with -s0), but
# multiple populations/chromosomes are allowed.
# Ascertained data automatically detected. 
# usage: python arelequin2eigenstrat.py -a arp_file.arp(.gz) -o out_root [other options]
# will generate out_root.[snp,ind,geno].
# By default will convert to genotype data, so use the -p option if you
# want to keep it unphased. 

from __future__ import division
import gzip, sys, getopt, gdc, pdb
import numpy as np

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "arp":None, "out":"out", "phased":False  }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "a:o:p", ["arp", "out", "phased"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-a","--arp"]:         options["arp"] = a
        if o in ["-o","--out"]:         options["out"] = a
        elif o in ["-p","--phased"]:    options["phased"] = True

    print "found options:"
    print options
    return options

################################################################################

def load_from_arp(arp):
    """
    Parse an arp file into an internal format, returning site data, which is a 
    list of lists, with the sites on each chromosome, and gt data which is a 
    dict with infomation on the samples and genotypes. 
    I'm really just guessing what the structure of the file is.
    """
    arp_file=gdc.open2(arp)

    N_chr=-1                                
    N_sites=-1
    site_data=[]
    gt_data={}
    ascertained=False                   
    for line in arp_file: 
        if line.startswith("#Number of independent chromosomes"):
            N_chr=int(line.split()[-1])
        if line.startswith("#Total number of polymorphic sites") and not ascertained:
            N_sites=int(line.split()[-1])
        if line.startswith("#ASCERTAINED DATA"):
            ascertained=True
            N_sites=-1
        if line.startswith("#Number of polym. sites meeting ascertainment criterion:"):
            N_sites=int(line.split()[-1])
        if "polymorphic positions on chromosome" in line and not ascertained:
            npos=int(line.split()[1])
            chrom=int(line.split()[-1])
            line=next(arp_file)
            site_data.append(None)
            site_data[chrom-1]=[int(x.replace(",", "")) for x in line[1:].split()]
        if line.startswith("#Ascertained polymorphic positions on chromosome") and ascertained:
            chrom=int(line.split()[-1])
            line=next(arp_file)
            npos=len(line.split())
            site_data.append(None)
            site_data[chrom-1]=[int(x.replace(",", "")) for x in line[1:].split()]
        if "SampleName=" in line:
            sname=line.split("\"")[-2]
            line=arp_file.next()
            ssize=int(line.split("=")[-1])
            gt=np.zeros((N_sites, ssize), dtype='int')
            arp_file.next()               # "SampleData= {" line
            for i in range(ssize):
                line=arp_file.next()
                gt_string=line.split()[2]
                if len(gt_string)!=N_sites:
                    raise Exception("Wrong number of sites")
                gt[:,i]=[int(int(x)>0) for x in gt_string]   # converting {1,2,3}->1
            gt_data[sname]={"size":ssize, "gt":gt}

    if not N_chr==len(site_data):
        raise Exception( "Number of chromosomes does not match site data" )
    if not N_sites==sum([len(x) for x in site_data]):
        raise Exception( "Total number of sites does not match" )

    return site_data, gt_data

################################################################################

def unphase(gt_data):
    """
    Halve the number of samples and unphase. 
    """
    for sample, data in gt_data.iteritems():
        oldrows, oldcols=data["gt"].shape
        if oldcols % 2:
            raise Exception("Can't unphase odd number of chromosomes")
        new_gt=np.zeros((oldrows, oldcols//2))
        for i in range(oldcols//2):
            new_gt[:,i]=data["gt"][:,2*i]+data["gt"][:,(2*i+1)]
        data["size"]=data["size"]//2
        data["gt"]=new_gt

################################################################################

def write_snp(site_data, options):
    """
    Write the .snp file - we only write 4 columns here, 
    because we only use snp output (fastsimcoal -s0 option). 
    Actually fastsimcoal can generate ACGT data, if we really
    wanted to. 
    """
    snp_file=open(options["out"]+".snp", "w")
    for i, poss in enumerate(site_data):
        chrom=str(i+1)
        for pos in poss:
            snp_file.write( "\t"+chrom+"_"+str(pos)+"\t"+chrom+"\t0.0\t"+str(pos)+"\n")
    snp_file.close()

################################################################################

def write_ind(gt_data, options):
    """
    Write the .ind file. 
    """
    ind_file=open(options["out"]+".ind", "w")

    for sample in sorted(gt_data):
        data=gt_data[sample]
        nsample=data["size"]
        sample_string=sample.replace(" ", "_")
        for i in range(1,nsample+1):
            ind_file.write("\t"+sample_string+"_"+str(i)+"\tU\t"+sample_string+"\n")
    ind_file.close()

################################################################################

def write_geno(gt_data, options):
    """
    Write the geno file.
    """
    geno_file=open(options["out"]+".geno", "w")

    pops=sorted(gt_data)

    joined_data=np.concatenate([gt_data[pop]["gt"] for pop in pops], axis=1)
    if options["phased"]:
        np.savetxt(geno_file, 1-joined_data, fmt="%d", delimiter="")
    else:
        np.savetxt(geno_file, 2-joined_data, fmt="%d", delimiter="")

    geno_file.close()
        
################################################################################

def main(options):
    """
    """
    site_data, gt_data=load_from_arp(options["arp"])

    if not options["phased"]:
       unphase(gt_data)
    
    write_snp(site_data, options)
    write_ind(gt_data, options)
    write_geno(gt_data, options)
       

################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
