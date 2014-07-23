# Converts shapeit haps/sample files to psmc or msmc input files.
# For msmc, If you supply 2 individuals, it will use both their haplotypes.
# If you supply 4, it will use one from each. 

from __future__ import division
import sys, getopt, gdc, gzip, pdb
import numpy as np

################################################################################

def parse_options():
    """
    Options:
    shapeit: root for .haps and .sample file.
    out: root for .msmc output
    individuals: either 2 or 4 samples from the sample file. 
    """
    options ={ "shapeit":None, "out":"out", "individuals":[], "psmc":False, "msmc": False }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:i:pm", ["shapeit", "out", "individuals", "psmc", "msmc"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-s","--shapeit"]:       options["shapeit"] = a
        elif o in ["-i","--individuals"]: options["individuals"] = a.split(",")
        elif o in ["-o","--out"]:         options["out"] = a
        elif o in ["-p","--psmc"]:        options["psmc"] = True
        elif o in ["-m","--msmc"]:        options["msmc"] = True

    print "found options:"
    print options
    return options

################################################################################

def read_sample(sample_file, individuals):
    """
    read the sample file and return the indices of the individuals (second col). 
    """
    names=[]

    file=open(sample_file)
    file.next()
    file.next()                       #First two lines are garbage.
    for i, line in enumerate(file):
        names.append(line.split()[1])

    index=[names.index(i) for i in individuals]
    return index

################################################################################

def index_to_index(index):
    """
    Depending on whether it's 1, 2 or 4, point to right columns
    """
    if len(index)==1:
        return np.array([2*index[0], 2*index[0]+1])
    elif len(index)==2:
        return np.array([2*index[0], 2*index[0]+1, 2*index[1], 2*index[1]+1])
    elif len(index)==4:
        return np.array([2*index[0], 2*index[1], 2*index[2], 2*index[3]])
    else:
        raise Exception("Must be either 1, 2 or 4 individuals")

################################################################################

def output_msmc(haps, chr, pos, alleles, options):
    """
    output a .msmc file. Assuming that there are 4 haplotypes
    """
    out=open(options["out"]+".msmc", "w")
    last_site=0
    for i in range(haps.shape[0]):
        if sum(haps[i,:]) not in [0,4] and pos[i]>last_site:
            this=alleles[i,:][haps[i,:]]
            out.write("%s\t%d\t%d\t%s\n"%(chr[i], pos[i], pos[i]-last_site, "".join(this)))
            last_site=pos[i]

    out.close()

################################################################################

def output_psmc(haps, chr, pos, options):
    """
    output a .psmc file. Assuming there are 2 or 4 haplotypes. If there
    are 4 then we use 0 and 2, assuming that (01) is one individual and (23)
    is another
    """

    used_haps=haps
    if haps.shape[1]==2:
        used_haps=haps[:,[0,2]]

    het_pos=pos[used_haps[:,0]!=used_haps[:,1]]
        
    out=open(options["out"]+".psmc", "w")

    out.write(">chr"+str(chr))
    
    this_block=0
    het_pos_iter=enumerate(het_pos)
    next_het_site,next_het_pos=het_pos_iter.next()
    
    while True:
        if not this_block % 60:
            out.write("\n")
        if next_het_pos > this_block*100+100:
            out.write("A")
        else:
            out.write("W")
            while next_het_pos <= this_block*100+100:   #Find the next het site after this block
                try:
                    next_het_site,next_het_pos=het_pos_iter.next()
                except StopIteration:
                    out.close()
                    return                
        this_block+=1

    
################################################################################

def main(options):
    """
    read and convert. 
    """

    if options["msmc"] and len(options["individuals"]) not in [2,4]:
        raise Exception("msmc output needs either 2 or 4 individuals")
    if options["psmc"] and len(options["individuals"]) not in [1,2]:
        raise Exception("psmc output needs either 1 or 2 individuals")
    
    index=read_sample(options["shapeit"]+".sample", options["individuals"])
    index=index_to_index(index)
    
    try:
        hap_file=gzip.open(options["shapeit"]+".haps.gz", "r")
    except IOError:
        hap_file=open(options["shapeit"]+".haps", "r")

    chr=np.genfromtxt(hap_file, dtype=str, usecols=0)
    hap_file.seek(0)
    pos=np.genfromtxt(hap_file, dtype=int, usecols=2)
    hap_file.seek(0)
    alleles=np.genfromtxt(hap_file, dtype=str, usecols=(3,4))
    hap_file.seek(0)
    haps=np.genfromtxt(hap_file, dtype=int, usecols=5+index)

    if options["msmc"]:
        output_msmc(haps, chr, pos, alleles, options)

    if options["psmc"]:
        output_psmc(haps, chr, pos, options)
            
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
