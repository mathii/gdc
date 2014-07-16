# Converts shapeit haps/sample files to (4-haplotype) msmc input files.
# If you supply 2 individuals, it will use both their haplotypes.
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
    options ={ "shapeit":None, "out":"out", "individuals":[] }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:i:", ["shapeit", "out", "individuals"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-s","--shapeit"]:     options["shapeit"] = a
        if o in ["-i","--individuals"]: options["individuals"] = a.split(",")
        elif o in ["-o","--out"]:       options["out"] = a

    print "found options:"
    print options
    return options

################################################################################

def read_sample(sample_file, individuals):
    """
    read the sample file and return the indices of the individuals (second col). 
    """
    if len(individuals) not in [2,4]:
        raise Exception("Must supply either 2 or 4 individuals")

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
    Depending on whether it's 2 or 4, point to right columns
    """
    if len(index)==2:
        return np.array([2*index[0], 2*index[0]+1, 2*index[1], 2*index[1]+1])
    elif len(index)==4:
        return np.array([2*index[0], 2*index[1], 2*index[2], 2*index[3]])
    else:
        raise Exception("Wrong length")
    
################################################################################

def main(options):
    """
    read and convert. 
    """

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

    out=open(options["out"]+".msmc", "w")
    last_site=0
    for i in range(haps.shape[0]):
        if sum(haps[i,:]) not in [0,4] and pos[i]>last_site:
            this=alleles[i,:][haps[i,:]]
            out.write("%s\t%d\t%d\t%s\n"%(chr[i], pos[i], pos[i]-last_site, "".join(this)))
            last_site=pos[i]

    out.close()
            
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
