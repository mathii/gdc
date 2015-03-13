# Converts a 4 haplotype macs output file (from msformatter)
# to msmc input format. Optionally add exponentially distributed phasing
# switch and flip errors between samples (0,1) and (2,3) 

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
    options ={ "input":None, "out":None, "switch_rate":0,  "flip_rate":0, "chr":0, "msmc":False, "psmc":False ,"macs":False, "length":None}
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:s:f:c:l:mpa", ["ms", "out", "switch_rate", "flip_rate", "chr", "msmc", "psmc", "macs"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-i","--ms"]:       options["ms"] = a
        elif o in ["-o","--out"]:      options["out"] = a
        elif o in ["-c","--chr"]:      options["chr"] = a
        elif o in ["-l","--length"]:   options["length"] = int(a)
        elif o in ["-s","--switch"]:   options["switch_rate"] = float(a)
        elif o in ["-f","--flip"]:     options["flip_rate"] = float(a)
        elif o in ["-m","--msmc"]:     options["msmc"] = True
        elif o in ["-p","--psmc"]:     options["psmc"] = True
        elif o in ["-a","--macs"]:     options["macs"] = True

    print "found options:"
    print options

    if not options["length"] and not options["macs"]:
        raise Exception("Must specify lenght unless using macs output")
    
    return options

################################################################################

def add_phasing_errors(haps, switch_rate, flip_rate):
    """
    Add phasing errors to 2 haplotypes with swtiches and flip rates per site
    """

    npos,nhap = haps.shape
    if nhap!=2:
        raise Exception("Can only add switch errors to 2 haplotypes")

    switches=np.random.uniform(size=npos)<switch_rate
    flips=np.random.uniform(size=npos)<flip_rate
    for i in range(npos):
        if switches[i]:
            tmp=haps[i:,0].copy()
            haps[i:,0]=haps[i:,1]
            haps[i:,1]=tmp
        if flips[i]:
            tmp=haps[i,0]
            haps[i,0]=haps[i,1]
            haps[i,1]=tmp
            
    return haps

################################################################################

def read_ms(ms_file, options):
    """ 
    Read a ms file and return positions and haplotypes
    """
    
    ms=gdc.open2(ms_file)
    if options["length"]:
        length=options["length"]           

    nhap=length=None
    line=ms.next()
    if options["macs"]:
        nhap, length = [int(x) for x in line.split()[1:3]]
    else:        #For example, srcm with -SC abs
        length=options["length"]

    for line in ms:
        if line.startswith("segsites:"):
            npos=int(line.split()[1])
        elif line.startswith("positions:"):
            pos=np.array([int(length*float(p)) for p in line.split()[1:]])
            if len(pos) != npos:
                raise Exception("Number of positions does not match segsites")
            break

    haps=np.genfromtxt(ms, dtype=int, delimiter=1)   #Assume the rest of the file is the haplotypes
    haps=np.transpose(haps)

    if haps.shape[0] != npos:
        raise Exception("Number of positions doesn't match")
    if nhap and haps.shape[1] != nhap:
        raise Exception("Number of haplotypes doesn't match")

    return length, pos, haps
    

################################################################################

def main(options):
    """
    read and convert. 
    """

    length,pos,haps=read_ms(options["ms"], options)

    npos,nhaps=haps.shape
    if nhaps not in [2,4,8]:
        raise Exception("Must have 2,4 or 8 haplotypes")
    
    if options["switch_rate"]>0 or options["flip_rate"]>0:
        site_switch_rate=options["switch_rate"]*length/len(pos)
        site_flip_rate=options["flip_rate"]*length/len(pos)        
        for i in range(int(nhaps/2)):
            haps[:,[2*i,2*i+1]]=add_phasing_errors(haps[:,[2*i,2*i+1]], site_switch_rate, site_flip_rate)

    alleles=np.zeros( (npos,2), dtype=str)
    alleles[:,0]="A"
    alleles[:,1]="T"
            
    if options["msmc"]:
        chrs=np.zeros(npos, dtype=str,)
        chrs[:]=options["chr"]
        gdc.output_msmc(haps, chrs, pos, alleles, options)

    if options["psmc"]:
        gdc.output_psmc(haps, options["chr"], pos, options)   #Assumes all on the same chromosome
            
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
