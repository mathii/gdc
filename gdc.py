# Shared functions
import gzip, pdb

################################################################################

def open2(file, mode="r"):
	"""
	Open a file, or a gzipped file if it ends in .gz
	"""
	if file[-3:]==".gz":
		return gzip.open(file, mode)
	else:
		return open(file, mode)

################################################################################

def output_msmc(haps, chr, pos, alleles, options):
    """
    output a .msmc file. Assuming that there are 4 or 8 haplotypes
    If there are 8, we use [0,2,4,6] assuming there are 4 individuals
    """
    out=open(options["out"]+".msmc", "w")
    last_site=0
    used_haps=haps
    if haps.shape[1]==8:
        used_haps=haps[:,[0,2,4,6]]

    for i in range(used_haps.shape[0]):
        if sum(used_haps[i,:]) not in [0,4] and pos[i]>last_site:
            this=alleles[i,:][used_haps[i,:]]
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
    if haps.shape[1]==4:
        used_haps=haps[:,[0,2]]

    het_pos=pos[used_haps[:,0]!=used_haps[:,1]]
        
    out=open(options["out"]+".psmc", "w")

    out.write(">chr"+str(chr))
    
    het_pos_iter=enumerate(het_pos)
    next_het_site,next_het_pos=het_pos_iter.next()
    # Start at the first het block. It's not the best thing to do,
    # but it's probably better than starting at 0.
    this_block=first_block=int(pos[0]/100)

    while True:
        if not (this_block-first_block) % 60:
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
