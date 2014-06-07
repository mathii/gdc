# Convert a chromopainter (?) file to an eigenstrat file.
# The haplotype file is a file of alleles, one haplotype per columm,
# one site per row. The columns are ordered by donor first and then
# by recipient.

# usage: python chromopainter2eigenstrat.py -v vcf_file.vcf(.gz) -o out_root
# will generate out_root.[snp,ind,geno].

from __future__ import division
import sys, getopt, gdc, pdb

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "hap":None, "donor":None, "recipient":None, "out":"out"  }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h:d:r:o:", ["hap", "out", "don", "rec"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        if o in ["-h","--hap"]:         options["hap"] = a
        if o in ["-d","--don"]:         options["donor"] = a
        if o in ["-r","--rec"]:         options["recipient"] = a
        elif o in ["-o","--out"]:       options["out"] = a

    print "found options:"
    print options
    return options

################################################################################

def write_ind_file(rec, don, ind):
    """
    Make the .ind file out of the recipient and donor files. 
    """
    total_individuals=0
    for what in [don, rec]:
        for line in what:
            if line=="\n":
                next
            pop,n=line[:-1].split()
            n=int(n)
            if n%2:
                raise Exception("%s has an odd number of haplotypes"%(pop))
            for i in range(n//2):
                ind.write("%s%d\tU\t%s\n"%(pop, i, pop))
                total_individuals+=1
    return total_individuals

################################################################################

def main(options):
    hap=gdc.open2(options["hap"])
    don=gdc.open2(options["donor"])
    rec=gdc.open2(options["recipient"])
    snp, ind, geno = [open(options["out"]+x, "w") for x in [".snp", ".ind", ".geno"]]

    total_individuals=write_ind_file(rec, don, ind)

    removed={"multiallelic":0}
    for line in hap:
        bits=line[:-1].split()
        chr,pos=bits[0:2]
        gts=bits[2:]
        alleles=list(set(gts))
        if len(alleles)>2:
            removed["multiallelic"]+=1
            next
        # Setting first allele to 0 - should check this is ok.
        if len(gts)!=2*total_individuals:
            raise Exception("N.genotypes!=N.individuals T %s %s"%(chr, pos))

        snp.write("\t".join([chr+":"+pos, chr, "0.0", pos, alleles[0], alleles[1]])+"\n")
        for i in range(total_individuals):
            this_gt=sum([g==alleles[0] for g in gts[(i*2):(i*2+2)]])
            geno.write(str(this_gt))
        geno.write("\n")

    [x.close() for x in [hap, don, rec, snp, ind, geno]]

################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
