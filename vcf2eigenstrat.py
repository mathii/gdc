# Convert a vcf file to eigenstrat format
# removes multi-alleleic and indel sites. 

from __future__ import division
import gzip, sys, getopt, pdb

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "out":"out"  }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:", ["vcf", "out"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-v","--vcf"]:         options["vcf"] = a
        elif o in ["-o","--out"]:       options["out"] = a

    print "found options:"
    print options
    return options

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

def main(options):
    """
    Convert vcf to eigenstrat format (ind, snp and geno files)
    """
    vcf=open2(options["vcf"])
    snp, ind, geno = [open(options["out"]+x, "w") for x in [".snp", ".ind", ".geno"]]
    removed={"multiallelic":0, "indel":0}
    count=0
    
    for line in vcf:
        if line[:2]=="##":				  # Comment line
            next
        elif line[:6]=="#CHROM":			  # Header line
            inds=line.split()[9:]
            for indi in inds:
                ind.write(indi+"\tU\tPOP\n")
        else:							  # data
            bits=line.split()
            if "," in bits[4]:
                removed["indel"]+=1
                continue
            if len(bits[3])!=1 or len(bits[4])!=1:
                removed["multiallelic"]+=1
                continue
            else:
                if bits[2]==".":
                    bits[2]=bits[0]+":"+bits[1]
                snp.write("    ".join([bits[2], bits[0], "0.0", bits[1], bits[3], bits[4]])+"\n")
                geno_string=""
                for gt in bits[9:]:
                    geno_string+=decode_gt_string(gt)
                geno.write(geno_string+"\n")
                count+=1

    [f.close for f in [ind, snp, geno]]

    print "Done. Wrote "+str(count) + " sites"
    print "Excluded " + str(sum(removed.values())) + " sites"
    for key in removed:
        print "Excluded " + str(removed[key]) + " " + key
    return

################################################################################

def decode_gt_string(gt_string):
    """
    Tries to work out the genotype from a vcf genotype entry. 9 for missing
    """
    gt=gt_string.split(":")[0]
    if len(gt)==1:
        if gt=="0":                       # haploid
            return "0"
        elif gt=="1":
            return "1"
        else:
            return "9"
    elif len(gt)==3:
        if gt[0]=="0" and gt[2]=="0":
            return "0"
        if gt[0]=="0" and gt[2]=="1":
            return "1"
        if gt[0]=="1" and gt[2]=="0":
            return "1"
        if gt[0]=="1" and gt[2]=="1":
            return "2"
        else:
            return "9"

    raise Exception("Unknown genotype: "+gt)
        
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
