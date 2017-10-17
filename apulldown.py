#Ancient DNA pulldown - pulldown from bam into pseudo-haploid data in vcf format, or readcount output
#Python3..
#Usage: python apulldown.py -b bamlist -s snps[.snp][.vcf][.vcf.gz] 

from __future__ import division, print_function
import argparse, gdc, pysam, random, pdb
from scipy.stats import binom

LOG10=2.302585092994046

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-b', '--bamlist', type=str, default="", help=
                        "list of bam files. Two columns, with sample name in first column and bam file in second column")
    parser.add_argument('-s', '--snps', type=str, default="", help=
                        "File with snps to include, either in eigenstrat .snp or .vcf format, must be sorted")
    parser.add_argument('-o', '--output', type=str, default="vcf", help=
                        "Output type, either vcf or readcounts")
    parser.add_argument('-m', '--method', type=str, default="random", help=
                        "Output type, either random or majority")
    parser.add_argument('--seed', type=int, default=12345, help=
                        "Random number generator seed")
    parser.add_argument('-e', '--error', type=float, default=0.01, help=
                        "Error rate (for likelihood calculations)")
    parser.add_argument('-d', '--damage', type=float, default=0.05, help=
                        "Damage rate (for likelihood calculations)")
    parser.add_argument('-q', '--basequal', type=int, default=30, help=
                        "Minimum base quality")
    parser.add_argument('-i', '--mapqual', type=int, default=0, help=
                        "Minimum map quality")
    parser.add_argument('--gl', dest='gl', action='store_true')
    parser.set_defaults(likelihoods=False)


    args=parser.parse_args()
    if args.output not in ["vcf", "readcounts"]:
        raise Exception("output must be either vcf or readcounts")
    if args.method not in ["random", "majority"]:
        raise Exception("output must be either random or majority")

    return parser.parse_args()

################################################################################

class bamcollection:
    """
    Hold all the bamfiles and access the pileups
    """

    def __init__(self, bamlist, options):
        self.sample_names=[]
        self.bamfiles=[]
        self.min_base_qual=options.basequal
        self.min_map_qual=options.mapqual

        bams=open(bamlist)
        for line in bams:
            bits=line.split()
            self.sample_names.append(bits[0])
            self.bamfiles.append(pysam.AlignmentFile(bits[1],"rb"))
            
    def bases(self, chrom, pos):
        """
        For each bam, return the bases at each position. 
        """
        bases=[]
        for bf in self.bamfiles:
            base=[]

            for pupcol in bf.pileup( chrom, pos-1, pos):
                for pupread in pupcol.pileups: 
                    if pupcol.pos == pos-1 and pupread.query_position !=None:
                        bq=pupread.alignment.query_qualities[pupread.query_position]>=self.min_base_qual
                        mq=pupread.alignment.mapq>=self.min_map_qual
                        if bq and mq:
                            base.append(pupread.alignment.query_sequence[pupread.query_position])
            bases.append("".join(base))
        return(bases)

################################################################################

class snplist:
    """
    An iterator over a list of snps
    """
    
    def __init__(self, file):
        self.file=file
        if file[-4:]==".snp":
            self.is_snp=True
            self.is_vcf=False
            self.data=gdc.open2(file)
        elif file[-4:]==".vcf" or file[-7:]==".vcf.gz":
            self.is_snp=False
            self.is_vcf=True
            self.data=gdc.open2(file)
            
            #Iterate through to get to the first data line
            line=next(self.data)
            if isinstance(line, bytes):
                line=line.decode()

            while line[0]!="#CHROM":
                line=next(self.data)
        else:
            print(file[-3:])
            raise Exception("Snp file must be .snp, .vcf or .vcf.gz format")
    
    def __iter__(self):
        return self
        
    def __next__(self):
        """
        Iterate over snp file:
        """
        line=next(self.data)
        if isinstance(line, bytes):
            line=line.decode()

        bits=line.split()
        if self.is_snp:
            return (bits[1], int(bits[3]), bits[0], bits[4], bits[5])
        elif self.is_vcf:
            while len(bits[3])!=1 or len(bits[4])!=1:
                line=self.data.next()
                bits=line.split()
            return (bits[0], int(bits[1]), bits[2], bits[3], bits[4])
            

################################################################################

def majority_call(counts):
    """
    (ref,alt) count - calls majority or "." if both missing
    """
    ref,alt=counts
    if ref==alt==0:
        return('.')
    elif ref==alt:
        if random.random()<0.5:
            return "0"
        else:
            return "1"
    elif ref>alt:
        return "0"
    elif alt>ref:
        return "1"
    else:
        raise Exception("Do not know how to deal with count "+str(counts))

################################################################################

def random_call(counts):
    """
    (ref,alt) count - calls alt with probability alt/(ref+alt)
    """
    ref,alt=counts
    if ref==alt==0:
        return('.')
    else:
        if random.random()< alt/(alt+ref):
            return "1"
        else:
            return "0"

################################################################################

def genotype_likelihoods(counts, ref_allele, alt_allele, options):
    """
    (ref,alt) count - return log10 scaled genotype likelihoods scaled so max likelihood is 0
    """
    ref,alt=counts

    CT=(ref_allele=="C" and alt_allele=="T") or (ref_allele=="G" and alt_allele=="A")
    TC=(ref_allele=="T" and alt_allele=="C") or (ref_allele=="A" and alt_allele=="G")

    #P(alt)|00
    if CT:
        eps=options.error+options.damage
    else:
        eps=options.error
    
    #P(alt)|01
    if CT:
        delta=0.5+options.damage
    elif TC:
        delta=0.5-options.damage
    else:
        delta=0.5

    #P(alt)|11
    if TC:
        gam=1-options.error-options.damage
    else:
        gam=1-options.error
    
    p00=binom.logpmf(alt, ref+alt, eps)/LOG10
    p01=binom.logpmf(alt, ref+alt, delta)/LOG10
    p11=binom.logpmf(alt, ref+alt, gam)/LOG10

    maxval=max([p00,p01,p11])
    return((p00-maxval, p01-maxval, p11-maxval))

################################################################################

def print_vcf_header(bams, options):
    print("##fileformat=VCFv4.3")
    print("##Source=apulldown.py")
    for k,v in vars(options).items():
        print("##Argument=<"+k+"="+str(v)+">")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    if  options.gl:
        print("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods\">")
    print("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allele depth\">")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(bams.sample_names))


################################################################################

def main(options):
    snps=snplist(options.snps)
    bams=bamcollection(options.bamlist, options)
    
    if options.output=="vcf":
        print_vcf_header(bams, options)
    
    for snp in snps:
        chrom,pos,id,ref,alt=snp
        bases=bams.bases(chrom, pos)
        counts=[(x.count(ref),x.count(alt)) for x in bases]
        if options.method=="majority":
            calls=[majority_call(x) for x in counts]
        if options.method=="random":
            calls=[random_call(x) for x in counts]
            
        if options.output=="readcounts":
            for i,sam in enumerate(bams.sample_names):
                print(" ".join([id,sam,str(counts[i][0]), str(counts[i][1])]))    
        
        if options.output=="vcf":
            if options.gl:
                likelihoods=[genotype_likelihoods(x, ref, alt, options) for x in counts]
                part1=[chrom, str(pos), id, ref, alt, ".", ".", ".", "GT:GL:AD"]
                part2=[x+"|"+x+":"+",".join([str(round(q,3)) for q in y])+":"+",".join([str(p) for p in z]) 
                       for x,y,z in zip(calls,likelihoods,counts)]
            else:
                part1=[chrom, str(pos), id, ref, alt, ".", ".", ".", "GT:AD"]
                part2=[x+"|"+x+":"+",".join([str(p) for p in z]) 
                   for x,z in zip(calls,counts)]

            print("\t".join(part1+part2))
                                       
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    random.seed(options.seed)
    main(options)

