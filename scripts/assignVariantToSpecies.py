__author__ = 'Max Salm'
__date__ = '10/08/2017'
__licence__ = "GPL"

### Given a VCF and a matching BAM, count the read origin of each SNV

## Run example
# python assignVariantToSpecies.py -bam a.bam -vcf a.vcf -out output

### Prerequisites
# conda install -c bioconda pysam



########################
### System libraries ###
########################

import os
import sys
import argparse
import pysam


#################
### Interface ###
#################

USE_CMDLINE = True ## Simple toggle to simplify debugging

if USE_CMDLINE:
    parser = argparse.ArgumentParser(prog='assignVariantToSpecies.py',
                                     description='Given a VCF and a matching BAM, count hte read origin of each SNV',
                                     usage='python assignVariantToSpecies.py --bam a.bam --vcf a.vcf --out output',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-bam', help = 'Input BAM', required=True)
    parser.add_argument('-vcf', help = 'Input VCF', required=True)
    parser.add_argument('-out', help = 'Output prefix', required=True)

    ## Gather arguments
    args = parser.parse_args()
    BAM = str(args.bam)
    VCF = str(args.vcf)
    OUTPUT_PREFIX = str(args.out)

else:
    root = "/home/maxsalm/pdxBlacklist/scripts/pdx_bcbio_40/final/pdx_synthetic_40/"
    BAM = root + "pdx_synthetic_40-ready.bam"
    VCF = root + "pdx_synthetic_40-vardict.vcf.gz"
    OUTPUT_PREFIX = "DEBUG"


########################
### Global Variables ###
########################


##################
### Gatekeeper ###
##################
## Test that files exist
if not os.path.exists(BAM):
    print "BAM not found:" + str(BAM)


if not os.path.exists(VCF):
    print "VCF not found:" + str(VCF)



if not input_bam.is_bam:
    print "File is not a VCF:" + str(VCF)
    sys.exit(1)

if not input_vcf.is_vcf:
    print "File is not a VCF:" + str(VCF)
    sys.exit(1)

#################
### Functions ###
#################

## Support functions
def isSNP(x):
    '''
    Given a VCF record x, determine whether this is a single-nucleotide bi-alleleic variant

    :return: True/False
    '''
    if len(x.alleles) == 2 and len(x.alts) == 1 and len(x.alts[0]) == 1 and len(x.alleles[0]) == 1:
        return True
    else:
        return False


def assignReadRefAlt(x, bam = input_bam):
    '''
    Given a VCF record and a bam, return the read names that support either the REF or ALT alleles

    :param x: VCF record python object
    :param bam: BAM file handle
    :return: dictionary containing two lists of strings, ref/alt.
    '''
    iter = input_bam.pileup(x.contig, x.start,
                            x.stop)  # Iterable mpileup, covering a much wider genomic range, from start of left-most read
    output_dict = {'ref': None, 'alt': None, 'nomatch': 0, 'coverage': None}  # Initialise output object
    for pileupcolumn in iter:
        ref_reads = []
        alt_reads = []
        if pileupcolumn.pos == x.start:
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    a_read = pileupread.alignment  # pysam.AlignedSegment
                    base_at_pos = a_read.query_sequence[pileupread.query_position]
                    #print base_at_pos, x.alleles, x.info["AF"]
                    if base_at_pos == x.ref[0]:
                        ref_reads.append(a_read.query_name)  # Read matches ref at pos
                    elif base_at_pos == x.alts[0]:
                        alt_reads.append(a_read.query_name)  # Read matches alt at pos
                    else:
                        output_dict['nomatch'] += 1  ## Allele matches neither ref/alt
            output_dict['ref'] = ref_reads
            output_dict['alt'] = alt_reads
            output_dict['coverage'] = float(pileupcolumn.n)
    return output_dict


def annotateVariant(x, human_sra = "ERR194147", mouse_sra = "ERR118255", normalise = True):
    '''
    Given the output of assignReadRefAlt(), count the human/mouse reads

    :param x: assignReadRefAlt() output dictionary
    :param human_sra: Identifier string for human reads
    :param mouse_sra: Identifier string for mouse reads
    :param normalise: Divide the counts by coverage at that position [boolean]

    :return: A dictionary of counts or proportions
    '''
    ## Set params
    if normalise and x['coverage'] != None:
        coverage_at_base = float(x['coverage'])
    else:
        coverage_at_base = 1.0
    ## Output object
    output = {'hs_ref_count': 0,
              'hs_alt_count': 0,
              'mm_ref_count': 0,
              'mm_alt_count': 0}
    if x['coverage'] == None:
        pass
    else:
        hs_ref_count = [1 for name in x['ref'] if human_sra in name]
        hs_alt_count = [1 for name in x['alt'] if human_sra in name]
        mm_ref_count = [1 for name in x['ref'] if mouse_sra in name]
        mm_alt_count = [1 for name in x['alt'] if mouse_sra in name]
        output = {'hs_ref_count': sum(hs_ref_count) / coverage_at_base,
                  'hs_alt_count': sum(hs_alt_count) / coverage_at_base,
                  'mm_ref_count': sum(mm_ref_count) / coverage_at_base,
                  'mm_alt_count': sum(mm_alt_count) / coverage_at_base}
    return output



## Main algo
def countHumanMouse():
    '''
    The main algorithm.
    For every SNV with non-zero coverage, count the reads of Human or Mouse origin

    :return: A list of lists (matrix)
    '''
    output = []
    header = ['contig', 'pos', 'REF', 'ALT',
              'hs_ref_count', 'hs_alt_count',
              'mm_ref_count', 'mm_alt_count']
    output.append(header)
    tally = {'hs' : 0, 'mm' : 0}
    for variant in input_vcf.fetch():
        if not isSNP(x = variant):
            pass
        else:
            read_assignments = assignReadRefAlt(x = variant, bam = input_bam)
            if read_assignments['coverage'] != None:
                read_counts = annotateVariant(x = read_assignments)
                ## Create output list object
                output_i = [variant.contig, variant.pos,
                          variant.ref[0], variant.alts[0],
                            read_counts['hs_ref_count'],
                            read_counts['hs_alt_count'],
                            read_counts['mm_ref_count'],
                            read_counts['mm_alt_count'] ]
                output.append(output_i)
                if sum(output_i[4:6]) > 0.0:
                    tally['hs'] += 1
                if sum(output_i[6:8]) > 0.0:
                    tally['mm'] += 1
    print tally
    return output

## Write output
def writeOutput(x):
    '''
    Write matrix to disk

    :return:
    '''
    OUTPUT_FILE = os.path.join(os.path.split(BAM)[0], OUTPUT_PREFIX + "_assigned.tsv")
    file_object  = open(OUTPUT_FILE, 'w')
    for line in x:
        arow = map(str, line) + ['\n']
        arow = '\t'.join(arow)
        file_object.write(arow)

    file_object.close()
    print 'File written to disk:' + OUTPUT_FILE




###########
### Run ###
###########

if __name__ == '__main__':
    ## Open BAM & VCF
    input_vcf = pysam.VariantFile(VCF)
    input_bam = pysam.AlignmentFile(BAM, "rb")
    tmp = countHumanMouse()
    writeOutput(x = tmp)
    ## Close BAM & VCF
    input_vcf.close()
    input_bam.close()




