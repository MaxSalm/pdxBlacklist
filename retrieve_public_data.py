__author__ = 'Max Salm'
__date__ = '23/06/2017'
__licence__ = "GPL"

### Generate an in silico PDx WGS run using SRA data

## Run example
# python retrieve_public_data.py --proportion

### Prerequisites
## Install cutadapt
# conda install -c bioconda cutadapt
## Install sratoolkit
# conda install sra-tools=2.8.1
## Install seqtk
# git clone https://github.com/lh3/seqtk.git;
# cd seqtk; make


########################
### System libraries ###
########################

import os
import sys
import argparse

USE_CMDLINE = True ## Simple toggle to simplify debugging

if USE_CMDLINE:
    parser = argparse.ArgumentParser(description='Generate an in silico PDx WGS run using SRA data',
                                     prog='retrieve data from SRA',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('proportion', help = 'Proportion of human reads in final output [0-1]', type=float)

    ## Gather arguments
    args = parser.parse_args()
    PROPORTION = args.proportion
else:
    PROPORTION = 0.5

print "\n\tProportion of human reads:" + str(PROPORTION)

########################
### Global Variables ###
########################
YIELD = 5000 # Number of reads to download from SRA

#################
### Functions ###
#################

def getSRA(srr = ["ERR194147", "ERR118255"][0], output = ["NA12878", "NODShiLtJ"][0], trim=True):
    """
    Download WGS data from SRA. Human, NA12878, https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP001960
    Download WGS data from SRA. Mouse, NOD/ShiLtJ, https://www.ncbi.nlm.nih.gov/sra/?term=ERP000927

    :return: None
    """

    f1 = output + "_1.fastq.gz"
    f2 = output + "_2.fastq.gz"
    srr_1 = srr + "_1.fastq.gz"
    srr_2 = srr + "_2.fastq.gz"

    ## Download data
    if not os.path.exists(f1) and not os.path.exists(f2):
        EXE = "fastq-dump"
        print "Processing:" + srr
        exe_h = [EXE, "--readids", "--origfmt", "-X " + str(YIELD), "--split-files", "--gzip", srr]
        print ' '.join(exe_h)
        os.system(' '.join(exe_h))

        ## Trim the additional base
        if trim:
            EXE = "cutadapt"
            exe_c = [EXE, "-u", "-1", "-o", f1, "-p", f2, srr_1, srr_2]
            os.system(' '.join(exe_c))
        else:
            os.rename(srr_1, f1)
            os.rename(srr_2, f2)

    return None



def createPDx():
    """
     Subsample and Merge the outputs, at pre-specified proportions.

     :return: None
     """

    ROOT = "pdx_synthetic_" + str(PROPORTION)

    ## Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):
    EXE = "~/seqtk/seqtk sample"
    exe_h1 = [EXE, "-s100", "NA12878_1.fastq.gz", str(int(YIELD * 2 * PROPORTION)), ">", "hsub1.fastq.gz"]
    exe_h2 = [EXE, "-s100", "NA12878_2.fastq.gz", str(int(YIELD * 2 * PROPORTION)), ">", "hsub2.fastq.gz"]
    os.system(' '.join(exe_h1))
    os.system(' '.join(exe_h2))

    exe_m1 = [EXE, "-s100", "NODShiLtJ_1.fastq.gz", str(int(YIELD * 2 * (1 - PROPORTION) )), ">", "msub1.fastq.gz"]
    exe_m2 = [EXE, "-s100", "NODShiLtJ_2.fastq.gz", str(int(YIELD * 2 * (1 - PROPORTION) )), ">", "msub2.fastq.gz"]
    os.system(' '.join(exe_m1))
    os.system(' '.join(exe_m2))

    ## Concatenate
    EXE = "cat"

    exe_c = [EXE, "hsub1.fastq.gz", "msub1.fastq.gz", ">", ROOT + "_1.fastq.gz"]
    os.system(' '.join(exe_c))

    exe_c = [EXE, "hsub1.fastq.gz", "msub1.fastq.gz", ">", ROOT + "_2.fastq.gz"]
    os.system(' '.join(exe_c))

    ## CLean-up
    os.remove("hsub1.fq.gz")
    os.remove("hsub2.fq.gz")
    os.remove("msub1.fq.gz")
    os.remove("msub2.fq.gz")

if __name__ == '__main__':
    getSRA(srr="ERR194147", output="NA12878", trim=True)
    getSRA(srr="ERR118255", output="NODShiLtJ", trim=False)
    createPDx()



