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
    parser = argparse.ArgumentParser(prog='retrieve_public_data.py',
                                     description='Generate an in silico PDx WGS run using SRA data',
                                     usage='python retrieve_public_data.py -p 0.4 -y 1',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', help = 'Percentage of human reads in final output [0-100]', required=True)
    parser.add_argument('-y', help = 'Number of reads in final output, in Millions', required=True)

    ## Gather arguments
    args = parser.parse_args()
    PROPORTION = float(args.p) / 100.0
    FINALYIELD = int(float(args.y) * 1000000)  # Number of reads in final file

else:
    PROPORTION = 0.5
    FINALYIELD = 10000


print "\n\tProportion of human reads:" + str(PROPORTION)
print "\n\tNumber of reads:" + str(FINALYIELD)

########################
### Global Variables ###
########################



##################
### Gatekeeper ###
##################
## Test for dependencies

def testDependencies():
    '''
    Trivial function to test for system dependencies

    :return:
    '''
    ## Dependencies
    from subprocess import Popen, PIPE

    ## On PATH
    EXES = ["fastq-dump", "cutadapt"]
    for i in range(len(EXES)):
        exe_h = ["which", EXES[i]]
        process = Popen(exe_h, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if stdout == '':
            print "System dependency not found:" + EXES[i]
            sys.exit(1)

    ## Off PATH
    EXES = ["seqtk"]
    for i in range(len(EXES)):
        home = os.path.expanduser('~')
        finder = ['find', home, '-name', EXES[i]]
        process = Popen(finder, stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if stdout == '':
            print "System dependency not found:" + EXES[i]
            sys.exit(1)

testDependencies()

def removeFileExists():
    '''
    Do the target output files exist on disk?
    :return:
    '''
    ROOT = "pdx_synthetic_" + str(args.p)
    f1 = ROOT + "_1.fastq.gz"
    f2 = ROOT + "_2.fastq.gz"
    if os.path.exists(f1) or os.path.exists(f2):
        print 'Previous run being erased from disk.\n'
        os.remove(f1)
        os.remove(f2)

removeFileExists()

#################
### Functions ###
#################

def getSRA(srr, output, trim=True):
    """
    Download WGS data from SRA. Human, NA12878, https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP001960
    Download WGS data from SRA. Mouse, NOD/ShiLtJ, https://www.ncbi.nlm.nih.gov/sra/?term=ERP000927

    :srr: SRA identifier for fastq-dump
    :output: FASTQ file prefix
    :trim: Trim the fastq data?

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
        exe_h = [EXE, "--split-files", "--clip", "--gzip", srr]
        print ' '.join(exe_h)
        os.system(' '.join(exe_h))

        ## Trim the additional base
        if trim:
            EXE = "cutadapt"
            exe_c = [EXE, "-u", "-1", "-o", f1, "-p", f2, srr_1, srr_2]
            os.system(' '.join(exe_c))
            os.remove(srr_1)
            os.remove(srr_2)
        else:
            os.rename(srr_1, f1)
            os.rename(srr_2, f2)
    else:
        print "SRA target file found on disk."

    return None



def createPDx():
    """
     Subsample and Merge the outputs, at pre-specified proportions.

     :return: None
     """

    ROOT = "pdx_synthetic_" + str(args.p)

    ## Subsample 10000 read pairs from two large paired FASTQ files (remember to use the same random seed to keep pairing):
    EXE = "~/seqtk/seqtk sample"
    exe_h1 = [EXE, "-s100", "NA12878_1.fastq.gz", str(int(FINALYIELD * 2 * PROPORTION)), ">", ROOT+"hsub1.fastq"]
    exe_h2 = [EXE, "-s100", "NA12878_2.fastq.gz", str(int(FINALYIELD * 2 * PROPORTION)), ">", ROOT+"hsub2.fastq"]
    os.system(' '.join(exe_h1))
    os.system(' '.join(exe_h2))

    exe_m1 = [EXE, "-s100", "NODShiLtJ_1.fastq.gz", str(int(FINALYIELD * 2 * (1 - PROPORTION) )), ">", ROOT+"msub1.fastq"]
    exe_m2 = [EXE, "-s100", "NODShiLtJ_2.fastq.gz", str(int(FINALYIELD * 2 * (1 - PROPORTION) )), ">", ROOT+"msub2.fastq"]
    os.system(' '.join(exe_m1))
    os.system(' '.join(exe_m2))

    ## Concatenate
    EXE = "cat"
    exe_c = [EXE, ROOT+"hsub1.fastq", ROOT+"msub1.fastq", ">", ROOT + "_1.fastq"]
    os.system(' '.join(exe_c))
    exe_c = [EXE, ROOT+"hsub2.fastq", ROOT+"msub2.fastq", ">", ROOT + "_2.fastq"]
    os.system(' '.join(exe_c))

    # Compress
    EXE = "gzip"
    exe_c = [EXE, ROOT + "_1.fastq"]
    os.system(' '.join(exe_c))
    exe_c = [EXE, ROOT + "_2.fastq"]
    os.system(' '.join(exe_c))

    ## CLean-up
    os.remove(ROOT+"hsub1.fastq")
    os.remove(ROOT+"hsub2.fastq")
    os.remove(ROOT+"msub1.fastq")
    os.remove(ROOT+"msub2.fastq")

if __name__ == '__main__':
    getSRA(srr="ERR194147", output="NA12878", trim=True)
    getSRA(srr="ERR118255", output="NODShiLtJ", trim=False)
    createPDx()
    print "retrieve_public_data.py complete."



