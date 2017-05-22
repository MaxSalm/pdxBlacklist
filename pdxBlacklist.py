__author__ = 'Max Salm'


### Notes
## Ruffus will check the time-stamp order of all intermediate files, and start
## the pipeline from files whose time-stamp is later than the final output file.
## To run the pipeline afresh, either delete all intermediate and final files in the
## run dir or specify within the pipeline manually.

########################
### System libraries ###
########################
import argparse
import sys
import os
import subprocess
import tempfile
from datetime import datetime
import timeit
import ftplib
# import psutil
import ruffus
from shutil import copyfile
import multiprocessing


#####################
### Define inputs ###
#####################

import ruffus.cmdline as cmdline
parser = cmdline.get_argparse(description='Align mouse WGS data from the Mouse Genome Project to the human genome and catalog artefactual variants',
                              version = "v. 1.0",
                              ignored_args = ["--jobs", "--log_file",
                                              "--draw_graph_horizontally",
                                              "--forced_tasks"])

# add your own command line options, to the existing ones
parser.add_argument('--strain',   help='Mouse Genome Project strain (see ftp://ftp-mouse.sanger.ac.uk/REL-1604-BAM/)[string]',
                    type=str)
parser.add_argument('--cores',   default = 1, help='Number of cores to use for bcbio [integer]',
                    type=int)
parser.add_argument('--debug',  default = False, help='Debugging mode [True/False]',
                    type=bool)
parser.add_argument('--config',  default = None, help='Full path to filename of an optional BCBIO config file: for details, see http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html',
                    type=str)

# Retrieve options from command line
args = parser.parse_args()

###################
### Gatekeepers ###
###################
## Available strains at Mouse Genome Project (2017)
strains = ['129P2_OlaHsd', '129S1_SvImJ', '129S5SvEvBrd', 'AKR_J', 'A_J', 'BALB_cJ', 'BTBR_T__Itpr3tf_J', 'BUB_BnJ',
           'C3H_HeH', 'C3H_HeJ', 'C57BL_10J', 'C57BL_6NJ', 'C57BR_cdJ', 'C57L_J', 'C58_J', 'CAST_EiJ', 'CBA_J',
           'DBA_1J', 'DBA_2J', 'FVB_NJ', 'I_LnJ', 'KK_HiJ', 'LEWES_EiJ', 'LP_J', 'MOLF_EiJ', 'NOD_ShiLtJ', 'NZB_B1NJ',
           'NZO_HlLtJ', 'NZW_LacJ', 'PWK_PhJ', 'RF_J', 'SEA_GnJ', 'SPRET_EiJ', 'ST_bJ', 'WSB_EiJ', 'ZALENDE_EiJ',
           'JF1_MsJ', 'LG_J', 'SJL_J', 'SM_J']

if args.strain not in strains and not args.debug:
    print 'Strain identifier not found...Please select one of the following:\n'
    print ', '.join(strains)
    sys.exit(1)

elif args.strain not in strains and args.debug:
    print '\n\n\nWARNING: DEBUGGING FLAG ACTIVE!!!\n\n\n'
    print 'Analysing test data: HG04093.chrom20'
    args.strain = 'HG04093.chrom20'

elif args.strain in strains and args.debug:
    print '\n\n\nWARNING: DEBUGGING FLAG ACTIVE!!!\n\n\n'
    print 'Analysing test data: HG04093.chrom20'
    args.strain = 'HG04093.chrom20'

else:
    print('Processing strain: ' + args.strain)

## Cores
if args.cores > multiprocessing.cpu_count():
    print('\n\n\n Option error, --cores: More cores requested than are detected on the machine:' + str(multiprocessing.cpu_count()) + '\n\n\n')
    sys.exit(1)

## Config
if args.config is not None:
    if not os.path.isfile(args.config):
        print('\n\n\n Option error, --config: Config file not found. \n\n\n')
        sys.exit(1)

    if args.config == os.path.basename(args.config):
        # Add path of dir to filename
        args.config = os.path.realpath(args.config)



########################
### Global variables ###
########################

WORKING_DIR = os.path.expanduser("~/pdxBlacklist_output")
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)

os.chdir(WORKING_DIR)
TMP_DIR = tempfile.gettempdir()  # identifies the current temporary directory

# A list of configurations
CONFIG = dict(STRAIN = args.strain,
                ALIGNER = 'bwa',
                GENOME = 'hg19',
                CORES= args.cores,
                bam2fastq = ['bedtools', 'BBMap', 'samtools', 'picard', 'biobambam2', 'bamUtil'][0],
                MAX_RAM='25G',
              USE_FTP = False)



####################
### Task logging ###
####################
import logging
import logging.handlers
time = datetime.now()
LOG_FILENAME = TMP_DIR + '/blacklist' +  str(time.year) + str(time.month) + str(time.day) + '.log' # Log to TMP to stop lock file issues
run_logger = logging.getLogger('My_Ruffus_logger') # Set up a specific logger with our desired output level
run_logger.setLevel(logging.INFO) ## so report messages to logger[debug, info, warning, error, critical]
handler = logging.handlers.RotatingFileHandler(LOG_FILENAME, maxBytes=2000, backupCount=5) # Add the log message handler to the logger
run_logger.addHandler(handler)

run_logger.info('\n\n------DATE------\n')
run_logger.info(time) ## Add time to log-info level; level options are [debug, info, warning, error, critical].
run_logger.info('\n\n------SETUP------\n')
run_logger.info('Mouse strain = ' + CONFIG['STRAIN'])
run_logger.info('Read aligner = ' + CONFIG['ALIGNER'])
run_logger.info('Target genome = ' + CONFIG['GENOME'])
run_logger.info('BAM to FASTQ converter = ' + CONFIG['bam2fastq'])
run_logger.info('\n\n------LOG------\n')

##################################
### Retrieve external WGS data ###
##################################
from ruffus import *  ## Pipeline manager, http://www.ruffus.org.uk/index.html
# TODO: Add file download step to Rufus
## Download the entire BAM, BAI and MD5
# TODO: locus specific queries - Use pysam/samtools



## Download section
if not args.debug:
    src_dir = 'ftp://ftp-mouse.sanger.ac.uk/REL-1604-BAM/'
    ftp = ftplib.FTP('ftp-mouse.sanger.ac.uk')  # Open FTP connection
    ftp.login()                                 # Anonymous login
    ftp.cwd('REL-1604-BAM')                     # Change to relevant directory
    bam_in = args.strain + '.bam'
    bai_in = args.strain + '.bam.bai'
    # files = ftp.retrlines('LIST')               # list directory content securely, returns to stdout

else:
    src_dir = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG04093/alignment/'
    ftp = ftplib.FTP('ftp.1000genomes.ebi.ac.uk')  # Open FTP connection
    ftp.login()                                 # Anonymous login
    ftp.cwd('vol1/ftp/phase3/data/HG04093/alignment/')                     # Change to relevant directory
    bam_in = 'HG04093.chrom20.ILLUMINA.bwa.ITU.low_coverage.20121211.bam'
    bai_in = 'HG04093.chrom20.ILLUMINA.bwa.ITU.low_coverage.20121211.bam.bai'

# Check file is on server
files_on_server = ftp.nlst()
if bam_in not in files_on_server:
    print bam_in + ' not found on FTP server.\n'
    sys.exit(1)


# Check file is not already on disk
ftp.sendcmd("TYPE i")    # Switch to Binary mode
file_size = ftp.size(bam_in)   # Get size of file on FTP server
ftp.sendcmd("TYPE i")    # Switch to ASCII mode

test = False
if os.path.isfile(bam_in):
    local_file_info = os.path.getsize(bam_in)
    if local_file_info == file_size:
        test = True

## Alternatively, test that the fastqcheck file is not on disk (i.e. that the FASTQ files have been generated)
fastqcheck = bam_in.replace(".bam", ".fastqcheck")
if os.path.isfile(fastqcheck):
    print fastqcheck + ' file found: BAM will not be re-downloaded.\n'
    test = True

## Continue using FTP for download
if CONFIG['USE_FTP']:
    # Use FTP to download files
    if test:
        print 'File on disk: will not be re-downloaded.\n'
    else:
        print 'Downloading: ' + bam_in + '\n'
        ## Get BAI
        bai_file_local = open(bai_in, "wb")         # Open connection on disk
        cmd = 'RETR ' + bai_in                      # Retrieve command
        ftp.retrbinary(cmd, bai_file_local.write)   # retrieve remote file, and write to disk
        bai_file_local.close()                      # Close connection on disk
        run_logger.info('Downloaded:' + bai_in)

        ## Get BAM
        bam_file_local = open(bam_in, "wb")         # Open connection on disk
        cmd = 'RETR ' + bam_in                      # Retrieve command
        ftp.retrbinary(cmd, bam_file_local.write)   # retrieve remote file, and write to disk
        bam_file_local.close()                      # Close connection on disk
        run_logger.info('Downloaded:' + bam_in)
        print 'Download complete\n'
else:
    # Use wget to download files - this won't hang for larger files
    if test:
        print 'File on disk: will not be re-downloaded.\n'
    else:
        print 'Downloading: ' + bam_in + '\n'
        ## Get BAI
        cmd = ['wget', src_dir + bai_in]
        subprocess.call(cmd)
        run_logger.info('Downloaded:' + bai_in)

        ## Get BAM
        cmd = ['wget', src_dir + bam_in]
        subprocess.call(cmd)
        run_logger.info('Downloaded:' + bam_in)
        print 'Download complete\n'


ftp.quit() # close FTP connection


### Misc helper functions
def installBBMap():
    '''
    Install the BBMap package
    :return:
    '''
    existing_installation = os.path.join(WORKING_DIR, 'bbmap')

    if not os.path.exists(existing_installation):
        target_file = 'BBMap_36.99.tar.gz'
        cmd = ['wget', 'https://sourceforge.net/projects/bbmap/files/' + target_file]
        subprocess.call(cmd)
        cmd = ['tar', '-xvzf', target_file]
        subprocess.call(cmd)
        os.remove(target_file)
    else:
        print 'BBMap installation found.\n'

def installFastQValidator():
    '''
    Install the FatqValidator package http://genome.sph.umich.edu/wiki/FastQValidator
    :return:
    '''
    existing_installation = os.path.join(WORKING_DIR, 'fastQValidator_0.1.1')

    if not os.path.exists(existing_installation):
        target_file = 'FastQValidator.0.1.1.tgz'
        cmd = ['wget', 'http://genome.sph.umich.edu/w/images/3/39/' +  target_file]
        subprocess.call(cmd)
        cmd = ['tar', '-xvzf', target_file]
        subprocess.call(cmd)
        os.remove(target_file)

    else:
        print 'FastQValidator installation found.\n'




def BAM2FASTQ(input_file = bam_in):
    '''
    Extract paired FASTQ records from a BAM.
    :return:
    '''
    import pysam
    samfile = pysam.AlignmentFile(bam_in, "rb") # Open file handle
    # Iterate through all reads


def blankAFile(file_path):
    '''
    truncate a file to zero bytes, and preserve its original modification time
    Adapted from 'Keeping Large intermediate files' (http://www.ruffus.org.uk/faq.html)

    :param file: Input file path
    :return: None
    '''
    if os.path.exists(file_path):
        timeInfo = os.stat(file_path) # retrieve current time stamp of the file
        try:
            f = open(file_path,'w')
        except IOError:
            pass
        else:
            f.truncate(0)
            f.close()
            # change the time of the file back to what it was
            os.utime(file_path,(timeInfo.st_atime, timeInfo.st_mtime))
            print file_path + ' blanked to save disk-space.'
    else:
        print 'blankAFile: ' + file_path + ' not found.'
        sys.exit(1)


########################
### Convert to FASTQ ###
########################
@transform(bam_in, suffix(".bam"), ".sorted.bam")
def sortBAM(input_file, output_file, cores=CONFIG['CORES']):
    '''
    # Re-order and filter BAM file to prepare for FASTQ extraction
    '''
    time_start = timeit.default_timer()
    cmd = ['sambamba', 'sort',
           '--tmpdir', './',
           '-n',
           '-t', str(cores),
           '-F', '"proper_pair"',
           '-l', '1',
           '-m', CONFIG['MAX_RAM'],
           '-p',
           '-o', output_file,
           input_file]
    os.system(' '.join(cmd))
    ## Log
    run_logger.info(' '.join(cmd))
    time_stop = timeit.default_timer()
    run_logger.info('Time(s):' + str(time_stop - time_start))
    ## Test output files exists
    if os.path.isfile(output_file):
        print 'Sorted BAM created.\n'
    else:
        print 'Sorted BAM not found.\n'
        sys.exit(1)

    ## Remove the contents of input file, but keep a "ghost"
    blankAFile(input_file)


@follows(sortBAM)
@transform(sortBAM, suffix(".sorted.bam"), ".fq")
def bam2fastq(input_file, output_files, method = CONFIG['bam2fastq']):
    '''
    Conversion utility for extracting FASTQ records from sequence alignments in BAM format.
    There are a host of different methods, only a few of which are implemented here.
    E.g.
    http://genome.sph.umich.edu/wiki/BamUtil:_bam2FastQ
    https://github.com/berguner/bam2fastq
    '''

    if method == 'bedtools':
        print 'Extracting FASTQ using BEDTools.\n'
        time_start = timeit.default_timer()
        cmd = ['bedtools', 'bamtofastq',
               '-i', input_file,
               '-fq', output_files.replace('.fq', '.1.fq'),
              '-fq2', output_files.replace('.fq', '.2.fq')]
        os.system(' '.join(cmd))
        ## Log
        run_logger.info(' '.join(cmd))
        time_stop = timeit.default_timer()
        run_logger.info('Time(s):' + str(time_stop - time_start))


    elif method == "BBMap":
        # Meant to be much faster: https://www.biostars.org/p/223625/
        print 'Extracting FASTQ using BBMap.\n'
        print 'Method aborted - the Mouse Genome Project files use an old SAM specification that is not supported by BBMap. Please ammend switch in config.'
        sys.exit(1)
        time_start = timeit.default_timer()
        installBBMap() # Test BBMap installation
        reformat_script = os.path.join(WORKING_DIR, 'bbmap', 'reformat.sh')
        p1 = [reformat_script,
         'in=' + input_file,
          'out=' + 'stdout.fq',
           'primaryonly']
        p2 = [reformat_script,
              'in=' + 'stdin.fq',
              'out1=' + output_files.replace('.fq', '.1.fq'),
              'out2=' + output_files.replace('.fq', '.2.fq'),
              'interleaved',
              'addcolon']
        p1 = subprocess.Popen(p1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(p2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
        output = p2.communicate()[0]
        ## Log
        #run_logger.info(' '.join(p1))
        #run_logger.info(' '.join(p2))
        time_stop = timeit.default_timer()
        run_logger.info('Time(s):' + str(time_stop - time_start))

    elif method == "BAM2FASTQ":
        print 'Still in development!\n'
        sys.exit(1)

    else:
        print 'No valid method selected for FASTQ extraction.\n'
        sys.exit(1)

    ## Test output files exists
    f1 = output_files.replace('.fq', '.1.fq')
    f2 = output_files.replace('.fq', '.2.fq')
    if os.path.isfile(f1) and os.path.isfile(f2):
        print 'FASTQ files created.\n'
    else:
        print 'FASTQ files not found.\n'
        sys.exit(1)

    ## Write a sentinel file to disk, rather than using split()
    text_file = open(output_files, "w")
    text_file.write("Placeholder file.\n")
    text_file.close()

    ## Remove the contents of input file, but keep a "ghost"
    blankAFile(input_file)


### A FASTQ record parser by Heng Li - https://github.com/lh3/readfq
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break
#
# if __name__ == "__main__":
#     import sys
#     n, slen, qlen = 0, 0, 0
#     for name, seq, qual in readfq(sys.stdin):
#         n += 1
#         slen += len(seq)
#         qlen += qual and len(qual) or 0
# print n, '\t', slen, '\t', qlen


@follows(bam2fastq)
@transform(bam2fastq, suffix(".fq"), ".fastqcheck")
def checkPairedFastq(input_file, output_file):
    '''
    TODO: Given a pair of four line fastq files, check:
    [ ]      Check qual encoding
    [ ]      Check DNA alphabet
    [ ]     Output to tmp file
    [ ]     Rename tmp file to input files

    :return:
    '''
    print 'Checking that paired FASTQ records match.'
    f1 = os.path.join( WORKING_DIR, input_file.replace('.fq', '.1.fq') )
    f2 = os.path.join( WORKING_DIR, input_file.replace('.fq', '.2.fq') )

    ## Test files are of equal length, overkill
    if False:
        t1 = ['wc', '-l', f1]
        t1 = subprocess.check_output(t1)
        t2 = ['wc', '-l', f2]
        t2 = subprocess.check_output(t2)
        if t1.split(' ')[0] != t2.split(' ')[0]:
            print 'Error in checkPairedFastq(): Paired FASTQs are of different length.\n'
            sys.exit(1)
        else:
            print 'FASTQ files are of equal length.'

    ## Test FASTQ read identifiers match line by line
    # pattern = "'{split($0,a," + ''"/"' + "); print a[1]}'"
    cmd1 = "awk 'NR % 4 == 1' " + f1 + " | awk '{split($0,a,\"/\"); print a[1]}' > fastq.test1"
    cmd2 = "awk 'NR % 4 == 1' " + f2 + " | awk '{split($0,a,\"/\"); print a[1]}' > fastq.test2"
    cmd3 = ["diff", "-q", "fastq.test1", "fastq.test2"]
    subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)
    test = subprocess.check_output(cmd3)
    if test != '':
        print test
        print "Error in checkPairedFastq(): FASTQ read identifiers do not match line by line"
        sys.exit(1)
    else:
        print "FASTQ read identifiers match line by line."
        os.remove("fastq.test1")
        os.remove("fastq.test2")

    # k = 0
    # record_count = 0
    # # Cycle through both files simultaneously
    # with open(f1) as file_one, open(f2) as file_two:
    #     for line_one, line_two in zip(file_one, file_two):
    #         k += 1
    #         if k % 4 == 1:
    #             record_count += 1
    #             if line_one.split("/")[0] != line_two.split("/")[0]:
    #                 print("{0}\t{1}".format(line_one.strip(), line_two.strip()))
    #                 sys.exit(10000)
    # print "All FASTQ names match in paired order.\n"
    fo = open(output_file, 'w')
    fo.writelines('Test passed.\n')
    fo.close()


#################
### Run BCBIO ###
#################
@follows(checkPairedFastq)
@transform(bam2fastq, suffix(".fq"), ".yml")
def bcbioConfig(input_file, output_file):
    '''
    Write BCBIO configuration file to disk - due to some fiddly formatting, formal YAML library not currently used.
    :return:
    '''
    # import yaml
    # algorithm = dict(aligner='bwa', mark_duplicates='false', realign='false', recalibrate='false',
    #                  variantcaller='vardict')
    #
    # data = dict(
    #     A='a',
    #     algorithm = algorithm
    # )

    # Config file
    f1 = os.path.join( WORKING_DIR, input_file.replace('.fq', '.1.fq') )
    f2 = os.path.join( WORKING_DIR, input_file.replace('.fq', '.2.fq') )

    adate = '/'.join([str(time.year), str(time.month), str(time.day)])
    yml = ['details:',
           '- algorithm:',
           '    aligner: ' + CONFIG['ALIGNER'],
           '    save_diskspace: true',
           '    mark_duplicates: true',
           '    realign: false',
           '    recalibrate: false',
           '    remove_lcr: false',
           '    platform: illumina',
           '    quality_format: standard',
           '    # svcaller: [cnvkit, lumpy, delly]',
           '    variantcaller: ' + ['mutect', 'freebayes', 'vardict', 'varscan'][2],
           '    indelcaller: ' + ['scalpel', 'pindel', 'sid', 'false'][3],
           '    phasing: false',
           '  analysis: variant2',
           '  description: ' + CONFIG['STRAIN'],
           '  files:',
           '  - ' + f1,
           '  - ' + f2,
           '  genome_build: ' + CONFIG['GENOME'],
           'fc_date: ' + adate,
           'fc_name: ' + args.strain,
           'resources:',
           '  tmp:',
           '    dir: ' + TMP_DIR,
           'upload:',
           '  dir: ../final',
           '']

    fo = open(output_file, 'w')
    fo.writelines('\n'.join(yml))
    fo.close()

    # If a config file has been specified...
    if args.config is not None:
        print 'Copying pre-specified config file:' + args.config + '\n'
        copyfile(args.config, output_file)  # copy file to target

    # Project directory structure
    subfolders = ['config', 'final', 'work']
    for f in subfolders:
        directory = os.path.join( WORKING_DIR, CONFIG['STRAIN'], f )
        if not os.path.exists(directory):
            os.makedirs(directory, 0755)
    ## Move/link config file
    src = os.path.join(WORKING_DIR, output_file)
    dest = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'config', output_file)
    if os.path.isfile(dest):
        os.remove(dest)

    os.symlink(src, dest)



@follows(checkPairedFastq)
@transform(bcbioConfig, suffix(".yml"), ".bcbio.log")
def bcbioRun(input_file, output_file, cores = CONFIG['CORES']):
    '''
    Execute the BCBIO pipeline

    :param input_file:
    :param output_file:
    :return:
    '''
    time_start = timeit.default_timer()

    # CD to work dir
    work_tmp = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'work')
    os.chdir(work_tmp)
    # Run
    dest = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'config', input_file)
    cmd = ['bcbio_nextgen.py', dest,
           '-n', str(cores)]
    os.system(' '.join(cmd))
    os.chdir(WORKING_DIR)
    ## Log
    run_logger.info(' '.join(cmd))
    time_stop = timeit.default_timer()
    run_logger.info('Time(s):' + str(time_stop - time_start))
    ## Write a sentinel file to disk, rather than using split()
    text_file = open(output_file, "w")
    text_file.write("BCBIO has run.\n")
    text_file.close()
    print 'bcbio run complete.\n'
    ## Test output files exists
    the_vcf = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'final', CONFIG['STRAIN'],
                           CONFIG['STRAIN'] + '-vardict.vcf.gz')
    the_bam = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'final', CONFIG['STRAIN'], CONFIG['STRAIN'] + '-ready.bam')
    if os.path.isfile(the_vcf) and os.path.isfile(the_bam):
        print 'Re-aligned BAM and VCF created.\n'
    else:
        print 'BCBIO output not found.\n'
        sys.exit(1)

################
### Clean-up ###
################

@follows(bcbioRun)
@transform(bcbioRun, suffix(".bcbio.log"), ".filter.log")
def filterOutputBAM(input_file, output_file):
    '''
    Filter the BCBIO generated BAM using the generated VCF.
    TODO: This could be improved using https://github.com/walaj/VariantBam

    :return:
    '''

    the_vcf = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'final', CONFIG['STRAIN'], CONFIG['STRAIN'] + '-vardict.vcf.gz' )
    the_bam = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'final', CONFIG['STRAIN'], CONFIG['STRAIN'] + '-ready.bam' )

    time_start = timeit.default_timer()
    cmd = ['bedtools', 'intersect',
           '-a', the_bam,
           '-b', the_vcf,
           '>', the_bam.replace('-ready.bam', '-filtered.bam')]
    os.system(' '.join(cmd))

    cmd = ['samtools', 'index',
           the_bam.replace('-ready.bam', '-filtered.bam') ]

    os.system(' '.join(cmd))

    ## Log
    run_logger.info(' '.join(cmd))
    time_stop = timeit.default_timer()
    run_logger.info('Time(s):' + str(time_stop - time_start))
    ## Write a sentinel file to disk, rather than using split()
    text_file = open(output_file, "w")
    text_file.write("BAM filtering complete\n")
    text_file.close()
    print 'BAM filtering complete.\n'
    ## Remove BCBIO BAM
    os.remove(the_bam)
    print 'BAM removal complete.\n'


@follows(filterOutputBAM)
def tidyUp():
    '''
    Remove all tmp files from BCBIO, and original BAM download
    '''

    ## Remove directories
    import shutil
    subfolders = ['config', 'work']
    for f in subfolders:
        directory = os.path.join(WORKING_DIR, CONFIG['STRAIN'], f)
        if os.path.exists(directory) and not args.debug:
            shutil.rmtree(directory)

    ## Remove files
    if not args.debug:
        files_in_dir = os.listdir(WORKING_DIR)
        for f in files_in_dir:
            if CONFIG['STRAIN'] in f and not os.path.isdir(f):
                targ = os.path.join(WORKING_DIR, f)
                print 'Deleting intermediate file:' + targ
                os.remove(targ)

    run_logger.info('All intermediate files deleted.')

####################
### Run Pipeline ###
####################
# After all pipelined functions, run commandline options
# cmdline.run(options)

if __name__ == '__main__':
    pipeline_run(logger=run_logger)
    # pipeline_printout_graph(open("flowchart.svg", "w"), "svg")
