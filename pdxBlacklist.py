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


#####################
### Define inputs ###
#####################

import ruffus.cmdline as cmdline
parser = cmdline.get_argparse(description='Align mouse WSG data to the human genome and call variants',
                              version = "v. 1.0",
                              ignored_args = ["--jobs", "--log_file",
                                              "--draw_graph_horizontally",
                                              "--forced_tasks"])

# add your own command line options, to the existing ones
parser.add_argument('--strain',   help='Mouse Genome Project strain (see ftp://ftp-mouse.sanger.ac.uk/REL-1604-BAM/)[string]',
                    type=str)
parser.add_argument('--cores',   default = 8, help='Number of cores to use for bcbio [integer]',
                    type=int)


# Retrieve options from command line
args = parser.parse_args()


## Available strains at Mouse Genome Project
strains = ['129P2_OlaHsd', '129S1_SvImJ', '129S5SvEvBrd', 'AKR_J', 'A_J', 'BALB_cJ', 'BTBR_T__Itpr3tf_J', 'BUB_BnJ',
           'C3H_HeH', 'C3H_HeJ', 'C57BL_10J', 'C57BL_6NJ', 'C57BR_cdJ', 'C57L_J', 'C58_J', 'CAST_EiJ', 'CBA_J',
           'DBA_1J', 'DBA_2J', 'FVB_NJ', 'I_LnJ', 'KK_HiJ', 'LEWES_EiJ', 'LP_J', 'MOLF_EiJ', 'NOD_ShiLtJ', 'NZB_B1NJ',
           'NZO_HlLtJ', 'NZW_LacJ', 'PWK_PhJ', 'RF_J', 'SEA_GnJ', 'SPRET_EiJ', 'ST_bJ', 'WSB_EiJ', 'ZALENDE_EiJ',
           'JF1_MsJ', 'LG_J', 'SJL_J', 'SM_J']

if args.strain not in strains:
    print 'Strain identifier not found...Please select one of the following:\n'
    print ', '.join(strains)
    sys.exit(1)
else:
    print('Processing strain: ' + args.strain)


########################
### Global variables ###
########################
DEBUG=True

if DEBUG:
    print '\n\n\nWARNING: DEBUGGING FLAG ACTIVE!!!\n\n\n'

WORKING_DIR = os.path.expanduser("~/pdxBlacklist")
if not os.path.exists(WORKING_DIR):
    os.makedirs(WORKING_DIR)

os.chdir(WORKING_DIR)
TMP_DIR = tempfile.gettempdir()  # identifies the current temporary directory

## System memory

CONFIG = dict(STRAIN = args.strain,
                ALIGNER = 'bwa',
                GENOME = 'hg19',
                CORES= args.cores,
                MAX_RAM='25G')

if DEBUG:
    CONFIG['STRAIN'] = 'HG04093.chrom20'


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
run_logger.info('\n\n------LOG------\n')

##################################
### Retrieve external WGS data ###
##################################
# TODO: Add file download step to Rufus
## Download the entire BAM, BAI and MD5
# TODO: locus specific queries - Use pysam/samtools



## Download section
if not DEBUG:
    ftp = ftplib.FTP('ftp-mouse.sanger.ac.uk')  # Open FTP connection
    ftp.login()                                 # Anonymous login
    ftp.cwd('REL-1604-BAM')                     # Change to relevant directory
    bam_in = args.strain + '.bam'
    bai_in = args.strain + '.bam.bai'
    # files = ftp.retrlines('LIST')               # list directory content securely, returns to stdout

else:
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
file_size = ftp.size(bam_in)   # Get size of file
ftp.sendcmd("TYPE i")    # Switch to ASCII mode

test = False
if os.path.isfile(bam_in):
    local_file_info = os.path.getsize(bam_in)
    if local_file_info == file_size:
        test = True

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


########################
### Convert to FASTQ ###
########################
from ruffus import *  ## Pipeline manager, http://www.ruffus.org.uk/index.html
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



@follows(sortBAM)
@transform(sortBAM, suffix(".sorted.bam"), ".fq")
def bam2fastq(input_file, output_files, method = ['bedtools', 'BBMap', 'samtools', 'picard', 'biobambam2', 'bamUtil'][1]):
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
        run_logger.info(' '.join(p1))
        run_logger.info(' '.join(p2))
        time_stop = timeit.default_timer()
        run_logger.info('Time(s):' + str(time_stop - time_start))

    else:
        print 'No valid method selected for FASTQ extraction.\n'
        sys.exit(1)

    ## Write a sentinel file to disk, rather than using split()
    text_file = open(output_files, "w")
    text_file.write("Placeholder file.\n")
    text_file.close()



#################
### Run BCBIO ###
#################
@follows(bam2fastq)
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
           '    remove_lcr: true',
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
    # Project directory structure
    subfolders = ['config', 'final', 'work']
    for f in subfolders:
        directory = os.path.join( WORKING_DIR, CONFIG['STRAIN'], f )
        if not os.path.exists(directory):
            os.makedirs(directory, 0755)
    ## Move/link config file
    src = os.path.join(WORKING_DIR, output_file)
    dest = os.path.join(WORKING_DIR, CONFIG['STRAIN'], 'config', output_file)
    os.symlink(src, dest )



@follows(bcbioConfig)
@transform(bcbioConfig, suffix(".yml"), ".tmp")
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

################
### Clean-up ###
################

@follows(bcbioRun)
@transform(bcbioRun, suffix(".tmp"), ".tmp1")
def filterOutputBAM(input_file, output_file):
    '''
    Filter the BCBIO generated BAM using the generated VCF

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
        if os.path.exists(directory) and not DEBUG:
            shutil.rmtree(directory)

    ## Remove files
    if not DEBUG:
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
    # pipeline_printout_graph(open("flowchart.svg", "w"), "svg", [sortBam])
