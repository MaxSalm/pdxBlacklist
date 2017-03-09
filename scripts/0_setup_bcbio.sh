## This simple script runs the automated bcbio installation script
# http://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html

## Variables
# Please ammend these as required
ROOT="/usr/local/share"
GENOME="hg19"
ALIGNER="bwa"

## Run 
wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
python bcbio_nextgen_install.py $ROOT/bcbio --genomes $GENOME --aligners $ALIGNER --distribution ubuntu --tooldir=$ROOT/bcbio
bcbio_nextgen.py upgrade -u stable 

### Ruffus
pip install ruffus
