pdxBlacklist
============

Patient-derived xenograft models are extremely valuable to immunology and oncology<sup>1</sup>. However, cross-contamination between host and graft tissue is inevitable<sup>2</sup>. Accordingly, sequencing data generated from these samples likely represent a mixture of mouse and human sequences.

*pdxBlacklist* catalogs genomic variation artifacts that arise by aligning mouse-derived short-read sequencing data to the human genome. To achieve this, mouse strain-specific short-reads generated by the [Mouse Genome Project](http://www.sanger.ac.uk/science/data/mouse-genomes-project)<sup>3</sup> are aligned to the human genome, before identifying likely artifacts (e.g. SNVs and indels). The resulting variation catalog enables annotation of existing variation datasets using tools such as *Vcfanno*<sup>4</sup>.

The *pdxBlacklist* pipeline is based on the bcbio-nextgen [pipeline](http://bcbio-nextgen.readthedocs.io/en/latest/), which affords great analytic flexibility in the choice of genome builds, aligners and variant callers. Currently, the defaults are set to:

-   Genome: `hg19`
-   Aligner: `BWA`
-   Caller: `VarDict`<sup>5</sup>

To run the *pdxBlacklist* pipeline, execute the following:

    python pdxBlacklist.py --strain NOD_ShiLtJ

where `--strain` specifies the mouse strain of interest. Please select from the following strain identifiers:

`129P2_OlaHsd, 129S1_SvImJ, 129S5SvEvBrd, AKR_J, A_J, BALB_cJ, BTBR_T__Itpr3tf_J, BUB_BnJ, C3H_HeH, C3H_HeJ, C57BL_10J, C57BL_6NJ, C57BR_cdJ, C57L_J, C58_J, CAST_EiJ, CBA_J, DBA_1J, DBA_2J, FVB_NJ, I_LnJ, KK_HiJ, LEWES_EiJ, LP_J, MOLF_EiJ, NOD_ShiLtJ, NZB_B1NJ, NZO_HlLtJ, NZW_LacJ, PWK_PhJ, RF_J, SEA_GnJ, SPRET_EiJ, ST_bJ, WSB_EiJ, ZALENDE_EiJ, JF1_MsJ, LG_J, SJL_J, SM_J`

Additional options include:

-   `--help` show a help message and exit
-   `--strain` Mouse Genome Project strain (see <ftp://ftp-mouse.sanger.ac.uk/REL-1604-BAM/>) \[string\]
-   `--cores` Number of cores to use for bcbio \[integer\]
-   `--debug` Debugging mode \[boolean\]
-   `--config` Full path to filename of an optional BCBIO config file: for details, see <http://bcbio-nextgen.readthedocs.io/en/latest/contents/configuration.html>

The automatically generated default BCBIO config file is:

    details:
    - algorithm:
        aligner: bwa
        save_diskspace: true
        mark_duplicates: true
        realign: false
        recalibrate: false
        remove_lcr: false
        platform: illumina
        quality_format: standard
        variantcaller: vardict
        indelcaller: false
        phasing: false
      analysis: variant2
      description: NOD_ShiLtJ
      files:
      - /home/maxsalm/pdxBlacklist_output/NOD_ShiLtJ.1.fq
      - /home/maxsalm/pdxBlacklist_output/NOD_ShiLtJ.2.fq
      genome_build: hg19
    fc_date: 2017/5/17
    fc_name: NOD_ShiLtJ
    resources:
      tmp:
        dir: /tmp
    upload:
      dir: ../final

The pipeline output is written to disk in VCF format to a directory labelled by strain within the `pdxBlacklist_output` directory in your home directory (see BCBIO for more details on output directory structure).

Features
--------

-   A choice of mouse strains (currently 40)
-   Diverse variant analysis pipelines

Recommended Minimum System Requirements
---------------------------------------

-   Operating Systems: Linux
-   Storage: 1Tb
-   CPU: 8 cores
-   Memory: 30Gb
-   Software: Python (2.7), git

Installation
------------

The install the *pdxBlacklist* pipeline, primarily relies on [Ruffus](http://www.ruffus.org.uk/) and [bcbio-nextgen](http://bcbio-nextgen.readthedocs.io/en/latest/). To simplify package management, the `bioconda` package [manager](https://bioconda.github.io/index.html) is recommended. After installing `conda` and setting up the `bioconda` channel, please execute the following at the commandline:

    ## Setup bioconda channel
    conda config --add channels conda-forge
    conda config --add channels defaults
    conda config --add channels r
    conda config --add channels bioconda

    ## BCBIO Variables
    # Please ammend these as required
    ROOT="/usr/local/share"
    GENOME="hg19"
    ALIGNER="bwa"

    ## Install bcbio 
    wget https://raw.github.com/chapmanb/bcbio-nextgen/master/scripts/bcbio_nextgen_install.py
    python bcbio_nextgen_install.py $ROOT/bcbio --genomes $GENOME --aligners $ALIGNER --distribution ubuntu --tooldir=$ROOT/bcbio
    bcbio_nextgen.py upgrade -u stable 

    ### Install Ruffus
    conda install ruffus=2.6.3

    ### Install SRA Toolkit
    conda install -c r r-xml=3.98_1.5 
    conda install bioconductor-sradb=1.32.0
    conda install sra-tools=2.8.1

Runtime
-------

The scripts in this project have been developed and tested on an 8-core server with 30GB RAM running Ubuntu 14.04.5 LTS (GNU/Linux 3.19.0-74-generic x86\_64). A full run of the *SM\_J* strain completed in approximately two days.

Pre-processed blacklist
-----------------------

A blacklist was created by processing whole genome sequencing data generated for the Mouse Genomes Project for 18 mouse strains<sup>3</sup>. Alignment to the hg19 genome build was performed by BWA-MEM, and variants were called using VarDict: this generated a false posi-tive dataset comprising 11,119,424 SNVs, 13,77,355 indels and 2,305,881 complex variants. The merged VCF file (merged\_blacklist.vcf.gz) can be downloaded from:

<https://download.genego.com/data/merged_blacklist.vcf.gz> <https://download.genego.com/data/merged_blacklist.vcf.gz.tbi>

Contribute
----------

-   Issue Tracker: <https://github.com/MaxSalm/pdxBlacklist/issues>
-   Source Code: <https://github.com/MaxSalm/pdxBlacklist>

Support
-------

If you are having issues, please let us know via the Issue Tracker and we will endeavour to resolve them as soon as possible.

Known Bugs
----------

1.  Mouse Genome Project BAMs are retrieved via FTP using the `ftplib` python package. At the end of file transfer of large files, `ftplib` occasionally [hangs](http://stackoverflow.com/questions/19692739/python-ftplib-hangs-at-end-of-transfer). If the BAM download is complete, then this problem can be resolved simply by stopping the the *pdxBlacklist* pipeline and re-starting it.

Planned improvements
--------------------

1.  An option to use a pre-specified BAM on disk
2.  An option to process an Sequence Read Archive (SRA) submission directly

License
-------

The project is licensed under the GPL license.

Citation
--------

pdxBlacklist: Identifying artefactual variants in patient-derived xenograft samples Max Salm, Sven-Eric Schelhorn, Lee Lancashire, Thomas Grombacher bioRxiv 180752; doi: <https://doi.org/10.1101/180752> <https://www.biorxiv.org/content/early/2017/08/25/180752>

References
----------

1. Townsend, E. C. *et al.* The public repository of xenografts enables discovery and randomized phase iI-like trials in mice. *Cancer cell* **29,** 574–586 (2016).

2. Conway, T. *et al.* Xenome–a tool for classifying reads from xenograft samples. *Bioinformatics (Oxford, England)* **28,** i172–i178 (2012).

3. Keane, T. M. *et al.* Mouse genomic variation and its effect on phenotypes and gene regulation. *Nature* **477,** 289–294 (2011).

4. Pedersen, B. S., Layer, R. M. & Quinlan, A. R. Vcfanno: Fast, flexible annotation of genetic variants. *Genome biology* **17,** 118 (2016).

5. Lai, Z. *et al.* VarDict: A novel and versatile variant caller for next-generation sequencing in cancer research. *Nucleic acids research* **44,** e108 (2016).
