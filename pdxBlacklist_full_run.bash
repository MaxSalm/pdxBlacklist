#!/bin/bash

## Clone the repository
# cd ~
# rm -rf pdxBlacklist
# git clone https://github.com/MaxSalm/pdxBlacklist.git
# cd pdxBlacklist
# git pull https://github.com/MaxSalm/pdxBlacklist.git

## Run pdxBlacklist for all avaialable strains (add NOD_ShiLtJ SM_J 129P2_OlaHsd)
# strains=(129S1_SvImJ 129S5SvEvBrd AKR_J A_J BALB_cJ BTBR_T__Itpr3tf_J BUB_BnJ C3H_HeH C3H_HeJ C57BL_10J C57BL_6NJ C57BR_cdJ C57L_J C58_J CAST_EiJ CBA_J DBA_1J DBA_2J FVB_NJ I_LnJ KK_HiJ LEWES_EiJ LP_J MOLF_EiJ NZB_B1NJ NZO_HlLtJ NZW_LacJ PWK_PhJ RF_J SEA_GnJ SPRET_EiJ ST_bJ WSB_EiJ ZALENDE_EiJ JF1_MsJ LG_J SJL_J)

strains=(BTBR_T__Itpr3tf_J BUB_BnJ C3H_HeH C3H_HeJ C57BL_10J C57BL_6NJ C57BR_cdJ C57L_J C58_J CAST_EiJ CBA_J DBA_1J DBA_2J FVB_NJ I_LnJ KK_HiJ LEWES_EiJ LP_J MOLF_EiJ NZB_B1NJ NZO_HlLtJ NZW_LacJ PWK_PhJ RF_J SEA_GnJ SPRET_EiJ ST_bJ WSB_EiJ ZALENDE_EiJ JF1_MsJ LG_J SJL_J)

for strain in ${strains[@]};
do 
echo "Processing strain: "$strain
python ~/pdxBlacklist/pdxBlacklist.py --strain $strain --cores 8
echo "Completed strain: "$strain
rm ~/pdxBlacklist_output/$strain.*
done

#################################
### pdxBlacklist_merge_vcf.sh ###
#################################

## Filter VCFs in protein_coding intervals 
cd pdxBlacklist_output/
mkdir ~/pdxBlacklist_filtered
to_merge=$(find . -name *-vardict.vcf.gz)
for strain in ${to_merge[@]};
do 
echo "Filtering strain: "$strain
output_name=$(basename $strain)
bcftools filter --threads 8 -O z -i 'ANN[*] ~ "protein_coding"' -o '/home/maxsalm/pdxBlacklist_filtered/'$output_name $strain
bcftools index --force --tbi '/home/maxsalm/pdxBlacklist_filtered/'$output_name 
echo "Filtered strain: "$strain
done

## Merge VCFs - Split multiallelics into seperate rows. reheader in preparation of merge. Filter to BED regions (-R) led to shortened files with UCSC BED (no idea what bug is)
# bedtools intersect -wb -a ~/refGene_ucsc_hg19.bed -b ./pdxBlacklist_output/SM_J/final/SM_J/SM_J-vardict.vcf.gz > test.vcf # Too Slow
cd ~
to_merge=$(find pdxBlacklist_filtered/ -name *-vardict.vcf.gz)
bcftools merge -m none --threads 8 -O z $to_merge > mgp_filtered.vcf.gz
bcftools index --force --tbi mgp_filtered.vcf.gz
mv mgp_* ~/pdxBlacklist_output
rm -r ~/pdxBlacklist_filtered


## Compare VCFs
A="/home/maxsalm/in_silico/single/final/pdx_synthetic_40/pdx_synthetic_40-vardict.vcf.gz"
B="/home/maxsalm/pdxBlacklist_output/A_J/final/A_J/A_J-vardict.vcf.gz" 
B="/home/maxsalm/pdxBlacklist_output/mgp_filtered.vcf.gz"
bcftools stats $A $B > vcf_comparison.txt
# Extract and write records from A shared by both A and B using exact allele match
# bcftools isec $A $B -p dir -n =2 -w 1


## Run
# nohup ./pdxBlacklist_full_run.sh > pdxBlacklist_full_run.log &
# nohup ./pdxBlacklist_merge_vcf.sh.sh > pdxBlacklist_merge_vcf.log &
## Kill 
# pstree -p 
# kill $(pstree <PID> -p -a -l | cut -d, -f2 | cut -d' ' -f1)

