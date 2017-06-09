#!/bin/bash

## CLone the repository
cd ~
git clone https://github.com/MaxSalm/pdxBlacklist.git


## Run pdxBlacklist for all avaialable strains
strains=(129P2_OlaHsd 129S1_SvImJ 129S5SvEvBrd AKR_J A_J BALB_cJ BTBR_T__Itpr3tf_J BUB_BnJ C3H_HeH C3H_HeJ C57BL_10J C57BL_6NJ C57BR_cdJ C57L_J C58_J CAST_EiJ CBA_J DBA_1J DBA_2J FVB_NJ I_LnJ KK_HiJ LEWES_EiJ LP_J MOLF_EiJ NZB_B1NJ NZO_HlLtJ NZW_LacJ PWK_PhJ RF_J SEA_GnJ SPRET_EiJ ST_bJ WSB_EiJ ZALENDE_EiJ JF1_MsJ LG_J SJL_J NOD_ShiLtJ SM_J)
for input in ${strains[@]};
do 
echo "Processing strain: "$input
python ~/pdxBlacklist/pdxBlacklist.py --strain $input --cores 8
done



# nohup python /home/maxsalm/pdxBlacklist/pdxBlacklist.py --strain 129P2_OlaHsd --cores 8 > 129P2_OlaHsd.log &
# nohup python ~/pdxBlacklist/pdxBlacklist.py --strain 129S1_SvImJ --cores 8 > 129P2_OlaHsd.log &