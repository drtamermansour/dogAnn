#!/bin/sh
script_path=$(dirname "${BASH_SOURCE[0]}")
###########################################################################################
## The map.chain file has the old genome as the target and the new genome as the query.
## The old genome (target) here should be the NCBI and the new (quary) is UCSC

# assembled chromosomes are too large and require very large memory. We can edit the psl files. All required data can be found in the fasta headers of ncbi files
wget --no-directories ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000002285.3.assembly.txt
grep -v "^#" GCF_000002285.3.assembly.txt | awk '{ sub("\r$", ""); print }' > GCF_000002285.3.assembly_noCommonts.txt
cat GCF_000002285.3.assembly_noCommonts.txt | awk -F "\t" '$10 != "na"' | awk -F "\t" -v OFS='\t' '{ print $5,$7,$9,"-",$10,$9,"0",$9 }' > chromosome_map.txt ## genBankAcc refSeqAcc(i.e.ncbiName) ncbilength ucscScaf ucscName ucsclenght start end

$script_path/NCBItoUCSCmapTopsl chromosome_map.txt > NCBItoUCSC_map.psl
$script_path/UCSC_kent_commands/pslToChain NCBItoUCSC_map.psl NCBItoUCSC_map.chain ## should I use "pslToChain" or "axtChain" ??
$script_path/UCSC_kent_commands/chainSort NCBItoUCSC_map.chain NCBItoUCSC_map.sorted.chain

