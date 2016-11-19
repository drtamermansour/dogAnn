## define pathes
mkdir -p /mnt/ls15/scratch/users/mansourt/Tamer/dogAnn/{scripts,refResources,track_hub,pubAssemblies}
work_dir=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogAnn"
script_path=$work_dir/scripts
refRes=$work_dir/refResources
track_hub=$work_dir/track_hub
pubAssemblies=$work_dir/pubAssemblies

dogSeq=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq"
genome_dir=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/refGenome"

cd $work_dir/refResources

## protein coding track (broad)
wget https://www.broadinstitute.org/ftp/pub/vgb/dog/trackHub/canFam3/annotation/canis_familiaris.protein_coding.bb
$script_path/UCSC_kent_commands/bigBedToBed canis_familiaris.protein_coding.bb canis_familiaris.protein_coding.bed_detail
cat canis_familiaris.protein_coding.bed_detail | cut -f1-12  > canis_familiaris.protein_coding.bed
$script_path/UCSC_kent_commands/bedToGenePred canis_familiaris.protein_coding.bed stdout | $script_path/UCSC_kent_commands/genePredToGtf -utr file stdin canis_familiaris.protein_coding.gtf

## refGene track
#wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/refGene.txt.gz
#ucscTable=$"refGene.txt.gz"
#zcat $ucscTable | cut -f2-16 | $script_path/UCSC_kent_commands/genePredToGtf file stdin refGene.gtf
#zcat $ucscTable | cut -f2-16 | $script_path/genePredToBed > refGene.bed
#refGenes from NCBI
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz
gunzip ref_CanFam3.1_top_level.gff3.gz
#  recognize RefSeq transcript changed compared to genomic sequence
#cat ref_CanFam3.1_top_level.gff3 | awk -F "\t" '$3=="mRNA"' | grep "Note=" > Note.gff3 ## the file has 1654 modified mRNA (they have a Note key)
grep -v "GeneID:486591" ref_CanFam3.1_top_level.gff3 > ref_CanFam3.1_top_level_edit.gff3
$script_path/UCSC_kent_commands/gff3ToGenePred -useName ref_CanFam3.1_top_level_edit.gff3 ref_CanFam3.1_top_level.gpred
#$script_path/UCSC_kent_commands/genePredToGtf -utr file ref_CanFam3.1_top_level.gpred ref_CanFam3.1_top_level.gtf
mkdir ncbi && cd ncbi
## map the genomes: create liftover files to allow mapping of NCBI annotation to UCSC tracks 
bash $script_path/mapGenome.sh $genome          ## ends by creating ncbi/NCBItoUCSC_map.sorted.chain
cd ../
$script_path/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level.gpred ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped.gpred unMapped -genePred
$script_path/UCSC_kent_commands/genePredToGtf file ref_CanFam3.1_top_level_mapped.gpred ref_CanFam3.1_top_level_mapped.gtf

## ensembl gene track
wget ftp://ftp.ensembl.org/pub/release-86/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.86.gtf.gz
gunzip Canis_familiaris.CanFam3.1.86.gtf.gz

## human protein track from CanFam2
wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam2/database/blastHg18KG.txt.gz
gunzip blastHg18KG.txt.gz
cat blastHg18KG.txt | cut -f2-22 | blastHg18KG.psl
wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam2/liftOver/canFam2ToCanFam3.over.chain.gz
gunzip canFam2ToCanFam3.over.chain.gz
$script_path/UCSC_kent_commands/liftOver blastHg18KG.psl canFam2ToCanFam3.over.chain  blastHg18KG_mapped.psl unMapped -pslT

#############################################################
## variant annotation
sed 's/>chr/>/' $genome_dir/canFam3.fa > $genome_dir/canFam3_ens.fa

#gtf2vep.pl -i Canis_familiaris.CanFam3.1.86.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canis_familiaris_Ensembl

#egrep "transcript_id \"NM_|transcript_id \"XM_" ref_CanFam3.1_top_level_mapped.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"NCBI\"; gene_biotype \"protein_coding\"; transcript_source \"NCBI\"; transcript_biotype \"protein_coding\";" }' > ref_CanFam3.1_coding_mapped_forVEP.gtf 
#while read line;do
# echo "$line" | awk '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_name/gene_name/'; 
# echo "$line" | awk '{if($3!="start_codon" && $3!="stop_codon")print;}' 
#done < ref_CanFam3.1_coding_mapped_forVEP.gtf > ref_CanFam3.1_coding_mapped_forVEP2.gtf
#sed 's/^chr//' ref_CanFam3.1_coding_mapped_forVEP2.gtf > ref_CanFam3.1_coding_mapped_forVEP3.gtf
#gtf2vep.pl -i ref_CanFam3.1_coding_mapped_forVEP3.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canis_familiaris_codingNCBI --no_transcripts --verbose &> gtf2vep_ncbi.log

#cat canis_familiaris.protein_coding.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"Broad\"; gene_biotype \"protein_coding\"; transcript_source \"Broad\"; transcript_biotype \"protein_coding\";" }' > canis_familiaris.protein_coding_forVEP.gtf
#while read line;do
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_name/gene_name/';
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3!="start_codon" && $3!="stop_codon")print;}'
#done < canis_familiaris.protein_coding_forVEP.gtf > canis_familiaris.protein_coding_forVEP2.gtf
#sed 's/^chr//' canis_familiaris.protein_coding_forVEP2.gtf > canis_familiaris.protein_coding_forVEP3.gtf
#gtf2vep.pl -i canis_familiaris.protein_coding_forVEP3.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canis_familiaris_Broad --no_transcripts --verbose &> gtf2vep_broad.log

###########################################################################################
### Initiate the basic structure for horse track hubs
UCSCgenome=canFam3
chromSizes=$genome_dir/$UCSCgenome.chrom.sizes
## fetch the UCSC database to get the chromosome sizes
bash ${script_path}/calcChromSizes.sh $UCSCgenome $chromSizes
## Create the basic directory structure of the track hubs
mkdir -p $track_hub/$UCSCgenome/BigBed
###########################################################################################
## Track for public assemblies
mkdir -p $pubAssemblies/NCBI && cd $pubAssemblies/NCBI
cp $work_dir/refResources/ref_CanFam3.1_top_level_mapped.gtf ncbiAnn.gtf

mkdir -p $pubAssemblies/Ensembl86 && cd $pubAssemblies/Ensembl86
grep -v "^#" $work_dir/refResources/Canis_familiaris.CanFam3.1.86.gtf | grep -v "^MT" | awk -F "\t" -v OFS='\t' '{ if(length($1)>2){split($1,a,"."); print "chrUn_"a[1],$2,$3,$4,$5,$6,$7,$8,$9;} else {print "chr"$0;} }' > ens86Ann.gtf

cp $work_dir/refResources/Canis_familiaris.CanFam3.1.86.gtf ens86Ann.gtf

## create list of public assemblies
rm -f $pubAssemblies/public_assemblies.txt
for tissue in $pubAssemblies/*; do
  echo "$pubAssemblies" "${tissue#$pubAssemblies/}" >> $pubAssemblies/public_assemblies.txt;
done

####################
## convert the gtf files into BigBed files & copy the BigBed files to the track hub directory
update=0    ## 0 means do not update Bigbed files & 1 means update
rm -f $work_dir/public_assemblies.txt
while read ass_path assembly; do
  echo $assembly
  cd $ass_path/$assembly
  if [[ ! -f "*.BigBed" || "$update" -eq 1 ]];then
    targetAss=$(ls *.gtf)
    if [ -f "$targetAss" ];then
      bash $script_path/gtfToBigBed.sh "$targetAss" "$genome_dir/$UCSCgenome.chrom.sizes" "$script_path"
    else echo "can not find target assembly"; break;fi
    if [ -f ${targetAss%.gtf}.BigBed ];then
      identifier=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
      cp ${targetAss%.gtf}.BigBed $track_hub/$UCSCgenome/BigBed/${identifier}.BigBed
    else echo "could not make merged.BigBed file"; break; fi
  fi
  echo $ass_path/$assembly >> $work_dir/public_assemblies.txt;
done < $pubAssemblies/public_assemblies.txt

## initiate a given track hub
hub_name=$"Davis_dogAnn"
shortlabel=$"Davis_dogAnn"
longlabel=$"Annotations of canFam3 assembly"
email=$"drtamermansour@gmail.com"
cd $track_hub
bash $script_path/create_trackHub.sh "$UCSCgenome" "$hub_name" "$shortlabel" "$longlabel" "$email"

## edit the trackDb
current_libs=$track_hub/current_libs_$shortlabel
current_tissues=$track_hub/current_tiss_$shortlabel
trackDb=$track_hub/$UCSCgenome/trackDb_$shortlabel.txt
lib_assemblies=$pubAssemblies/public_assemblies.txt
tiss_assemblies=$horse_trans/emptyTemp.txt
bash $script_path/edit_trackDb.sh $current_libs $current_tissues $trackDb $lib_assemblies $tiss_assemblies

######################
mkdir -p $pubAssemblies_psl/humanPtn && cd $pubAssemblies_psl/humanPtn
cp $work_dir/refResources/blastHg18KG_mapped.psl .

