## define pathes
mkdir -p /mnt/ls15/scratch/users/mansourt/Tamer/dogAnn/{scripts,refResources,track_hub,pubAssemblies}
work_dir=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogAnn"
script_path=$work_dir/scripts
refRes=$work_dir/refResources
track_hub=$work_dir/track_hub
pubAssemblies=$work_dir/pubAssemblies

dogSeq=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq"
genome_dir=$"/mnt/ls15/scratch/users/mansourt/Tamer/dogSeq/refGenome"
#######################################################################
cd $refRes

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

##refGenes from NCBI
wget ftp://ftp.ncbi.nih.gov/genomes/Canis_lupus_familiaris/GFF/ref_CanFam3.1_top_level.gff3.gz
gunzip ref_CanFam3.1_top_level.gff3.gz
#  recognize RefSeq transcript changed compared to genomic sequence
#cat ref_CanFam3.1_top_level.gff3 | awk -F "\t" '$3=="mRNA"' | grep "Note=" > Note.gff3 ## the file has 1654 modified mRNA (they have a Note key)
grep -v "GeneID:486591" ref_CanFam3.1_top_level.gff3 > ref_CanFam3.1_top_level_edit.gff3
$script_path/UCSC_kent_commands/gff3ToGenePred -useName ref_CanFam3.1_top_level_edit.gff3 ref_CanFam3.1_top_level.gpred
#$script_path/UCSC_kent_commands/genePredToGtf -utr file ref_CanFam3.1_top_level.gpred ref_CanFam3.1_top_level.gtf
mkdir $refRes/ncbi && cd $refRes/ncbi
## map the genomes: create liftover files to allow mapping of NCBI annotation to UCSC tracks 
bash $script_path/mapGenome.sh $genome          ## ends by creating ncbi/NCBItoUCSC_map.sorted.chain
cd $refRes/
$script_path/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level.gpred $refRes/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped.gpred unMapped -genePred
$script_path/UCSC_kent_commands/genePredToGtf file ref_CanFam3.1_top_level_mapped.gpred ref_CanFam3.1_top_level_mapped.gtf

cd $refRes
grep -v "Curated Genomic" ref_CanFam3.1_top_level_edit.gff3 > ref_CanFam3.1_top_level_edit2.gff3
$script_path/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level_edit2.gff3 $refRes/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped.gff3 unMapped2 -gff

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
#wget http://hgdownload.soe.ucsc.edu/goldenPath/canFam3/liftOver/canFam3ToCanFam2.over.chain.gz
#gunzip canFam3ToCanFam2.over.chain.gz
#$script_path/UCSC_kent_commands/liftOver blastHg18KG.psl canFam3ToCanFam2.over.chain  blastHg18KG_mapped.psl unMapped -pslT
###########################################################################################
### Initiate the basic structure for horse track hubs
UCSCgenome=canFam3
chromSizes=$genome_dir/$UCSCgenome.chrom.sizes
## fetch the UCSC database to get the chromosome sizes
bash ${script_path}/calcChromSizes.sh $UCSCgenome $chromSizes
## Create the basic directory structure of the track hubs
mkdir -p $track_hub/$UCSCgenome/{BigBed,BigPsl}
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
pubAssemblies_psl=$work_dir/pubAssemblies_psl
mkdir -p $pubAssemblies_psl/humanPtn && cd $pubAssemblies_psl/humanPtn
cp $work_dir/refResources/blastHg18KG_mapped.psl .

cat blastHg18KG_mapped.psl | $script_path/mypslToBigPsl | sort -k1,1 -k2,2n > blastHg18KG_mapped.bigPsl
#cat blastHg18KG_mapped.psl | awk -F "\t" -v OFS='\t' '{split($9,a,""); split($21,b,","); print $14,$16,$17,$10,1000,a[2],$16,$17,0,$18,$19,$21,$12,$13,a[1],$11,$20,"","",$15,$1,$2,$3,$4;}' | sort -k1,1 -k2,2n > blastHg18KG_mapped.bigPsl
#$script_path/UCSC_kent_commands/pslToBigPsl blastHg18KG_mapped.psl stdout | sort -k1,1 -k2,2n > blastHg18KG_mapped.bigPsl
$script_path/UCSC_kent_commands/bedToBigBed -type=bed12+12 -tab -as=bigPsl.as blastHg18KG_mapped.bigPsl "$genome_dir/$UCSCgenome.chrom.sizes" blastHg18KG_mapped.BigBed
identifier=$(echo "humanPtn" | sed 's/\//_/g' | sed 's/_output//g')
cp blastHg18KG_mapped.BigBed $track_hub/$UCSCgenome/BigPsl/${identifier}.BigBed

## create non-composite entry of the assembly
priority=1
t="humanPtn"
assembly="humanPtn"
filename=$(echo $assembly | sed 's/\//_/g' | sed 's/_output//g')
echo "track $t" >> $trackDb
echo "bigDataUrl BigPsl/$filename.BigBed" >> $trackDb
echo "shortLabel $t" >> $trackDb
echo "longLabel $assembly" >> $trackDb
echo "type bigPsl" >> $trackDb
echo "colorByStrand 255,0,0 0,0,255" >> $trackDb
echo "visibility dense" >> $trackDb
echo "priority $priority" >> $trackDb
echo "html $filename" >> $trackDb
echo " " >> $trackDb
#echo $assembly >> $current_libs

#############################################################
## variant annotation
#module load annovar/20140409
#annotate_variation.pl -buildver canFam3 --downdb --webfrom ucsc refGene dogdb/ ## failed

## http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html
wget https://github.com/Ensembl/ensembl-tools/archive/release/85.zip
gunzip 85.gz
cd ~/ensembl-tools-release-85/scripts/variant_effect_predictor/
perl INSTALL.pl
## install local cache using VEP-INSTALL.pl (http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#offline)
## for API installation say "Yes". but you can skip to download cache only
## for cache files choose 12 : a merged file of RefSeq and Ensembl transcripts. Remember to use --merged when running the VEP with this cache
## for FASTA files choose 9: Canis_familiaris.CanFam3.1.dna.toplevel.fa.gz. The FASTA file should be automatically detected by the VEP when using --cache or --offline. If it is not, use "--fasta /mnt/home/mansourt/.vep/canis_familiaris/81_CanFam3.1/Canis_familiaris.CanFam3.1.dna.toplevel.fa"
## for plugins choose 0 for all
#variant_effect_predictor.pl -i GenotypeGVCFs_output_max50.raw_SNPs.vcf --cache --offline -species canis_familiaris --merged

## alterantive approach by installing the cache manually
mkdir $HOME/.vep && cd $HOME/.vep
#wget ftp://ftp.ensembl.org/pub/current_variation/VEP/canis_familiaris_merged_vep_86_CanFam3.1.tar.gz
#tar xfz canis_familiaris_merged_vep_86_CanFam3.1.tar.gz 
wget ftp://ftp.ensembl.org/pub/current_variation/VEP/canis_familiaris_vep_86_CanFam3.1.tar.gz
tar xfz canis_familiaris_vep_86_CanFam3.1.tar.gz
rm canis_familiaris_vep_86_CanFam3.1.tar.gz

## create additional local databases
cd $refRes
sed 's/>chr/>/' $genome_dir/canFam3.fa > $genome_dir/canFam3_ens.fa
module load VEP/85
## Ensembl
#wget ftp://ftp.ensembl.org/pub/release-86/gff3/canis_familiaris/Canis_familiaris.CanFam3.1.86.gff3.gz
#gunzip Canis_familiaris.CanFam3.1.86.gff3.gz
#gtf2vep.pl -i Canis_familiaris.CanFam3.1.86.gff3 -f $genome_dir/canFam3_ens.fa -d 85 -s canFam.Ensemblgff --verbose &> gtf2vep_all_ensGFF.log
gtf2vep.pl -i Canis_familiaris.CanFam3.1.86.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canFam.Ensembl --verbose &> gtf2vep_all_ens.log

## NCBI
grep -v "Curated Genomic" ref_CanFam3.1_top_level_edit.gff3 > ref_CanFam3.1_top_level_edit2.gff3
$script_path/UCSC_kent_commands/liftOver ref_CanFam3.1_top_level_edit2.gff3 $refRes/ncbi/NCBItoUCSC_map.sorted.chain ref_CanFam3.1_top_level_mapped.gff3 unMapped2 -gff
sed 's/^chr//' ref_CanFam3.1_top_level_mapped.gff3 > ref_CanFam3.1_top_level_mapped_VEP.gff3
gtf2vep.pl -i ref_CanFam3.1_top_level_mapped_VEP.gff3 -f $genome_dir/canFam3_ens.fa -d 85 -s canFam.NCBIgff --verbose &> gtf2vep_all_ncbiGFF.log
#egrep "transcript_id \"NM_|transcript_id \"XM_" ref_CanFam3.1_top_level_mapped.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"NCBI\"; gene_biotype \"protein_coding\"; $
#while read line;do
# echo "$line" | awk '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_name/gene_name/';
# echo "$line" | awk '{if($3!="start_codon" && $3!="stop_codon")print;}'
#done < ref_CanFam3.1_coding_mapped_forVEP.gtf > ref_CanFam3.1_coding_mapped_forVEP2.gtf
#sed 's/^chr//' ref_CanFam3.1_coding_mapped_forVEP2.gtf > ref_CanFam3.1_coding_mapped_forVEP3.gtf
#gtf2vep.pl -i ref_CanFam3.1_coding_mapped_forVEP3.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canis_familiaris_codingNCBI --no_transcripts --verbose &> gtf2vep_ncbi.log

## Broad
#cat canis_familiaris.protein_coding.gtf | awk -F "\t" -v OFS='\t' '{ print $0" gene_source \"Broad\"; gene_biotype \"protein_coding\"; transcript_source \"Broad\"; transcript_bi$
#while read line;do
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3=="start_codon" || $3=="stop_codon")print;}' | sed 's/exon_id ".*"; gene_name/gene_name/';
# echo "$line" | awk -F "\t" -v OFS='\t' '{if($3!="start_codon" && $3!="stop_codon")print;}'
#done < canis_familiaris.protein_coding_forVEP.gtf > canis_familiaris.protein_coding_forVEP2.gtf
#sed 's/^chr//' canis_familiaris.protein_coding_forVEP2.gtf > canis_familiaris.protein_coding_forVEP3.gtf
#gtf2vep.pl -i canis_familiaris.protein_coding_forVEP3.gtf -f $genome_dir/canFam3_ens.fa -d 85 -s canis_familiaris_Broad --no_transcripts --verbose &> gtf2vep_broad.log

### Notes in VEP:
## When you use annotation file to generate local cache, then try to annotate a variant out side the annotated boundries, the variant_effect_predictor.pl will give you a warning says for example: "WARNING: Could not find cache for 3:27000001-28000000; correct assembly used?"
## If the annotation file is missing a chromosome, the variant_effect_predictor.pl will give you a warning says for example: "WARNING: Chromosome X not found in cache on line 3773"
## When using GFF to generate caches, gtf2vep.pl will ignore: exons from genes without transcripts, cDNA_match, match, region  && wil generate warnings for: exons or CDS without transcripts
## when using GTF to generate caches, gtf2vep.pl will generate warnings for: start_codon, stop_codon, five_prime_utr, three_prime_utr,
## when using GTF to generate caches, some formats induce gtf2vep.pl to connect to the Ensembl MySQL server (This is my problem with NCBI gtf and it happens with one feature in the whole Ensembl gtf
