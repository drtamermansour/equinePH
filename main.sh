myRoot="/mnt/home/mansourt/Tamer"
cd myRoot

mkdir equinePH && cd equinePH
equPH=$(pwd)
mkdir -p $equPH/{scripts,prepdata}

script_path=$equPH/scripts
prepdata=$equPH/prepdata

horse_trans=$HOME/Tamer/horse_trans
genome_dir=$horse_trans/refGenome
genome=$genome_dir/Bowtie2Index/genome.fa


##1. download fastq data
mkdir -p $prepdata/fastq_data && cd $prepdata/fastq_data
gsutil cp "gs://2017_swidstorbucket/USDA2015 RNA sequencing/*" .

##2. check data quality 
#bash ${script_path}/run_fastqc.sh "$prepdata/fastq_data"

##3. define the the list samples for trimming
for f in $prepdata/fastq_data/*_R1_*.fastq.gz;do echo $f; done > $prepdata/fastq_data/sample_list.txt
sample_list=$prepdata/fastq_data/sample_list.txt

##4. Mild trimming with Trimmomatic using sliding window
mkdir -p $prepdata/trimmed_RNA_reads && cd $prepdata/trimmed_RNA_reads
bash ${script_path}/run_adapter_trimmer.sh $sample_list "PE" "HPC"

##5. concatenate replicates 
#for f in *-RM_*_R1_*.pe.fq;do echo ${f%-RM_*_R1_*.pe.fq};done | uniq > RM_sampleIDs
#for f in *-DZ_*_R1_*.pe.fq;do echo ${f%-DZ_*_R1_*.pe.fq};done | uniq > DZ_sampleIDs
#mkdir -p $prepdata/concat_trimmed_reads && cd $prepdata/concat_trimmed_reads
#while read id;do 
# cat $prepdata/trimmed_RNA_reads/${id}-RM_*_R1_*.pe.fq > ${id}-RM_R1_001.pe.fq;
# cat $prepdata/trimmed_RNA_reads/${id}-RM_*_R2_*.pe.fq > ${id}-RM_R2_001.pe.fq;
#done < $prepdata/trimmed_RNA_reads/RM_sampleIDs
#while read id;do  
# cat $prepdata/trimmed_RNA_reads/${id}-DZ_*_R1_*.pe.fq > ${id}-DZ_R1_001.pe.fq;
# cat $prepdata/trimmed_RNA_reads/${id}-DZ_*_R2_*.pe.fq > ${id}-DZ_R2_001.pe.fq;
#done < $prepdata/trimmed_RNA_reads/DZ_sampleIDs
######################
## prepare salmon index for differential expression analysis
mkdir $equPH/prepSalmonIndex_ncbi && cd $equPH/prepSalmonIndex_ncbi
ncbiGTF_file=$genome_dir/ref_EquCab2.0_top_level_rna.gtf
bash $script_path/run_genome_to_cdna_fasta.sh "$ncbiGTF_file" "$genome" "ncbi_transcripts.fa" "$script_path/genome_to_cdna_fasta.sh"
bash $script_path/run_salmonIndex.sh "ncbi_horse_index" "ncbi_transcripts.fa" ${script_path}/salmonIndex.sh
#awk -F"[\t\"]" 'BEGIN{OFS="\t"}{print $10,$12}' $ncbiGTF_file | sort | uniq > gene_transcript_map
ncbi_horse_index=$equPH/prepSalmonIndex_ncbi/ncbi_horse_index
#ncbi_gene_transcript_map=$equPH/prepSalmonIndex/ncbi_gene_transcript_map

## differential expression analysis
DE_dir=$equPH/DE_ncbi
mkdir -p $DE_dir && cd $DE_dir
for f in $prepdata/trimmed_RNA_reads/*_R1_001.pe.fq; do if [ -f $f ]; then
 identifier=$(basename ${f%_*_*_R1_001.pe.fq}); echo $identifier;
fi;done | uniq > $prepdata/trimmed_RNA_reads/identifiers.txt
identifiers=$prepdata/trimmed_RNA_reads/identifiers.txt

while read label;do
 identifier=$prepdata/trimmed_RNA_reads/$label;
 ls ${identifier}_*_R1_001.pe.fq ${identifier}_*_R2_001.pe.fq; #${identifier}_*_R1_001.se.fq;
 qsub -v index="$ncbi_horse_index",identifier=$identifier ${script_path}/salmonQuant_PE.sh;
 #qsub -v index="$horse_index",identifier=$identifier ${script_path}/salmonQuant_SE.sh;
done < $identifiers

#find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
for f in salmonQuant_[PS]E.e*;do grep "\[ output \]" $f; grep "Automatically detected most likely library type as" $f;done > libStats
sed -i 's/.*Automatically/Automatically/g' libStats
for f in *.quant/logs/salmon_quant.log;do grep "\[ output \]" $f;grep "Mapping rate" $f;done > salmonQuant_summary.txt
sed -i 's/.*Mapping/Mapping/g' salmonQuant_summary.txt

## create expression files suitable for edgeR input
python $script_path/gather-counts2.py

## run edgeR (repeat after exclusion of Blue samples)
module load R/3.0.1
R CMD BATCH --no-save --no-restore ${script_path}/edgeR_isoform.R
## Extracting and clustering differentially expressed genes
cd $DE_dir/effDisOnSeasonalChange
${script_path}/trinity_util/analyze_diff_expr.pl --matrix ../counts.tab -P 1e-2 -C 2 --samples samples_file
cd $DE_dir/effSummerInDis
${script_path}/trinity_util/analyze_diff_expr.pl --matrix ../counts.tab -P 1e-2 -C 2 --samples samples_file
cd $DE_dir/effSummerInHealthy
${script_path}/trinity_util/analyze_diff_expr.pl --matrix ../counts.tab -P 0.05 -C 2 --samples samples_file

cd $DE_dir/
R CMD BATCH --no-save --no-restore ${script_path}/edgeR_isoform2.R
cd $DE_dir/summerDisVsSummerHealth
${script_path}/trinity_util/analyze_diff_expr.pl --matrix ../summer_counts.tab -P 1e-2 -C 2 --samples samples_file

cd $DE_dir/
R CMD BATCH --no-save --no-restore ${script_path}/edgeR_isoform3.R
cd $DE_dir/winterDisVsSummerHealth
${script_path}/trinity_util/analyze_diff_expr.pl --matrix ../winter_counts.tab -P 1e-2 -C 2 --samples samples_file



$horse_trans/resources/NCBI/GFF/ref_EquCab2.0_top_level.gff3
## B) Extend NCBI annotation
GFF=$genome_dir/ref_EquCab2.0_top_level.gff3;
grep "ID=rna" $GFF | awk -F "\t" '{print $9}' | awk -F "[,;]" -v OFS="\t" '{ delete vars; for(i = 1; i <= NF; ++i) { n = index($i, "="); if(n) { vars[substr($i, 1, n - 1)] = substr($i, n + 1) } } Dbxref = vars["Dbxref"]; Name = vars["Name"]; gene = vars["gene"]; product = vars["product"]; } { print Dbxref,Name,gene,product }' > NCBI_TransInfo.txt;
sed -i 's/%2C/,/g' NCBI_TransInfo.txt; sed -i 's/%3B/;/g' NCBI_TransInfo.txt; sed -i 's/^GeneID://g' NCBI_TransInfo.txt; 
printf '1\ni\nGeneID\tRefSeqID\tGeneSymbol\tDescription\n.\nw\n' | ed -s NCBI_TransInfo.txt

Rscript -e 'args=(commandArgs(TRUE));data1=read.delim(args[1],quote = "");data2=read.delim(args[2],quote = "");'\
'dataMerge=merge(data1,data2,by.x="ID",by.y="RefSeqID",all.x=T,all.y=F);'\
'write.csv(dataMerge,"Summary_analysis_annotated.csv",row.names =F);' Summary_analysis.txt NCBI_TransInfo.txt



########################
## re-run differential expression with refined transcriptome
mkdir $equPH/prepSalmonIndex_davis && cd $equPH/prepSalmonIndex_davis
wget http://de.cyverse.org/dl/d/8D105FFF-C9E7-437F-8710-45FF0C98377D/refined.fasta
bash $script_path/run_salmonIndex.sh "davis_horse_index" "refined.fasta" ${script_path}/salmonIndex.sh

