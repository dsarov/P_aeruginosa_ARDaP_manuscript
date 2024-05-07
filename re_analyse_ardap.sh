#Script to reinterrogate variant and bam files against a newly created database.

#set variables

#home to database

#ref genome

source /home/dsarovich/bin/gatk_source_card.config

#Option to repeat analysis against the CARD database
repeat_card=0

resfinder_stringency=3
#identify which files need to be created

#test for bam and bam.bai

#test for snps and indels annotated

repeat_resfinder=1

#recreate coverage and delly information

#activate environ
echo "creating environment"
source ~/miniconda3/etc/profile.d/conda.sh
conda activate ardap
echo "created"
#create all the links
echo "linking files"
reference="Pa_PAO1.fasta"
ref_base="Pa_PAO1"
reference_location="/home/dsarovich/bin/ardap/Databases/Pseudomonas_aeruginosa_pao1/Pa_PAO1.fasta"
pt_metadata="/home/dsarovich/bin/ardap/Reports/data/patientMetaData.csv"
resistance_db="/home/dsarovich/bin/ardap/Databases/Pseudomonas_aeruginosa_pao1/Pseudomonas_aeruginosa_pao1.db"
baseDir="/home/dsarovich/bin/ardap"
snpeff="Pseudomonas_aeruginosa_pao1"
#resistance_db_card="/home/dsarovich/bin/Pseudomonas_aeruginosa_pao1_CARD.db"
#card_ref="/home/dsarovich/bin/nucleotide_fasta_protein_homolog_model.fasta"
#card_ref_name="nucleotide_fasta_protein_homolog_model.fasta"
#card_ref_label="nucleotide_fasta_protein_homolog_model"

echo "linked"
seq_id=$1
echo "seq to analyse is $1"
id="${seq_id/_1_sequence.fastq.gz/}"
echo "id is $id"

#bam and bam bai test
echo "running file tests"
if [ ! -s ./Outputs/bams/"$id".dedup.bam ]; then
echo "bam file for $id not found. ardap needs to be re-run or generate bam manually. Exiting"
exit 1
else
echo "file found"
fi

if [ ! -s ./Outputs/bams/"$id".dedup.bam.bai ]; then
echo "bam file for $id not found. ardap needs to be re-run or generate bam manually. Exiting"
exit 1
else
echo "file found"
fi

if [ ! -s ./Outputs/Variants/Annotated/"$id".PASS.indels.annotated.vcf ]; then
echo "indel vcf file for $id not found. ardap needs to be re-run or generate vcf manually. Exiting"
exit 1
else
echo "file found"
fi

if [ ! -s ./Outputs/Variants/Annotated/"$id".PASS.snps.annotated.vcf ]; then
echo "snps vcf file for $id not found. ardap needs to be re-run or generate snps manually. Exiting"
exit 1
else
echo "file found"
fi

#if [ ! -e ./Outputs/CARD/"$id".CARD_primary_output.txt ]; then
#echo "CARD file for $id not found. ardap needs to be re-run or generate CARD manually. Exiting"
#exit 1
#else
#echo "file found"
#fi
#echo "finished"
echo "ceating directory structure"
#make a new temp directory for processing and link existing files
if [ ! -d ./tmp ]; then
  mkdir ./tmp
fi

sleep 1

cd ./tmp

if [ ! -d "$id" ]; then
  mkdir "$id"
  cd "$id"
else
  cd "$id"
  rm ./*
fi

sleep 1


echo "finished"
echo "linking files"
ln -s ../../Outputs/bams/"$id".dedup.bam
ln -s ../../Outputs/bams/"$id".dedup.bam.bai
ln -s ../../Outputs/Variants/Annotated/"$id".PASS.indels.annotated.vcf
ln -s ../../Outputs/Variants/Annotated/"$id".PASS.snps.annotated.vcf
ln -s ../../"$id"_1_sequence.fastq.gz
ln -s ../../"$id"_2_sequence.fastq.gz
ln -s "$reference_location"
ln -s "$pt_metadata"
echo "done"

if [ ! -d ../../Outputs/Resfinder/ ]; then
mkdir ../../Outputs/Resfinder/
fi


#run resfinder
if [ "$repeat_resfinder" == 1 ]; then
python3 ~/.nextflow/assets/dsarov/ardap/bin/resfinder/run_resfinder.py -ifq "$id"_1_sequence.fastq.gz "$id"_2_sequence.fastq.gz -acq -db_res_kma /home/dsarovich/bin/ardap/Databases/Resfinder_general -db_res /home/dsarovich/bin/ardap/Databases/Resfinder_general

mv pheno_table.txt ${id}_resfinder.txt
cp ${id}_resfinder.txt ../../Outputs/Resfinder/${id}_resfinder.txt
fi

echo "re-running variant recalculations"
recalculate_variants () {
echo "Calculating duplication and deletion events"

echo "indexing ref"
samtools faidx ${reference}
bedtools makewindows -g ${reference}.fai -w 100000000 > ${reference}.bed
echo "done"
echo "recalculating depth"
mosdepth --by ${reference}.bed output ${id}.dedup.bam
sum_depth=$(zcat output.regions.bed.gz | awk '{print $4}' | awk '{s+=$1}END{print s}')
total_chromosomes=$(zcat output.regions.bed.gz | awk '{print $4}' | wc -l)
echo "$sum_depth/$total_chromosomes" | bc > ${id}.depth.txt


echo -e "Chromosome\tStart\tEnd\tInterval" > tmp.header
zcat output.per-base.bed.gz | awk '$4 ~ /^0/ { print $1,$2,$3,$3-$2 }' > del.summary.tmp
cat tmp.header del.summary.tmp > ${id}.deletion_summary.txt
echo "depth summary created"
covdep=$(head -n 1 ${id}.depth.txt)
DUP_CUTOFF=$(echo $covdep*3 | bc)
echo "dup cutoff is $DUP_CUTOFF"

zcat output.per-base.bed.gz | awk -v DUP_CUTOFF="$DUP_CUTOFF" '$4 >= DUP_CUTOFF { print $1,$2,$3,$3-$2 }' > dup.summary.tmp

i=$(head -n1 dup.summary.tmp | awk '{ print $2 }')
k=$(tail -n1 dup.summary.tmp | awk '{ print $3 }')
chr=$(head -n1 dup.summary.tmp | awk '{ print $1 }')

awk -v i="$i" -v k="$k" -v chr="$chr" 'BEGIN {printf "chromosome " chr " start " i " "; j=i} {if (i==$2 || i==$2-1 || i==$2-2 ) {
i=$3;
}
else {
  print "end "i " interval " i-j;
  j=$2;
  i=$3;
  printf "chromosome " $1 " start "j " ";
}} END {print "end "k " interval "k-j}' < dup.summary.tmp > dup.summary.tmp1

sed -i 's/chromosome\|start \|end \|interval //g' dup.summary.tmp1
echo -e "Chromosome\tStart\tEnd\tInterval" > dup.summary.tmp.header
cat dup.summary.tmp.header dup.summary.tmp1 > ${id}.duplication_summary.txt
echo "dup summary created"
awk '{
  if (match($0,"ANN=")){print substr($0,RSTART)}
  }' ${id}.PASS.indels.annotated.vcf > indel.effects.tmp

awk -F "|" '{ print $4,$10,$11,$15 }' indel.effects.tmp | sed 's/c\.//' | sed 's/p\.//' | sed 's/n\.//'> ${id}.annotated.indel.effects

awk '{
  if (match($0,"ANN=")){print substr($0,RSTART)}
  }' ${id}.PASS.snps.annotated.vcf > snp.effects.tmp
awk -F "|" '{ print $4,$10,$11,$15 }' snp.effects.tmp | sed 's/c\.//' | sed 's/p\.//' | sed 's/n\.//' > ${id}.annotated.snp.effects

echo 'Identifying high consequence mutations'

grep 'HIGH' snp.effects.tmp  | awk -F"|" '{ print $4,$11 }' >> ${id}.Function_lost_list.txt
grep 'HIGH' indel.effects.tmp | awk -F"|" '{ print $4,$11 }' >> ${id}.Function_lost_list.txt

sed -i 's/p\.//' ${id}.Function_lost_list.txt


echo "identifying inversions"
delly call -q 5 -o ${id}.delly.bcf -g ${reference} ${id}.dedup.bam
bcftools view ${id}.delly.bcf > ${id}.delly.vcf
grep "#" ${id}.delly.vcf > delly.header
grep "<INV>" ${id}.delly.vcf > ${id}.delly.inv.vcf
grep -v "LowQual" ${id}.delly.inv.vcf > ${id}.delly.inv.vcf.tmp
cat delly.header ${id}.delly.inv.vcf.tmp > ${id}.delly.inv.vcf
echo "annotating inversions"
snpEff eff -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff ${snpeff} ${id}.delly.inv.vcf > ${id}.delly.inv.annotated.vcf
echo "filtering inversions"
if [ -s ${id}.delly.inv.vcf.tmp ]; then
  bcftools query -f '%CHROM %POS[\t%GT\t%GL]\n' ${id}.delly.inv.vcf > likelihoods.delly
  while read line; do
    echo "$line" > line.desc;
    awk '{print $4}' line.desc > geno.likelihoods;
    genotype_RR=$(awk -F"," '{print $1}' geno.likelihoods);
    genotype_RA=$(awk -F"," '{print $2}' geno.likelihoods);
    genotype_AA=$(awk -F"," '{print $3}' geno.likelihoods);
    genotype=$(awk '{print $3}' line.desc);
    if [ "$genotype" == "0/1"  ]; then
      if [ "$genotype_RR" == 0 ]; then
      echo "Genotype ignored";
      fi;
    if [ "$genotype_AA" == 0 ]; then
      echo "Genotype included";
      chromosome=$(awk '{print $1}' line.desc);
      location=$(awk '{print $2}' line.desc);
      echo -e "$chromosome\t$location" >> filtered.inversions;
    fi;
    if [ "$genotype_RA" == 0 ]; then
      alt_ref_check=0;
      alt_ref_check=$(awk -v a="$genotype_RR" -v b="$genotype_AA" 'BEGIN {if (a < b) {print "1" }}');
      if [ "$alt_ref_check" == 1 ]; then
        echo "calculating log likelihood";
        log_genotype_AA=$(awk -v a="$genotype_AA" 'BEGIN {print (10^a)}');
        log_genotype_RA=$(awk -v a="$genotype_RA" 'BEGIN {print (10^a)}');
        log_genotype_RR=$(awk -v a="$genotype_RR" 'BEGIN {print (10^a)}');
        sum_AA_RR=$(awk -v a="$log_genotype_AA" -v b="$log_genotype_RR" 'BEGIN {print (a+b)}' );
        likelihood_ratio=$(awk -v a="$log_genotype_RA" -v b="$sum_AA_RR" 'BEGIN {print (a/b)}');
        echo -e "$log_genotype_AA\t$log_genotype_RA\t$log_genotype_RR" >> likelihood.ratios.2
        echo -e "$likelihood_ratio\n" >> likelihood.ratios.2
        likelihood_ratio_test=$(awk -v a="$likelihood_ratio" 'BEGIN {if (a < 100000) {print "1" }}')
        if [ "$likelihood_ratio_test" == 1 ]; then
          echo "changing genotype to 1/1";
          chromosome=$(awk '{print $1}' line.desc);
          location=$(awk '{print $2}' line.desc);
          echo -e "$chromosome\t$location" >> filtered.inversions;
        else
          echo "Ignoring variant due to poor quality"
        fi;
	  else
	    echo "Ignoring variant due to likely reference allele"
      fi;
    fi;
  fi;
  if [ "$genotype" == "1/1"  ]; then
    echo "Genotype included";
    chromosome=$(awk '{print $1}' line.desc);
    location=$(awk '{print $2}' line.desc);
    echo -e "$chromosome\t$location" >> filtered.inversions;
  fi;
  done < likelihoods.delly
fi;

if [ -s filtered.inversions ]; then
  while read line; do grep -w "$line" ${id}.delly.inv.annotated.vcf >> ${id}.delly.inv.annotated.vcf.tmp ; done < filtered.inversions
  cat delly.header ${id}.delly.inv.annotated.vcf.tmp > ${id}.delly.inv.annotated.vcf
  awk -F"|" '/HIGH/ {f=NR} f&&NR-1==f' RS="|" ${id}.delly.inv.annotated.vcf > delly.tmp
  sed -i '/^\s*$/d' delly.tmp
  cat delly.tmp ${id}.Function_lost_list.txt > ${id}.Function_lost_list.txt.tmp
  mv ${id}.Function_lost_list.txt.tmp ${id}.Function_lost_list.txt
fi;

}

recalculate_variants
echo "finished variant recreation"

awk '!seen[$0]++' ${id}.Function_lost_list.txt > ${id}.Function_lost_list.txt.temp
mv ${id}.Function_lost_list.txt.temp ${id}.Function_lost_list.txt

echo "running SQL queries"
run_SQL () {

 #required
# ${id}.annotated.indel.effects
# ${id}.annotated.snp.effects
# ${id}.Function_lost_list.txt
# ${id}.CARD_primary_output.txt
# ${id}.duplication_summary.txt
#  ${id}.deletion_summary.txt
#  patientMetaData.csv
#  resistance_db
echo "running snps"
bash /home/dsarovich/bin/ardap/bin/SQL_queries_SNP_indel.sh ${id} ${resistance_db}
echo "running del and dup"
bash /home/dsarovich/bin/ardap/bin/SQL_queries_DelDup.sh ${id} ${resistance_db}
echo "creating reports"
bash /home/dsarovich/bin/ardap/bin/AbR_reports.sh ${id} ${resistance_db}


bash /home/dsarovich/bin/ardap/bin/Report_html.sh Pseudomonas_aeruginosa


}


run_SQL

#calculating individual drug scores
echo "Calculating cumulative scores"
awk -F "," '{ print $3 }' drug.table.txt.backup > antibiotic.loop.file

while read Antibiotic; do
  awk -F "|" -v Ab="$Antibiotic" '$4 ~ Ab { print $0 }' "$id".AbR_output.final.txt >temp.file.loop

#return just the lines containing that antibiotic
  cumulative_score=0
  while read line; do
    echo "$line" >tmp.file
    #return column matching Ab
    Ab_col=$(awk -F"|" '{print $4 }' tmp.file | sed 's/,/\n/g' | grep -n "$Antibiotic" |  awk -F: '{ print $1 }' )
    score_col=$(awk -F"|" '{print $5 }' tmp.file | awk -F"," -v col=$Ab_col '{ print $col }')
    cumulative_score=$(echo "$cumulative_score + $score_col" | bc)
	#causing a syntax error here if no antibiotic score is included.
    #for each line matching antibiotic
  done < temp.file.loop


echo "Cumulative resistance score for $Antibiotic = $cumulative_score" >> "$id".AbR_output.final.cumulative.txt

done <antibiotic.loop.file



#mv_files
echo "moving files"

if [ -s ../../Outputs/AbR_reports/"$id".AbR_output.final.txt ]; then
 date=$(date | sed 's/ /_/g')
 mv ../../Outputs/AbR_reports/"$id".AbR_output.final.txt ../../Outputs/AbR_reports/"$id".AbR_output.backup."$date"
 cp "$id".AbR_output.final.txt ../../Outputs/AbR_reports/"$id".AbR_output.final.txt
 if [ "$repeat_resfinder" == 0 ]; then
   grep "Resfinder" ../../Outputs/AbR_reports/"$id".AbR_output.backup."$date" >> ../../Outputs/AbR_reports/"$id".AbR_output.final.txt
 fi
 cp "$id"_report.html ../../Outputs/AbR_reports/"$id"_report.html
fi

if [ ! -s ../../Outputs/AbR_reports/"$id".AbR_output.final.txt ]; then
 cp "$id".AbR_output.final.txt ../../Outputs/AbR_reports/"$id".AbR_output.final.txt
 cp "$id"_report.html ../../Outputs/AbR_reports/"$id"_report.html
fi

if [ -s ../../Outputs/AbR_reports/"$id".AbR_output.final.cumulative.txt ]; then
date=$(date | sed 's/ /_/g')
 mv ../../Outputs/AbR_reports/"$id".AbR_output.final.cumulative.txt ../../Outputs/AbR_reports/"$id".AbR_output.backup.cumulative."$date"
 cp "$id".AbR_output.final.cumulative.txt ../../Outputs/AbR_reports/"$id".AbR_output.cumulative.final.txt
fi

if [ ! -s ../../Outputs/AbR_reports/"$id".AbR_output.final.cumulative.txt ]; then
 cp "$id".AbR_output.final.cumulative.txt ../../Outputs/AbR_reports/"$id".AbR_output.cumulative.final.txt
fi
