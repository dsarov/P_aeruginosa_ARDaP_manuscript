#Script that takes the output of two strains from ARDaP and reports all differences
#Usage information
#Must be rrrn from the ARDaP analysis directory as all paths are relative
#var 1 = strain 1
#var 2 = strain 2

strain_1=$1
strain_2=$2

snp_dir="$PWD/Outputs/Variants/Annotated/"
card_dir="$PWD/Outputs/CARD/"
bams="$PWD/Outputs/bams/"
compare_dir="${strain_1}_${strain_2}"
snps_1="${snp_dir}/${strain_1}.PASS.snps.annotated.vcf"
snps_2="${snp_dir}/${strain_2}.PASS.snps.annotated.vcf"
indels_1="${snp_dir}/${strain_1}.PASS.indels.annotated.vcf"
indels_2="${snp_dir}/${strain_2}.PASS.indels.annotated.vcf"
card_1="$PWD/Outputs/CARD/${strain_1}.card.unfiltered.bedcov"
card_2="$PWD/Outputs/CARD/${strain_2}.card.unfiltered.bedcov"
card_aro="/home/dsarovich/bin/ardap/Databases/CARD/aro_index.tsv"


if [ ! -d "$compare_dir" ]; then
  mkdir "$compare_dir"
  cd "$compare_dir"
else
  cd "$compare_dir"  
fi

#SNPS and indels
#find snps in strain 1 and not strain 2
while read line; do
echo "$line" > tmp.file
  chromosome=$(awk '{print $1}' tmp.file)
  snp_loc=$(awk '{print $2}' tmp.file)
  awk -v v1="$chromosome" -v v2="$snp_loc" '$0 ~ v1"\t"v2"\t" {print $0; exit 1}' "$snps_2"  &> /dev/null
  status=$? 
  if [ "$status" == 0 ]; then #SNP not found
    #echo "snp not found"
	#echo "$line"
    echo "$line" >> ${strain_1}_snps.txt
  fi
done < <(grep -v '#' "$snps_1")

#find snps in strain 2 and not strain 1
while read line; do
echo "$line" > tmp.file
  chromosome=$(awk '{print $1}' tmp.file)
  snp_loc=$(awk '{print $2}' tmp.file)
  awk -v v1="$chromosome" -v v2="$snp_loc" '$0 ~ v1"\t"v2"\t" {print $0; exit 1}' "$snps_1"  &> /dev/null
  status=$? 
  if [ "$status" == 0 ]; then #SNP not found
    #echo "snp not found"
	#echo "$line"
    echo "$line" >> ${strain_2}_snps.txt
  fi
done < <(grep -v '#' "$snps_2")

#find indels in strain 1 and not strain 2
while read line; do
echo "$line" > tmp.file
  chromosome=$(awk '{print $1}' tmp.file)
  snp_loc=$(awk '{print $2}' tmp.file)
  awk -v v1="$chromosome" -v v2="$snp_loc" '$0 ~ v1"\t"v2"\t" {print $0; exit 1}' "$indels_2"  &> /dev/null
  status=$? 
  if [ "$status" == 1 ]; then #SNP not found
    #echo "snp not found"
	#echo "$line"
    echo "$line" >> ${strain_1}_indels.txt
  fi
done < <(grep -v '#' "$indels_1")

#find indels in strain 1 and not strain 2
while read line; do
echo "$line" > tmp.file
  chromosome=$(awk '{print $1}' tmp.file)
  snp_loc=$(awk '{print $2}' tmp.file)
  awk -v v1="$chromosome" -v v2="$snp_loc" '$0 ~ v1"\t"v2"\t" {print $0; exit 1}' "$indels_1"  &> /dev/null
  status=$? 
  if [ "$status" == 1 ]; then #SNP not found
    #echo "snp not found"
	#echo "$line"
    echo "$line" >> ${strain_2}_indels.txt
  fi
done < <(grep -v '#' "$indels_2")



#compare CARD
#card_dir="$PWD/Outputs/CARD/"
#refilter and annotate CARD hits that pass the filter
#card_1="$PWD/Outputs/CARD/${strain_1}.card.unfiltered.bedcov"

#refilter and annotate CARD
#strain 1
awk '{ if ( $7 > 0.97) print $0 }' "$card_1" > "$strain_1"_98p.bedcov
awk ' { print $1 } ' "$strain_1"_98p.bedcov | while read f; do grep -w "$f" "$card_aro" >> "$strain_1"_98p.genes; done

#strain 2
awk '{ if ( $7 > 0.97) print $0 }' "$card_2" > "$strain_2"_98p.bedcov
awk ' { print $1 } ' "$strain_2"_98p.bedcov | while read f; do grep -w "$f" "$card_aro" >> "$strain_2"_98p.genes; done

#find card genes in strain 1 and not strain 2
while read line; do
echo "$line" > tmp.file
  aro=$(awk '{print $1}' tmp.file)
  awk -v v1="$aro" '$0 ~ v1 {print}' "$strain_2"_98p.genes  &> /dev/null
  status=$? 
  if [ "$status" == 1 ]; then #CARD gene not found
    echo "card gene not found"
	echo "$line"
    echo "$line" >> "${strain_1}"_card.txt
  fi
done < "$strain_1"_98p.genes

#find card genes in strain 2 and not strain 1
while read line; do
echo "$line" > tmp.file
  aro=$(awk '{print $1}' tmp.file)
  awk -v v1="$aro" '$0 ~ v1 {print}' "$strain_1"_98p.genes  &> /dev/null
  status=$? 
  if [ "$status" == 1 ]; then #CARD gene not found
    echo "card gene not found"
	echo "$line"
    echo "$line" >> "${strain_2}"_card.txt
  fi
done < "$strain_2"_98p.genes


