
abbrev=$1

antibiotic=$2


while read file; do

echo "formatting file $file"

file="$file".AbR_output.GWAS_SNP.txt
#make file tab seperated
sed -i 's/Genomic_coord/Chromosome|coord/g' "$file"
sed -i 's/Chromosome_/Chromosome|/g' "$file"
sed -i 's/|/\t/g' "$file"
sed -i 's/ I /_I_/g' "$file"
sed -i 's/ R /_R_/g' "$file"
sed -i 's/_p-value//g' "$file"
sed -i 's/ /\t/g' "$file"

#find antibiotic in gwas table for resistant values
column=$(head -n1 $file | sed 's/\t/\n/g' | nl | grep "$antibiotic"_R | awk '{print $1}')

cat "$file" | sort -t$'\t' -gk"$column" | head -n21 >"$file"."$antibiotic"_R

done < <(awk -F"," '{print $1 }' Pipeline_False_Negatives."$abbrev".txt)


echo "creating summary"

cat *."$antibiotic"_R > "$antibiotic"_R_summary.txt

file=$(head -n1 Pipeline_False_Negatives."$abbrev".txt | awk -F"," '{ print $1}' )
file="$file".AbR_output.GWAS_SNP.txt

echo "counting hits"
gene_name_col=$(head -n1 $file | sed 's/\t/\n/g' | nl | grep Gene_name | awk '{print $1}')
echo "gene column is $gene_name_col"
mutation_col=$(head -n1 $file | sed 's/\t/\n/g' | nl | grep -w Mutation | awk '{print $1}')
echo "mutation col is $mutation_col"
grep -v "GWASID" "$antibiotic"_R_summary.txt | awk -v x="$gene_name_col" -v y="$mutation_col" '{ print $x, $y }' | sort | uniq -c > number_of_hits.txt

sort -nrk1 number_of_hits.txt | tail -n+2 >number_of_hits.txt.tmp

mv number_of_hits.txt.tmp number_of_hits.txt

if [ -s variant_counts.txt ]; then
  rm variant_counts.txt
fi  

echo "looking through resistant/senstive strains"

while read line; do

  echo $line > line.tmp
  gene=$(awk '{print $2 }' line.tmp)
  var=$(awk '{print $3 }' line.tmp)

  #count variant in resistant strains
  count=0 
  while read strain; do
    
    cat ../../../Variants/Annotated/$strain.PASS.snps.annotated.vcf | grep "$gene" | grep "$var" > /dev/null
	status=$?
	if [[ "$status" -eq 0 ]]; then
	count=$((count+1))
    fi
	
  done < ../../"$abbrev".R.txt
  
  echo "$gene $var found in $count resistant strains" >> variant_counts.txt
  
  
  #count variant in sensitive strains
  count=0 
  while read strain; do
    
    cat ../../../Variants/Annotated/$strain.PASS.snps.annotated.vcf | grep "$gene" | grep "$var" > /dev/null
	status=$?
	if [[ "$status" -eq 0 ]]; then
	count=$((count+1))
    fi
	
  done < ../../"$abbrev".S.txt
  
  echo "$gene $var found in $count sensitive strains" >> variant_counts.txt
  
done < number_of_hits.txt

exit 0