#! /bin/bash


Ab=$1
ref="/home/dsarovich/bin/Pa_PAO1.fasta"
RESISTANCE_DB="Pseudomonas_aeruginosa_pao1.db"

current_dir=$(pwd)

##rescore

score () {
if [ -f "$current_dir"/ARDaP_full_summary.txt ]; then
  rm "$current_dir"/ARDaP_full_summary.txt
fi


"$current_dir"/summarise_all_ardap.sh

##Score antibiotic
#false pos
false_pos=0
while read line; do
false_pos=$(echo "$false_pos + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false positive" | awk '{print $6}')

#echo "False pos = $false_pos"

false_pos_c_int=0
while read line; do
false_pos_c_int=$(echo "$false_pos_c_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "Sensitive strains identified as intermediate" | awk '{print $7}')


#echo "False pos c int = $false_pos_c_int"

false_neg=0
while read line; do
false_neg=$(echo "$false_neg + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false negatives" | awk '{print $6}')

#echo "False neg = $false_neg"


true_neg=0
while read line; do
true_neg=$(echo "$true_neg + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "true negatives" | awk '{print $6}')

#echo "True neg = $true_neg"


true_pos=0
while read line; do
true_pos=$(echo "$true_pos + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "true positives" | awk '{print $6}')

#echo "True pos = $true_pos"

true_pos_c_int=0
while read line; do
true_pos_c_int=$(echo "$true_pos_c_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "Resistant strains identified as intermediate" | awk '{print $7}')

#echo "True pos c int = $true_pos_c_int"

#APV
true_pos_true_neg=$(echo "$true_pos + $true_neg + $true_pos_c_int" | bc)
false_pos_false_neg=$(echo "$false_pos + $false_neg + $false_pos_c_int" | bc)
total=$(echo "$true_pos_true_neg + $false_pos_false_neg" | bc)

APV_current=$(echo "$true_pos_true_neg / $total" | bc -l)
echo "True ex int = $true_pos_true_neg"
echo "False ex int = $false_pos_false_neg"
echo "APV excluding intermediates is $APV_current" >> APV.${Ab}.indel.log
echo "APV excluding intermediates is $APV_current" 

#intermediates
true_int=0
while read line; do
true_int=$(echo "$true_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "true intermediates" | awk '{print $6}')

echo "True int = $true_int" 
echo "True int = $true_int" >> APV.${Ab}.indel.log



false_int=0
while read line; do
false_int=$(echo "$false_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false intermediates" | awk '{print $6}')

echo "False int = $false_int" 
echo "False int = $false_int" >> APV.${Ab}.indel.log

true_pos_true_neg_true_int=$(echo "$true_pos + $true_neg + $true_int" | bc)
false_pos_false_neg_false_int=$(echo "$false_pos + $false_neg + $false_int" | bc)
total_int=$(echo "$true_pos_true_neg_true_int + $false_pos_false_neg_false_int" | bc)
APV_current_int=$(echo "$true_pos_true_neg_true_int / $total_int" | bc -l)
echo "APV including intermediates is $APV_current_int" >> APV.${Ab}.indel.log
echo "True c int = $true_pos_true_neg_true_int"
echo "False c int = $false_pos_false_neg_false_int"
echo "APV including intermediates is $APV_current_int" 

}


score

#/home/dsarovich/analyses/pa_ardap/"$i"/Outputs/Variants/Annotated/*.vcf

#for each variant, go through each directory and add to strains if identified
#for snps, format of line =
#parC|PA4964|Ser87Trp|CIPr|100|

 #barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sherrard tsang wardel
cat indel.scores.full.$Ab | while read line; do 
  echo "$line" >> APV.${Ab}.indel.log
  echo "$line" > tmp.var
  gene=$(awk '{print $1}' tmp.var)
  var=$(awk '{print $2}' tmp.var)
  for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sherrard tsang wardel; do
    if [ -f ./"$i"/Outputs/AbR_reports/${Ab}.R.txt ]; then
     for strain in ./"$i"/Outputs/AbR_reports/*.AbR_output.final.txt; do
      strain_clean=${strain/.\/"$i"\/Outputs\/AbR_reports\//}
      #echo $strain
      grep "$gene" ./$i/Outputs/Variants/Annotated/${strain_clean/AbR_output.final.txt/PASS.indels.annotated.vcf} | grep "$var" &> /dev/null         
      status=$? 
     if [ "$status" == 0 ]; then
        echo "$gene|$gene|$var|${Ab}r" >> $strain
        #echo "found"
       # echo "$strain_clean"
      #else
       #echo "not found"
        #echo "$strain_clean"
      fi
     done
	fi 
  done
  score
done
