#! /bin/bash

for Ab in MEM FEP CAZ TZP PIP CST CIP AMK TOB; do
ref="/home/dsarovich/bin/Pa_PAO1.fasta"
RESISTANCE_DB="Pseudomonas_aeruginosa_pao1.db"

current_dir=$(pwd)

##rescore
if [ -f "$current_dir"/ARDaP_full_summary.txt ]; then
  rm "$current_dir"/ARDaP_full_summary.txt
fi


score () {
"$current_dir"/summarise_all_ardap.sh

##Score antibiotic
#false pos
false_pos=0
while read line; do
false_pos=$(echo "$false_pos + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false positive" | awk '{print $6}')

#echo "False pos = $false_pos"

true_neg_c_int=0
while read line; do
true_neg_c_int=$(echo "$true_neg_c_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "Sensitive strains identified as intermediate" | awk '{print $7}')


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

false_neg_c_int=0
while read line; do
false_neg_c_int=$(echo "$false_neg_c_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "Resistant strains identified as intermediate" | awk '{print $7}')

#echo "True pos c int = $true_pos_c_int"

#APV
true_pos_true_neg=$(echo "$true_pos + $true_neg + $true_neg_c_int" | bc)
false_pos_false_neg=$(echo "$false_pos + $false_neg + $false_neg_c_int" | bc)
total=$(echo "$true_pos_true_neg + $false_pos_false_neg" | bc)
APV_current=$(echo "$true_pos_true_neg / $total" | bc -l)

echo "True ex int = $true_pos_true_neg"  >> APV.${Ab}_Resfinder.log
echo "False ex int = $false_pos_false_neg"  >> APV.${Ab}_Resfinder.log
echo "APV excluding intermediates is $APV_current" >> APV.${Ab}_Resfinder.log
echo "APV excluding intermediates is $APV_current"

if [ ! -z "$APV_previous" ]; then
  echo "RUNNING COMPARISON TEST"
  improved=$(echo "scale=4; $APV_current / $APV_previous > 1" | bc) #0 = no improvement or decreased APV
  if [ "$improved" -eq 1 ]; then
    echo -e "APV has improved - ACCEPT"
    echo -e "APV has improved - ACCEPT\n" >> APV.${Ab}_Resfinder.log
    APV_previous="$APV_current"
  else
    echo -e "APV hasn't improved - REJECT"
    echo -e "APV hasn't improved - REJECT\n" >> APV.${Ab}_Resfinder.log
  	APV_previous="$APV_current"
  fi
else
  APV_previous="$APV_current"
fi



#intermediates
true_int=0
while read line; do
true_int=$(echo "$true_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "true intermediates" | awk '{print $6}')

echo "True int = $true_int"
echo "True int = $true_int" >> APV.${Ab}_Resfinder.log



false_int=0
while read line; do
false_int=$(echo "$false_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false intermediates" | awk '{print $6}')

echo "False int = $false_int"
echo "False int = $false_int" >> APV.${Ab}_Resfinder.log

true_pos_true_neg_true_int=$(echo "$true_pos + $true_neg + $true_int" | bc)
false_pos_false_neg_false_int=$(echo "$false_pos + $false_neg + $false_int" | bc)
total_int=$(echo "$true_pos_true_neg_true_int + $false_pos_false_neg_false_int" | bc)
APV_current_int=$(echo "$true_pos_true_neg_true_int / $total_int" | bc -l)


echo -e "APV including intermediates is $APV_current_int\n" >> APV.${Ab}_Resfinder.log
echo "True c int = $true_pos_true_neg_true_int"
echo "False c int = $false_pos_false_neg_false_int"
echo -e "APV including intermediates is $APV_current_int"

if [ ! -z "$APV_previous_int" ]; then
  echo "RUNNING COMPARISON TEST"
  improved=$(echo "scale=4; $APV_current_int / $APV_previous_int > 1" | bc) #0 = no improvement or decreased APV
  if [[ "$improved" -eq 1 ]]; then
    echo -e "APV int has improved - ACCEPT"
    echo -e "APV int has improved - ACCEPT\n" >> APV.${Ab}_Resfinder.log
    APV_previous_int="$APV_current_int"
  else
    echo -e "APV int hasn't improved - REJECT"
    echo -e "APV int hasn't improved - REJECT\n" >> APV.${Ab}_Resfinder.log
	  APV_previous_int="$APV_current_int"
  fi
else
    APV_previous_int="$APV_current_int"
fi

}

#Resfinder optimize
echo "dump output of resfinder"
if [ -s Resfinder.variants.txt ]; then
  rm Resfinder.variants.txt
fi

for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sherrard sun cortez tsang wardell; do
grep 'Resfinder' ./"$i"/Outputs/AbR_reports/*.AbR_output.final.txt | grep "$Ab" | awk -F"|" '{print $2}' | sed 's/,/\n/g' >> Resfinder.variants.txt
done
sort Resfinder.variants.txt > Resfinder.variants.txt.tmp
uniq Resfinder.variants.txt.tmp Resfinder.variants.txt




while read variant; do
   found_R_strains=0
   not_found_R_strains=0
   found_S_strains=0
   not_found_S_strains=0
   found_I_strains=0
   not_found_I_strains=0
   score=0

   for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sherrard sun cortez tsang wardell; do
   cd ./"$i"/Outputs/AbR_reports/

   if [ -f ${Ab}.R.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_R_strains=$((found_R_strains+1))
         score=$((score+2)) #Variant found in resistant
       else
         not_found_R_strains=$((not_found_R_strains+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.R.txt
   fi
   if [ -f ${Ab}.S.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_S_strains=$((found_S_strains+1))
         score=$((score-8)) #Variant found in sensitive
		 echo "$line" >> "$current_dir"/"$Ab".sensitive.list.txt
       else
         not_found_S_strains=$((not_found_S_strains+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.S.txt
   fi
   if [ -f ${Ab}.I.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_I_strains=$((found_I_strains+1))
         score=$((score+1)) #Variant found in intermediate
       else
         not_found_I_strains=$((not_found_I_strains+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.I.txt
   fi
   cd "$current_dir"
done

 echo "$variant" >tmp.variant
 short_var=$(awk -F"|" '{print $1,$2,$3,$4}' tmp.variant)

echo "$variant found in $found_R_strains resistant strains" >> ${Ab}_Resfinder.txt
echo "$variant not found in $not_found_R_strains resistant" >> ${Ab}_Resfinder.txt
echo "$variant found in $found_S_strains sensitive strains" >> ${Ab}_Resfinder.txt
echo -e "$variant not found in $not_found_S_strains sensitive strains" >> ${Ab}_Resfinder.txt
echo "$variant found in $found_I_strains intermediate strains" >> ${Ab}_Resfinder.txt
echo -e "$variant not found in $not_found_I_strains intermediate strains\n" >> ${Ab}_Resfinder.txt
echo -e "$score score for $variant"  >> ${Ab}_database_score_Resfinder.txt
done < Resfinder.variants.txt

sort -n -k1 ${Ab}_database_score_Resfinder.txt > ${Ab}_database_score_Resfinder.txt.tmp
mv ${Ab}_database_score_Resfinder.txt.tmp ${Ab}_database_score_Resfinder.txt


score

##Remove false positive hits

#removing false pos
echo "Removing false positives and rescoring"
echo "Removing false positives and rescoring" >> APV.${Ab}_Resfinder.log
#sed 's/, /,/g' ${Ab}_database_score_Resfinder.txt
head ${Ab}_database_score_Resfinder.txt | awk -F" " ' $1 < 0 ' |  awk -F" " '{print $5 }' | while read line; do
  if [ ! -z "$line" ]; then
    echo "removing $line due to flagging as a false positive"
    echo "removing $line due to flagging as a false positive" >> APV.${Ab}_Resfinder.log
    for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
      cd "$current_dir"/"$i"/Outputs/AbR_reports/
        for f in *.AbR_output.final.txt; do
          grep "$line" "$f" | sed "s/${Ab}r/None/" >> "$f".tmp
          grep -v "$line" "$f" >> "$f".tmp
          mv "$f".tmp "$f"
        done
    done
    cd "$current_dir"
    #reanalyse
	rm "$current_dir"/ARDaP_full_summary.txt

	score #run score function

  fi
done




done
