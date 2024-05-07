#! /bin/bash


Ab=$1
ref="/home/dsarovich/bin/Pa_PAO1.fasta"
RESISTANCE_DB="Pseudomonas_aeruginosa_pao1.db"

current_dir=$(pwd)

##rescore
if [ -f "$current_dir"/ARDaP_full_summary.txt ]; then
  rm "$current_dir"/ARDaP_full_summary.txt
fi

if [ -s APV."$Ab".log ]; then
  rm APV."$Ab".log
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

echo "True ex int = $true_pos_true_neg"  >> APV.${Ab}.log
echo "False ex int = $false_pos_false_neg"  >> APV.${Ab}.log
echo "APV excluding intermediates is $APV_current" >> APV.${Ab}.log
echo "APV excluding intermediates is $APV_current"

if [ ! -z "$APV_previous" ]; then
  echo "RUNNING COMPARISON TEST"
  improved=$(echo "scale=4; $APV_current / $APV_previous > 1" | bc) #0 = no improvement or decreased APV
  if [ "$improved" -eq 1 ]; then
    echo -e "APV has improved - ACCEPT"
    echo -e "APV has improved - ACCEPT\n" >> APV.${Ab}.log
    APV_previous="$APV_current"
  else
    echo -e "APV hasn't improved - REJECT"
    echo -e "APV hasn't improved - REJECT\n" >> APV.${Ab}.log
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
echo "True int = $true_int" >> APV.${Ab}.log



false_int=0
while read line; do
false_int=$(echo "$false_int + $line" | bc)
done< <(grep "$Ab" $current_dir/ARDaP_full_summary.txt | grep "false intermediates" | awk '{print $6}')

echo "False int = $false_int"
echo "False int = $false_int" >> APV.${Ab}.log

true_pos_true_neg_true_int=$(echo "$true_pos + $true_neg + $true_int" | bc)
false_pos_false_neg_false_int=$(echo "$false_pos + $false_neg + $false_int" | bc)
total_int=$(echo "$true_pos_true_neg_true_int + $false_pos_false_neg_false_int" | bc)
APV_current_int=$(echo "$true_pos_true_neg_true_int / $total_int" | bc -l)


echo -e "APV including intermediates is $APV_current_int\n" >> APV.${Ab}.log
echo "True c int = $true_pos_true_neg_true_int"
echo "False c int = $false_pos_false_neg_false_int"
echo -e "APV including intermediates is $APV_current_int"

if [ ! -z "$APV_previous_int" ]; then
  echo "RUNNING COMPARISON TEST"
  improved=$(echo "scale=4; $APV_current_int / $APV_previous_int > 1" | bc) #0 = no improvement or decreased APV
  if [[ "$improved" -eq 1 ]]; then
    echo -e "APV int has improved - ACCEPT"
    echo -e "APV int has improved - ACCEPT\n" >> APV.${Ab}.log
    APV_previous_int="$APV_current_int"
  else
    echo -e "APV int hasn't improved - REJECT"
    echo -e "APV int hasn't improved - REJECT\n" >> APV.${Ab}.log
	  APV_previous_int="$APV_current_int"
  fi
else
    APV_previous_int="$APV_current_int"
fi

}

score

echo "scoring variants"
echo "Creating queries"

declare -A SQL_loss_report=()
STATEMENT_GENE_LOSS_COV () {
COUNTER=1
SQL_loss_report[$COUNTER]=$(
cat << EOF
SELECT
    Coverage.Gene,
    Coverage.Locus_tag,
    Coverage.Upregulation_or_loss,
    Coverage.Antibiotic_affected,
	Coverage.Threshold,
	Coverage.Comments
FROM
    Coverage
WHERE
    Coverage.Antibiotic_affected LIKE '%${Ab}r%';
EOF
)
}

STATEMENT_GENE_LOSS_COV
sqlite3 "$RESISTANCE_DB" "${SQL_loss_report[1]}" > db.R.dump

declare -A SQL_snps_report=()
STATEMENT_GENE_SNPS () {
COUNTER=1
SQL_snps_report[$COUNTER]=$(
cat << EOF
SELECT
    Variants_SNP_indel.Gene_name,
    Variants_SNP_indel.Locus_tag,
    Variants_SNP_indel.Variant_annotation,
    Variants_SNP_indel.Antibiotic_affected,
	Variants_SNP_indel.Threshold,
	Variants_SNP_indel.Comments
FROM
    Variants_SNP_indel
WHERE
    Variants_SNP_indel.Antibiotic_affected LIKE '%${Ab}r%';
EOF
)
}
STATEMENT_GENE_SNPS
sqlite3 "$RESISTANCE_DB" "${SQL_snps_report[1]}" >> db.R.dump

declare -A SQL_loss_report=()
STATEMENT_GENE_LOSS_COV () {
COUNTER=1
SQL_loss_report[$COUNTER]=$(
cat << EOF
SELECT
    Coverage.Gene,
    Coverage.Locus_tag,
    Coverage.Upregulation_or_loss,
    Coverage.Antibiotic_affected,
	Coverage.Threshold,
	Coverage.Comments
FROM
    Coverage
WHERE
    Coverage.Antibiotic_affected LIKE '%${Ab}i%';
EOF
)
}

STATEMENT_GENE_LOSS_COV
sqlite3 "$RESISTANCE_DB" "${SQL_loss_report[1]}" > db.I.dump

declare -A SQL_snps_report=()
STATEMENT_GENE_SNPS () {
COUNTER=1
SQL_snps_report[$COUNTER]=$(
cat << EOF
SELECT
    Variants_SNP_indel.Gene_name,
    Variants_SNP_indel.Locus_tag,
    Variants_SNP_indel.Variant_annotation,
    Variants_SNP_indel.Antibiotic_affected,
	Variants_SNP_indel.Threshold,
	Variants_SNP_indel.Comments
FROM
    Variants_SNP_indel
WHERE
    Variants_SNP_indel.Antibiotic_affected LIKE '%${Ab}i%';
EOF
)
}
STATEMENT_GENE_SNPS
sqlite3 "$RESISTANCE_DB" "${SQL_snps_report[1]}" >> db.I.dump


#
#
#   Need to pull out everything that isn't tagged with the antibiotic not none!!!!
#


declare -A SQL_loss_report=()
STATEMENT_GENE_LOSS_COV () {
COUNTER=1
SQL_loss_report[$COUNTER]=$(
cat << EOF
SELECT
    Coverage.Gene,
    Coverage.Locus_tag,
    Coverage.Upregulation_or_loss,
    Coverage.Antibiotic_affected,
	Coverage.Threshold,
	Coverage.Comments
FROM
    Coverage
WHERE
    Coverage.Antibiotic_affected NOT LIKE '%${Ab}%';
EOF
)
}

STATEMENT_GENE_LOSS_COV
sqlite3 "$RESISTANCE_DB" "${SQL_loss_report[1]}" > db.nat_var.dump

declare -A SQL_snps_report=()
STATEMENT_GENE_SNPS () {
COUNTER=1
SQL_snps_report[$COUNTER]=$(
cat << EOF
SELECT
    Variants_SNP_indel.Gene_name,
    Variants_SNP_indel.Locus_tag,
    Variants_SNP_indel.Variant_annotation,
    Variants_SNP_indel.Antibiotic_affected,
	Variants_SNP_indel.Threshold,
	Variants_SNP_indel.Comments
FROM
    Variants_SNP_indel
WHERE
    Variants_SNP_indel.Antibiotic_affected NOT LIKE '%${Ab}%';
EOF
)
}
STATEMENT_GENE_SNPS
sqlite3 "$RESISTANCE_DB" "${SQL_snps_report[1]}" >> db.nat_var.dump

#score variants

if [ -f ${Ab}_database_score.txt ]; then
 rm ${Ab}_database_score.txt ${Ab}_database_score_additional.txt
fi

echo "scoring performance of resistance markers"
#
# Each variant with "r" eg CIPr is analysed across all strains
#       and scored +2 if found in a resistant isolate
#           scored +1 if found in an intermediate isolate
#           and -8 if found in a sensitive isolate
#
#
#


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

echo "$short_var found in $found_R_strains resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_R_strains resistant" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_I_strains intermediate strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_I_strains intermediate strains\n" >> ${Ab}_database_score_additional.txt
echo -e "$score score for $short_var"  >> ${Ab}_database_score.txt
done < db.R.dump


#
# Each variant with "None" is analysed across all strains
#       and scored +2 if found in a resistant isolate
#           scored +1 if found in an intermediate isolate
#           and -8 if found in a sensitive isolate
#
#
##
# Each variant with "None" is analysed across all strains
#       and scored +1 if found in a resistant isolate
#           scored +2 if found in an intermediate isolate
#           and -4 if found in a sensitive isolate
#
#
#

echo "Looking for natural variation that may cause resistance or intermediate resistance"
while read variant; do
   found_R_nat_var=0
   found_S_nat_var=0
   found_I_nat_var=0
   not_found_R_nat_var=0
   not_found_S_nat_var=0
   not_found_I_nat_var=0
   nat_var_score=0

   found_R_nat_var_I=0
   found_S_nat_var_I=0
   found_I_nat_var_I=0
   #not_found_R_nat_var_I=0
   #not_found_S_nat_var_I=0
   #not_found_I_nat_var_I=0
   nat_var_I_score=0

   for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sherrard sun cortez tsang wardell; do
   cd ./"$i"/Outputs/AbR_reports/

   if [ -f ${Ab}.R.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_R_nat_var=$((found_R_nat_var+1))
		 found_R_nat_var_I=$((found_R_nat_var_I+1))
         nat_var_score=$((nat_var_score+2))
		 nat_var_I_score=$((nat_var_I_score+1))
      else
         not_found_R_nat_var=$((not_found_R_nat_var+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.R.txt
   fi

   if [ -f ${Ab}.S.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_S_nat_var=$((found_S_nat_var+1))
		 found_S_nat_var_I=$((found_S_nat_var_I+1))
         nat_var_score=$((nat_var_score-8))
		 nat_var_I_score=$((nat_var_I_score-4))
       else
         not_found_S_nat_var=$((not_found_S_nat_var+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.S.txt
   fi

   if [ -f ${Ab}.I.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_I_nat_var=$((found_I_nat_var+1))
		 found_I_nat_var_I=$((found_I_nat_var_I+1))
         nat_var_score=$((nat_var_score+1))
		 nat_var_I_score=$((nat_var_I_score+2))
       else
         not_found_I_nat_var=$((not_found_I_nat_var+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.I.txt
   fi
cd "$current_dir"

done

 echo "$variant" >tmp.variant
 short_var=$(awk -F"|" '{print $1,$2,$3,$4}' tmp.variant)

echo "$short_var found in $found_R_nat_var resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_R_nat_var resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_S_nat_var sensitive strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_S_nat_var sensitive strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var found in $found_I_nat_var intermediate strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_I_nat_var intermediate strains\n" >> ${Ab}_database_score_additional.txt
echo -e "$nat_var_score score for $short_var" >> ${Ab}_database_score.txt
echo -e "$nat_var_I_score score for $short_var" >> ${Ab}_database_int_score.txt
done < db.nat_var.dump

sort -nr -k1 ${Ab}_database_score.txt > ${Ab}_database_score.txt.tmp
mv ${Ab}_database_score.txt.tmp ${Ab}_database_score.txt

sort -nr -k1 ${Ab}_database_int_score.txt > ${Ab}_database_int_score.txt.tmp
mv ${Ab}_database_int_score.txt.tmp ${Ab}_database_int_score.txt


# Each variant with "i" e.g. CIPi is analysed across all strains
#       and scored +1 if found in a resistant isolate
#           scored +2 if found in an intermediate isolate
#           and -4 if found in a sensitive isolate
#
#
#

#
#
#  Intermediate variants are scored +3 if found in an intermediate strain
#                            scored -4 if found in a resistant strain
#                            scored -1 if found in a sensitive strain
#
echo "Looking at intermediate markers"
while read variant; do
   found_I_strains=0
   not_found_I_strains=0
   found_S_strains=0
   not_found_S_strains=0
   found_R_strains=0
   not_found_R_strains=0
   score_int=0

   for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
   cd ./"$i"/Outputs/AbR_reports/

   if [ -f ${Ab}.I.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt"  &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_I_strains=$((found_I_strains+1))
         score_int=$((score_int+3))
       else
         not_found_I_strains=$((not_found_I_strains+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.I.txt
   fi
   if [ -f ${Ab}.R.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_R_strains=$((found_R_strains+1))
         score_int=$((score_int-4))
		 #echo "found in $line"
       else
         not_found_R_strains=$((not_found_R_strains+1))

       fi
    done < ${Ab}.R.txt
   fi
   if [ -f ${Ab}.S.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_S_strains=$((found_S_strains+1))
         score_int=$((score_int-1))
       else
         not_found_S_strains=$((not_found_S_strains+1))
         #echo "not found in $strain"
       fi
    done < ${Ab}.S.txt
   fi
   cd "$current_dir"
done

 echo "$variant" >tmp.variant
 short_var=$(awk -F"|" '{print $1,$2,$3,$4}' tmp.variant)

echo "$short_var found in $found_R_strains resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_R_strains resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_I_strains intermediate strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_I_strains intermediate strains\n" >> ${Ab}_database_score_additional.txt
echo -e "$score_int score for $short_var"  >> ${Ab}_database_int_score.txt
done < db.I.dump

sort -nr -k1 ${Ab}_database_int_score.txt > ${Ab}_database_int_score.txt.tmp
mv ${Ab}_database_int_score.txt.tmp ${Ab}_database_int_score.txt

#
#     Scores resistance markers +2 if found in intermediate strains
#          scores resistance markers -4 if found in resistant isolates
#          scores resistance markers -1 if found in sensitive isolates
#
#

echo "Looking at resistance markers that might be intermediate markers"

while read variant; do
   found_I_strains=0
   not_found_I_strains=0
   found_S_strains=0
   not_found_S_strains=0
   found_R_strains=0
   not_found_R_strains=0
   score_int=0

   for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
   cd ./"$i"/Outputs/AbR_reports/

   if [ -f ${Ab}.I.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt"  &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_I_strains=$((found_I_strains+1))
         score_int=$((score_int+3))
       else
         not_found_I_strains=$((not_found_I_strains+1))
         #echo "found in $strain"
       fi
    done < ${Ab}.I.txt
   fi
   if [ -f ${Ab}.R.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_R_strains=$((found_R_strains+1))
         score_int=$((score_int-4))
		 #echo "found in $line"
       else
         not_found_R_strains=$((not_found_R_strains+1))

       fi
    done < ${Ab}.R.txt
   fi
   if [ -f ${Ab}.S.txt ]; then
    while read line; do
       grep "$variant" "${line}.AbR_output.final.txt" &> /dev/null
       status=$?
       if [ "$status" == 0 ]; then
         found_S_strains=$((found_S_strains+1))
         score_int=$((score_int-1))
       else
         not_found_S_strains=$((not_found_S_strains+1))
         #echo "not found in $strain"
       fi
    done < ${Ab}.S.txt
   fi
   cd "$current_dir"
done

 echo "$variant" >tmp.variant
 short_var=$(awk -F"|" '{print $1,$2,$3,$4}' tmp.variant)

echo "$short_var found in $found_R_strains resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_R_strains resistant strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo "$short_var not found in $not_found_S_strains sensitive strains" >> ${Ab}_database_score_additional.txt
echo "$short_var found in $found_I_strains intermediate strains" >> ${Ab}_database_score_additional.txt
echo -e "$short_var not found in $not_found_I_strains intermediate strains\n" >> ${Ab}_database_score_additional.txt
echo -e "$score_int score for $short_var"  >> ${Ab}_database_int_score_for_R_markers.txt
done < db.R.dump

sort -nr -k1 ${Ab}_database_int_score_for_R_markers.txt > ${Ab}_database_int_score_for_R_markers.txt.tmp
mv ${Ab}_database_int_score_for_R_markers.txt.tmp ${Ab}_database_int_score_for_R_markers.txt


#running database cleanup

#top ten hits if they exist

grep -v "${Ab}r" ${Ab}_database_score.txt | grep -v "${Ab}i" | head -n10  | awk -F" " ' $1 > 5 '  |  awk '{print $4,$5,$6 }' | while read line; do
grep "$line" ${Ab}_database_score_additional.txt >> ${Ab}_natural_variation_to_add.txt;
done

grep -e "${Ab}r" -e "${Ab}i" ${Ab}_database_score.txt | tail -n10 | awk -F" " ' $1 < 0 ' |  awk '{print $4,$5,$6 }' | while read line; do
grep "$line" ${Ab}_database_score_additional.txt >> ${Ab}_variants_causing_false_pos.txt;
done

grep -v 'None' ${Ab}_database_int_score.txt  | tail -n10 | awk -F" " ' $1 < 0 ' |  awk '{print $4,$5,$6 }' | while read line; do
grep "$line" ${Ab}_database_score_additional.txt >> ${Ab}_variants_causing_false_pos_intermediate.txt;
done


#rescore

#removing false pos
echo "Removing false positives and rescoring"
echo "Removing false positives and rescoring" >> APV.${Ab}.log
sed -i 's/, /,/g' ${Ab}_database_score.txt
grep "${Ab}r" ${Ab}_database_score.txt | tail -n10 | awk -F" " ' $1 < 0 ' |  tac | awk -F" " '{print $4"|"$5"|"$6"|" }' | while read line; do
  if [ ! -z "$line" ]; then
    echo "removing $line due to flagging as a false positive"
    echo "removing $line due to flagging as a false positive" >> APV.${Ab}.log
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


#Adding natural variation
echo "Adding natural variation and rescoring"
echo "Adding natural variation and rescoring"  >> APV.${Ab}.log
grep -v "${Ab}r" ${Ab}_database_score.txt | grep -v "${Ab}i"  | head -n10  | awk -F" " ' $1 > 5 '  |  awk '{print $4"|"$5"|"$6"|"$7 }' | while read line; do
  echo "Adding $line"
  echo "Adding $line" >> APV.${Ab}.log
  if [ ! -z "$line" ]; then
    for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
      cd "$current_dir"/"$i"/Outputs/AbR_reports/
      for f in *.AbR_output.final.txt; do
        grep "$line" "$f" | awk -F"|" -v Ab="$Ab" 'BEGIN { OFS="|" }  {print $1,$2,$3,Ab"r,"$4,$5,$6}' >> "$f".tmp
        grep -v "$line" "$f" >> "$f".tmp
        mv "$f".tmp "$f"
      done
    done
    cd "$current_dir"
    rm "$current_dir"/ARDaP_full_summary.txt
    score #run score function
  fi
done


#____________________________________
#
#       Intermediate opt
#_____________________________________
echo "removing false pos intermediate markers"
echo -e "\nremoving false pos intermediate markers"  >> APV.${Ab}.log
sed -i 's/, /,/g' ${Ab}_database_int_score.txt
grep "${Ab}i" ${Ab}_database_int_score.txt | tail -n10 | awk -F" " ' $1 < 0 ' |  awk -F" " '{print $4"|"$5"|"$6"|" }' | while read line; do
if [ ! -z "$line" ]; then
    echo "removing $line due to flagging as a false positive intermediate"
    echo "removing $line due to flagging as a false positive intermediate" >> APV.${Ab}.log
    for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
      cd "$current_dir"/"$i"/Outputs/AbR_reports/
        for f in *.AbR_output.final.txt; do
          grep "$line" "$f" | sed "s/${Ab}i/None/" >> "$f".tmp
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


echo -e "\nChanging resistance markers to intermediate markers" >> APV.${Ab}.log
sed -i 's/, /,/g' ${Ab}_database_int_score.txt
grep "${Ab}r" ${Ab}_database_int_score_for_R_markers.txt | head -n10 | awk -F" " ' $1 > 0 ' |  awk -F" " '{print $4"|"$5"|"$6"|" }' | while read line; do
if [ ! -z "$line" ]; then
    echo "changing $line into intermediate"
    echo "changing $line into intermediate" >> APV.${Ab}.log
    for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
      cd "$current_dir"/"$i"/Outputs/AbR_reports/
        for f in *.AbR_output.final.txt; do
          grep "$line" "$f" | sed "s/${Ab}r/${Ab}i/" >> "$f".tmp
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

#Adding natural variation
echo "Adding natural variation for intermediate resistance and rescoring"
echo "Adding natural variation for intermediate resistance and rescoring"  >> APV.${Ab}.log
grep -v "${Ab}r" ${Ab}_database_int_score.txt | grep -v "${Ab}i"  | head -n10  | awk -F" " ' $1 > 5 '  |  awk '{print $4"|"$5"|"$6"|" }' | while read line; do
  echo "Adding $line"
  echo "Adding $line" >> APV.${Ab}.log
  if [ ! -z "$line" ]; then
    for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan sun cortez sherrard tsang wardell; do
      cd "$current_dir"/"$i"/Outputs/AbR_reports/
      for f in *.AbR_output.final.txt; do
        grep "$line" "$f" | awk -F"|" -v Ab="$Ab" 'BEGIN { OFS="|" } {print $1,$2,$3,Ab"i,"$4,$5,$6}' >> "$f".tmp
        grep -v "$line" "$f" >> "$f".tmp
        mv "$f".tmp "$f"
      done
    done
    cd "$current_dir"
    rm "$current_dir"/ARDaP_full_summary.txt
    score #run score function
  fi
done
