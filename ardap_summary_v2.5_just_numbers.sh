#! /bin/bash

#echo -e "This script will determine the rate of false positives, false negatives, true positive and true negatives following a complete ARDaP run.\n\n"
#echo -e "The false positive assessment for coverage requires the reference genome in the same directory and specified in the command line"

# The script requires the following files

DATASHEET=/home/dsarovich2/bin/ARDaP/Databases/Pseudomonas_aeruginosa/Pa_res_data.txt

ref=/home/dsarovich/bin/Pa_PAO1.fasta
RESISTANCE_DB=/home/dsarovich/bin/Pseudomonas_aeruginosa_pao1.db
false_pos_check=0
false_neg_check=0


Antibiotic=$1

#declare variables

FalseNegative=0
FalsePositive=0
TrueNegative=0
TruePositive=0
FalseIntermediate=0
TrueIntermediate=0
FalseIntermediate_res_FN=0
FalseIntermediate_sens_FP=0
FalseIntermediate_int_FN=0
FalseIntermediate_int_FP=0
## Logic flow

# for each antibiotic; determine the following numbers

# false negative is a strain that is in the resistant list but called as sensitive

# false positive is a strain that is in the sensitive list but called as resistant

# true negative is a strain that is in the sensitive list and called as sensitive

# true positive is a strain that is in the resistant list and called as resistant

# Load the list of sensitive strains for each antibiotic



if [ ! -s "$Antibiotic".R.txt ]; then
   echo "Cannot find file listing resistant strains"
   echo "I am looking for this file --> ${Antibiotic}.R.txt"
   echo "Please check and re-run if there should be resistant strains in this dataset"
  # exit 1
fi
if [ ! -s "$Antibiotic".S.txt ]; then
   echo "Cannot find file listing resistant strains"
   echo "I am looking for this file --> ${Antibiotic}.S.txt"
   echo "Please check and re-run if there should be sensitive strains in this dataset"
  # exit 1
fi

if [ ! -s ./summary ]; then
  mkdir ./summary
fi

cd ./summary

#cleanup previous run
if [ -s Pipeline_True_Positives.${Antibiotic}.txt ]; then
  rm Pipeline_True_Positives.${Antibiotic}.txt
fi
if [ -s Pipeline_True_Negatives.${Antibiotic}.txt ]; then
  rm Pipeline_True_Negatives.${Antibiotic}.txt
fi
if [ -s Pipeline_True_Intermediate.${Antibiotic}.txt ]; then
  rm Pipeline_True_Intermediate.${Antibiotic}.txt
fi
if [ -s Pipeline_False_Positives.${Antibiotic}.FULL.txt ]; then
  rm Pipeline_False_Positives.${Antibiotic}.FULL.txt
fi
if [ -s Pipeline_False_Positives.${Antibiotic}.COUNTS.txt ]; then
  rm Pipeline_False_Positives.${Antibiotic}.COUNTS.txt
fi
if [ -s Pipeline_False_Positives.${Antibiotic}.txt ]; then
  rm Pipeline_False_Positives.${Antibiotic}.txt
fi
if [ -s Pipeline_False_Negatives.${Antibiotic}.txt ]; then
  rm Pipeline_False_Negatives.${Antibiotic}.txt
fi
if [ -s Pipeline_False_Intermediate.${Antibiotic}.txt ]; then
  rm Pipeline_False_Intermediate.${Antibiotic}.txt
fi
if [ -s False.positive.loss.variation.${Antibiotic}.txt ]; then
  rm False.positive.loss.variation.${Antibiotic}.txt
fi
 
#echo -e "Running pipeline assessment for $Antibiotic\n\n"
while read SensitiveStrain; do
    if [ ! -s ../${SensitiveStrain}.AbR_output.final.txt ]; then
      echo -e "Cannot find the AbR output for the sensitive strain ${SensitiveStrain}. Did ARDaP complete normally?"
      echo -e "Please check and re-run"
      echo "Exiting"
      exit 1
    else 
      # echo -e "Testing $SensitiveStrain\n"
	  Antibiotic_res=${Antibiotic}r
      awk -F "|" -v Ab=$Antibiotic_res '$4 ~ Ab { exit 1 }' ../${SensitiveStrain}.AbR_output.final.txt &> /dev/null
	  status=$? 
	  if [ "$status" == 1 ]; then #sensitive strain but called as resistant
	     FalsePositive=$((FalsePositive+1))
	     echo -e "$SensitiveStrain,False Positive for $Antibiotic" >> Pipeline_False_Positives.${Antibiotic}.txt
		 echo -e "$SensitiveStrain,False Positive for $Antibiotic" >> Pipeline_False_Positives.${Antibiotic}.FULL.txt
		 cat ../${SensitiveStrain}.AbR_output.final.txt >> Pipeline_False_Positives.${Antibiotic}.FULL.txt
		 echo -e "\n" >> Pipeline_False_Positives.${Antibiotic}.FULL.txt
	  else #look for intermediate
	     Antibiotic_res=${Antibiotic}i
	     awk -F "|" -v Ab=$Antibiotic_res '$4 ~ Ab { exit 1 }' ../${SensitiveStrain}.AbR_output.final.txt &> /dev/null
		 status=$?
		 if [ "$status" == 1 ]; then #sensitive strain but called as intermediate
           FalseIntermediate=$((FalseIntermediate+1)) #Sensitive strain called as Intermediate	
		   FalseIntermediate_sens_FP=$((FalseIntermediate_sens_FP+1))
	       echo -e "$SensitiveStrain, Sensitive strain identified as intermediate for $Antibiotic" >> Pipeline_False_Intermediate.sensitive.FP.${Antibiotic}.txt
		   echo -e "$SensitiveStrain, Sensitive strain, called as intermediate for $Antibiotic" >> Pipeline_False_Intermediate.sensitive.FP.${Antibiotic}.FULL.txt
		   cat ../${SensitiveStrain}.AbR_output.final.txt >> Pipeline_False_Intermediate.sensitive.FP.${Antibiotic}.FULL.txt
		   echo -e "\n" >> Pipeline_False_Intermediate.sensitive.FP.${Antibiotic}.FULL.txt
		 else 
		   TrueNegative=$((TrueNegative+1)) #Sensitive strain not called as resistant or intermediate
           echo -e "$SensitiveStrain,True Negative for $Antibiotic" >> Pipeline_True_Negatives.${Antibiotic}.txt		 
		 fi
	   fi
    fi
done < ../"$Antibiotic".S.txt


if [ -s ../"$Antibiotic".R.txt ]; then
  while read ResistantStrain; do
    if [ ! -s ../${ResistantStrain}.AbR_output.final.txt ]; then
	  echo -e "Cannot find AbR output for the resistant strain ${ResistantStrain}"
      echo -e "Please check and re-run"
      echo "Exiting"
      exit 1
	else
     # echo -e "Testing $ResistantStrain\n"	
	  Antibiotic_res=${Antibiotic}r
      awk -F "|" -v Ab=$Antibiotic_res '$4 ~ Ab { exit 1} ' ../${ResistantStrain}.AbR_output.final.txt &> /dev/null
	  status=$? 
	  if [ "$status" == 1 ]; then  #Resistant strain called as resistant
	   TruePositive=$((TruePositive+1))
	   echo -e "$ResistantStrain,True Positive for $Antibiotic" >> Pipeline_True_Positives.${Antibiotic}.txt
	  else
	   Antibiotic_int=${Antibiotic}i
	   awk -F "|" -v Ab="$Antibiotic_int" '$4 ~ Ab { exit 1} ' ../${ResistantStrain}.AbR_output.final.txt &> /dev/null
	   status=$? 
	   if [ "$status" == 1 ]; then
	     FalseIntermediate=$((FalseIntermediate+1)) #resistant strain identified as intermediate
		 FalseIntermediate_res_FN=$((FalseIntermediate_res_FN+1))
	     echo -e "$ResistantStrain, Resistant strain identified as intermediate for $Antibiotic" >> Pipeline_False_Intermediate.resistant.FN.${Antibiotic}.txt
		 echo -e "$ResistantStrain, Resistant strain identified as intermediate for $Antibiotic" >> Pipeline_False_Intermediate.resistant.FN.${Antibiotic}.FULL.txt
		 cat ../${ResistantStrain}.AbR_output.final.txt >> Pipeline_False_Intermediate.resistant.FN.${Antibiotic}.FULL.txt
		 echo -e "\n" >> Pipeline_False_Intermediate.resistant.FN.${Antibiotic}.FULL.txt
       else
	     FalseNegative=$((FalseNegative+1)) #Resistant strain not called as resistant	
	     echo -e "$ResistantStrain,False Negative for $Antibiotic" >> Pipeline_False_Negatives.${Antibiotic}.txt
		 echo -e "$ResistantStrain,False Negative for $Antibiotic" >> Pipeline_False_Negatives.${Antibiotic}.FULL.txt
		 cat ../${ResistantStrain}.AbR_output.final.txt >> Pipeline_False_Negatives.${Antibiotic}.FULL.txt
		 echo -e "\n" >> Pipeline_False_Negatives.${Antibiotic}.FULL.txt
       fi	  
	  fi 
	fi  
  done < ../"$Antibiotic".R.txt
fi

if [ -s ../"$Antibiotic".I.txt ]; then
    echo "found list of intermediate strains"
    while read IntermediateStrain; do
        if [ ! -s ../${IntermediateStrain}.AbR_output.final.txt ]; then
	        echo -e "Cannot find AbR output for the Intermediate strain ${IntermediateStrain}"
            echo -e "Please check and re-run"
            echo "Exiting"
            exit 1
	    else
            # echo -e "Testing $IntermediateStrain\n"	
			Antibiotic_int=${Antibiotic}i
            awk -F "|" -v Ab=$Antibiotic_int '$4 ~ Ab { exit 1}' ../${IntermediateStrain}.AbR_output.final.txt &> /dev/null
	        status=$? 
	        if [ "$status" == 1 ]; then  #Intermediate strain called as Intermediate
	            TrueIntermediate=$((TrueIntermediate+1))
	            echo -e "$IntermediateStrain,True Intermediate for $Antibiotic" >> Pipeline_True_Intermediate.${Antibiotic}.txt
				echo -e "$IntermediateStrain,True Intermediate for $Antibiotic" >> Pipeline_True_Intermediate.${Antibiotic}.FULL.txt
		        cat ../${IntermediateStrain}.AbR_output.final.txt >> Pipeline_True_Intermediate.${Antibiotic}.FULL.txt
				echo -e "\n" >> Pipeline_True_Intermediate.${Antibiotic}.FULL.txt
	        else
                Antibiotic_int=${Antibiotic}r
				awk -F "|" -v Ab=$Antibiotic_int '$4 ~ Ab { exit 1}' ../${IntermediateStrain}.AbR_output.final.txt &> /dev/null
				status=$?
				if [ "$status" == 1 ]; then  #Intermediate strain called as resistant
				  FalseIntermediate=$((FalseIntermediate+1)) #Intermediate strain identified as resistant	
				  FalseIntermediate_int_FP=$((FalseIntermediate_int_FP+1))
	              echo -e "$IntermediateStrain, Intermediate strain identified as resistant for $Antibiotic" >> Pipeline_False_Intermediate.intermediate.FP.${Antibiotic}.txt
				  echo -e "$IntermediateStrain, Intermediate strain, called as resistance for $Antibiotic" >> Pipeline_False_Intermediate.intermediate.FP.${Antibiotic}.FULL.txt
		          cat ../${IntermediateStrain}.AbR_output.final.txt >> Pipeline_False_Intermediate.intermediate.FP.${Antibiotic}.FULL.txt
				  echo -e "\n" >> Pipeline_False_Intermediate.intermediate.FP.${Antibiotic}.FULL.txt
				else
				  FalseIntermediate=$((FalseIntermediate+1))
				  FalseIntermediate_int_FN=$((FalseIntermediate_int_FN+1))
				  echo -e "$IntermediateStrain, Intermediate strain, identified as sensitive for $Antibiotic" >> Pipeline_False_Intermediate.intermediate.FN.${Antibiotic}.txt
                  echo -e "$IntermediateStrain, Intermediate strain, called as sensitive for $Antibiotic" >> Pipeline_False_Intermediate.intermediate.FN.${Antibiotic}.FULL.txt
		          cat ../${IntermediateStrain}.AbR_output.final.txt >> Pipeline_False_Intermediate.intermediate.FN.${Antibiotic}.FULL.txt	
                  echo -e "\n" >> Pipeline_False_Intermediate.intermediate.FN.${Antibiotic}.FULL.txt			  
				fi
	        fi 
	    fi  
    done < ../"$Antibiotic".I.txt
fi


#false positive interrogation set to 1 in header to run
if [ "$false_pos_check" == 1 ]; then
sort Pipeline_False_Positives.${Antibiotic}.FULL.txt | uniq -c | grep -w "$Antibiotic"r > Pipeline_False_Positives.${Antibiotic}.COUNTS.txt

if [ -s table.queries.${Antibiotic} ]; then
 rm table.queries.${Antibiotic}
fi

if [ -s Variant_ignore.txt ]; then
rm Variant_ignore.txt
fi
   
cat << _EOF_ >  Variant_ignore_Q.txt  
SELECT
Variants_SNP_indel.Gene_name,
Variants_SNP_indel.Variant_annotation
FROM
	Variants_SNP_indel
WHERE
	Variants_SNP_indel.Antibiotic_affected LIKE 'none';
_EOF_

sqlite3 "$RESISTANCE_DB" < Variant_ignore_Q.txt >> Variant_ignore.txt;
sed -i 's/|/ /g' Variant_ignore.txt 



while read line; do 
 #number of hit
 echo "$line" > hit.file
 hitnum=$(awk '{print $1}' hit.file)
 hitdesc=$(awk '{$1=""; print $0}' hit.file)
 awk '{$1=""; print $0}' hit.file > hit.desc.file
 gene_name=$(awk -F"|" '{print $1}' hit.desc.file)
 var_annot=$(awk -F"|" '{print $3}' hit.desc.file)
 echo -e "Here is a potential false positive, $hitdesc. This occurred in $hitnum strains" >> table.queries.${Antibiotic}
 echo -e "If you want to delete it, run" >> table.queries.${Antibiotic}
 echo -e "sqlite3 resistance_database \"DELETE FROM Variants_SNP_indel WHERE Gene_name = '$gene_name' AND Variant_annotation = '$var_annot';\"\n\n" >> table.queries.${Antibiotic}
 rm hit.file hit.desc.file
done < Pipeline_False_Positives.${Antibiotic}.COUNTS.txt
 
while read line; do
 echo "$line" > loss.file
 
 awk '{$1=""; print $0}' loss.file > loss.desc.file
 gene_name=$(awk -F"|" '{print $1}' loss.desc.file)
 rm loss.file loss.desc.file
 while read strain; do
   echo "$strain" >> False.positive.loss.variation.${Antibiotic}.txt
   echo "$gene_name" >> False.positive.loss.variation.${Antibiotic}.txt
    #echo "$gene_name" > gene.name.tmp
    #gene_name=$(cat gene.name.tmp)
    #rm gene.name.tmp
   
    if [ -s False.positive.loss.variation.${Antibiotic}.tmp ]; then
      rm False.positive.loss.variation.${Antibiotic}.tmp
    fi
   
    grep $gene_name ../../Variants/Annotated/${strain}.PASS.indels.annotated.vcf | grep HIGH | awk -F"|" '{ print $4,$11 }' >> False.positive.loss.variation.${Antibiotic}.tmp
	grep $gene_name ../../Variants/Annotated/${strain}.PASS.snps.annotated.vcf | grep HIGH | awk -F"|" '{ print $4,$11 }' >> False.positive.loss.variation.${Antibiotic}.tmp
    sed -i 's/p\.//' False.positive.loss.variation.${Antibiotic}.tmp


    while read f; do 
	  grep -vw "$f" False.positive.loss.variation.${Antibiotic}.tmp > False.positive.loss.variation.${Antibiotic}.tmp.1
	  mv False.positive.loss.variation.${Antibiotic}.tmp.1 False.positive.loss.variation.${Antibiotic}.tmp
	  #rm False.positive.loss.variation.${Antibiotic}.tmp False.positive.loss.variation.${Antibiotic}.tmp.1
    done < Variant_ignore.txt

    cat False.positive.loss.variation.${Antibiotic}.tmp >> False.positive.loss.variation.${Antibiotic}.txt
    echo -e "\n\n" >> False.positive.loss.variation.${Antibiotic}.txt
   
   #reference file check
    if [ ! -s ${ref} ]; then
		echo "This script requires a reference file before the false positive assessment can be run"
		exit 1
	fi
	
    if [ ! -s ${ref}.fai ]; then
		samtools faidx ${ref}
	fi   
	   
	if [ ! -s ${ref}.bed ]; then
		bedtools makewindows -g ${ref}.fai -w 1000 > ${ref}.bed
	fi	
	
	mosdepth --by ${ref}.bed output ../../bams/${strain}.dedup.bam
	sum_depth=$(zcat output.regions.bed.gz | awk '{print $4}' | awk '{s+=$1}END{print s}')
    total_chromosomes=$(zcat output.regions.bed.gz | awk '{print $4}' | wc -l)
    echo "$sum_depth/$total_chromosomes" | bc > ${strain}.depth.txt
	
	echo -e "Chromosome\tStart\tEnd\tInterval" > tmp.header
    zcat output.per-base.bed.gz | awk '$4 ~ /^0/ { print $1,$2,$3,$3-$2 }' > del.summary.tmp
    cat tmp.header del.summary.tmp > ${strain}.FP.deletion_summary.txt
	
	declare -A SQL_loss_report=()
STATEMENT_GENE_LOSS_COV () {
COUNTER=1
while read line; do 
chromosome=$(echo "$line" | awk '{print $1}')
q_start=$(echo "$line" | awk '{print $2}')
q_end=$(echo "$line" | awk '{print $3 }')
SQL_loss_report[$COUNTER]=$(
cat << EOF
SELECT
    Coverage.Chromosome,
    Coverage.Start_coords,
    Coverage.End_coords,
    Coverage.Gene,
    Coverage.Locus_tag,
    Coverage.Upregulation_or_loss,
    Coverage.Antibiotic_affected,
    Coverage.Known_combination,
	Coverage.Comments
FROM
    Coverage
WHERE
    Coverage.Chromosome = '$chromosome'
    AND Coverage.Start_coords < '$q_end'
    AND Coverage.End_coords > '$q_start'
    AND Coverage.Upregulation_or_loss LIKE 'Loss';
EOF
)

COUNTER=$((COUNTER+1))
done < <(tail -n +2 ${strain}.FP.deletion_summary.txt)
LOSS_COUNT="$COUNTER"
}
	STATEMENT_GENE_LOSS_COV
for (( i=1; i<"$LOSS_COUNT"; i++ )); do sqlite3 "$RESISTANCE_DB" "${SQL_loss_report[$i]}" >> ${strain}.AbR_output_del_dup.txt; done

 done < <(awk -F"," '{print $1}' Pipeline_False_Positives.${Antibiotic}.txt)
done < <(grep -w loss Pipeline_False_Positives.${Antibiotic}.COUNTS.txt)

#include logic here to automatically drop line from table if it appears more than X times? Probably have to do the loss manually?
fi


#false negative interrogation
if [ "$false_neg_check" == 1 ]; then
  while read strain; do
	mosdepth --by ${ref}.bed output ../../bams/${strain}.dedup.bam
	sum_depth=$(zcat output.regions.bed.gz | awk '{print $4}' | awk '{s+=$1}END{print s}')
    total_chromosomes=$(zcat output.regions.bed.gz | awk '{print $4}' | wc -l)
    echo "$sum_depth/$total_chromosomes" | bc > ${strain}.depth.txt
	
	echo -e "Chromosome\tStart\tEnd\tInterval" > tmp.header
    zcat output.per-base.bed.gz | awk '$4 ~ /^0/ { print $1,$2,$3,$3-$2 }' > del.summary.tmp
    cat tmp.header del.summary.tmp > ${strain}.FN.deletion_summary.txt
  done < <(awk -F"," '{print $1}' Pipeline_False_Negatives.${Antibiotic}.txt)
fi





echo -e "Summary results for $Antibiotic"
echo -e "Number of false negatives($Antibiotic) = $FalseNegative"
echo -e "Number of false positive($Antibiotic) = $FalsePositive"
echo -e "Number of true positives($Antibiotic) = $TruePositive"
echo -e "Number of true negatives($Antibiotic) = $TrueNegative"
echo -e "Number of true intermediates($Antibiotic) = $TrueIntermediate"
echo -e "Number of false intermediates($Antibiotic) = $FalseIntermediate"
echo -e "            Resistant strains identified as intermediate = $FalseIntermediate_res_FN" 
echo -e "            Sensitive strains identified as intermediate = $FalseIntermediate_sens_FP" 
echo -e "            Intermediate strains identified as resistant = $FalseIntermediate_int_FP" 
echo -e "            Intermediate strains identified as sensitive = $FalseIntermediate_int_FN" 
TotalStrains=$(echo "$FalseNegative + $FalsePositive + $TruePositive + $TrueNegative + $TrueIntermediate + $FalseIntermediate" | bc)
echo -e "Total number of strains($Antibiotic) = $TotalStrains\n"
exit 0
