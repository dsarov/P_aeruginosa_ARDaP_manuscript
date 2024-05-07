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
good_int=0
bad_int=0
## Logic flow

# for each antibiotic; determine the following numbers

# false negative is a strain that is in the resistant list but called as sensitive

# false positive is a strain that is in the sensitive list but called as resistant

# true negative is a strain that is in the sensitive list and called as sensitive

# true positive is a strain that is in the resistant list and called as resistant

# Load the list of sensitive strains for each antibiotic



if [ ! -f "$Antibiotic".R.txt ]; then
   echo "Cannot find file listing resistant strains"
   echo "I am looking for this file --> ${Antibiotic}.R.txt"
   echo "Please check and re-run if there should be resistant strains in this dataset"
  # exit 1
fi
if [ ! -f "$Antibiotic".S.txt ]; then
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
    if [ ! -s ../${SensitiveStrain}.AbR_output.cumulative.final.txt ]; then
      echo -e "Cannot find the AbR output for the sensitive strain ${SensitiveStrain}. Did ARDaP complete normally?"
      echo -e "Please check and re-run"
      echo "Exiting"
      exit 1
    else 
      # echo -e "Testing $SensitiveStrain\n"
	#  Antibiotic_res=${Antibiotic}r
      res_score=$(grep "$Antibiotic" ../${SensitiveStrain}.AbR_output.cumulative.final.txt  | awk -F" " '{print $7}')
	  if [[ "$res_score" -gt 99 ]]; then #sensitive strain but called as resistant
	     FalsePositive=$((FalsePositive+1))
		 echo "FP for ${SensitiveStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
	  else #look for intermediate
        if [[ "$res_score" -gt 45 ]]; then #sensitive strain but called as intermediate
           FalseIntermediate=$((FalseIntermediate+1))
           FalseIntermediate_sens_FP=$((FalseIntermediate_sens_FP+1))  #Sensitive strain called as Intermediate	
		   echo "S_as_I_FP for ${SensitiveStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
	    else 
		   TrueNegative=$((TrueNegative+1))
           echo "TN for ${SensitiveStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log		   #Sensitive strain not called as resistant or intermediate
 
		fi
	  fi
    fi
done < ../"$Antibiotic".S.txt


if [ -f ../"$Antibiotic".R.txt ]; then
  while read ResistantStrain; do
    if [ ! -s ../${ResistantStrain}.AbR_output.cumulative.final.txt ]; then
	  echo -e "I Cannot find AbR output for the resistant strain ${ResistantStrain}"
      echo -e "Please check and re-run"
      echo "Exiting"
      exit 1
	else
     # echo -e "Testing $ResistantStrain\n"	
      res_score=$(grep "$Antibiotic" ../${ResistantStrain}.AbR_output.cumulative.final.txt | awk -F" " '{print $7}')
	   if [[ "$res_score" -gt 99 ]]; then #Resistant strain called as resistant
	     TruePositive=$((TruePositive+1))
		 echo "TP for ${ResistantStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
	   else
        if [[ "$res_score" -gt 45 ]]; then 
	     FalseIntermediate=$((FalseIntermediate+1)) #resistant strain identified as intermediate
		 FalseIntermediate_res_FN=$((FalseIntermediate_res_FN+1))
		 echo "R_as_I_FN for ${ResistantStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
       else
	     FalseNegative=$((FalseNegative+1)) #Resistant strain not called as resistant	
		 echo "FN for ${ResistantStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
       fi	  
	  fi 
	fi  
  done < ../"$Antibiotic".R.txt
fi

if [ -f ../"$Antibiotic".I.txt ]; then
    echo "found list of intermediate strains"
    while read IntermediateStrain; do
        if [ ! -s ../${IntermediateStrain}.AbR_output.cumulative.final.txt ]; then
	        echo -e "Cannot find AbR output for the Intermediate strain ${IntermediateStrain}"
            echo -e "Please check and re-run"
            echo "Exiting"
            exit 1
	    else
            # echo -e "Testing $IntermediateStrain\n"	
            res_score=$(grep "$Antibiotic" ../${IntermediateStrain}.AbR_output.cumulative.final.txt  | awk -F" " '{print $7}')
	        if [[ "$res_score" -gt 99 ]]; then #Intermediate strain called as resistant
	            FalseIntermediate=$((FalseIntermediate+1)) #Intermediate strain identified as resistant	
				FalseIntermediate_int_FP=$((FalseIntermediate_int_FP+1))
				echo "I_as_R_FP for ${IntermediateStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
				#echo "Found $IntermediateStrain" 
	        else
                if [[ "$res_score" -gt 45 ]]; then   #Intermediate strain called as intermediate
				  TrueIntermediate=$((TrueIntermediate+1))
				  echo "I_true for ${IntermediateStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
				else
				  FalseIntermediate=$((FalseIntermediate+1))
				  FalseIntermediate_int_FN=$((FalseIntermediate_int_FN+1)) 
				  echo "I_as_S_FN for ${IntermediateStrain} and $Antibiotic. Score = $res_score" >> ~/Int_check.log
				fi
	        fi 
	    fi  
    done < ../"$Antibiotic".I.txt
fi


echo -e "Summary results for $Antibiotic"
echo -e "Number of false negatives($Antibiotic) = $FalseNegative"
echo -e "Number of false positive($Antibiotic) = $FalsePositive"
echo -e "Number of false intermediates($Antibiotic) = $FalseIntermediate"
echo -e "Number of true positives($Antibiotic) = $TruePositive"
echo -e "Number of true negatives($Antibiotic) = $TrueNegative"
echo -e "Number of true intermediates($Antibiotic) = $TrueIntermediate"

echo -e "            Resistant strains identified as intermediate($Antibiotic) = $FalseIntermediate_res_FN" 
echo -e "            Sensitive strains identified as intermediate($Antibiotic) = $FalseIntermediate_sens_FP" 
echo -e "            Intermediate strains identified as resistant($Antibiotic) = $FalseIntermediate_int_FP" 
echo -e "            Intermediate strains identified as sensitive($Antibiotic) = $FalseIntermediate_int_FN"
good_int=$(echo "$FalseIntermediate_int_FP+$FalseIntermediate_res_FN" | bc)
bad_int=$(echo "$FalseIntermediate_int_FN+$FalseIntermediate_sens_FP" | bc)
#echo -e "            Bad intermediate($Antibiotic) = $bad_int"
#echo -e "            Good intermediate($Antibiotic) = $good_int"           
total_strains_no_int=$(echo "$FalseNegative + $FalsePositive + $TruePositive + $TrueNegative + $FalseIntermediate_res_FN + $FalseIntermediate_sens_FP" | bc )
TotalStrains=$(echo "$FalseNegative + $FalsePositive + $TruePositive + $TrueNegative + $TrueIntermediate + $FalseIntermediate" | bc)
echo -e "Total number of strains($Antibiotic) = $TotalStrains"
echo -e "Total number of strains($Antibiotic) without intermediate = $total_strains_no_int\n"
exit 0



if [ "$res_score" -gt 99 ]; then 
   echo "res_score is gt 99"
else
   echo "res_score is lt 99"
   if [ $res_score -lt 0 ]; then
     echo "res_score lt 0"
   fi
fi
 