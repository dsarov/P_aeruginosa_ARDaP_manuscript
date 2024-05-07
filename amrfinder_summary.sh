#! /bin/bash

#echo -e "This script will determine the rate of false positives, false negatives, true positive and true negatives following a complete ARDaP run.\n\n"

#to run this script
# e.g. amrfinder_summary.sh MEROPENEM CARBAPENEM

#The first variable is the anitibiotic to be interrogated (and the files containing the list of sensitive strains and resistant strains should follow the convention MEROPENEM.R.txt and MEROPENEM.S.txt in the above example.
#The second variable is the search term used to look through the AMRFinder outputs.


false_pos_check=0
false_neg_check=0


Antibiotic=$1

search_term=$2
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
    if [ ! -s ${SensitiveStrain}.output ]; then
      echo -e "Cannot find the AbR output for the sensitive strain ${SensitiveStrain}. Did resfinder complete normally?"
      echo -e "Please check and re-run"
	  echo "${SensitiveStrain}" >> missing_isolates.txt
      #sleep 5
    else 
    #  echo -e "Testing sensitive strain $SensitiveStrain\n"
	  #find res line
	  if [ -s tmp.txt ]; then 
	    rm tmp.txt
	  fi
	  #grep -w ${search_term} ${SensitiveStrain}.output &> /dev/null
	  awk -v search=${search_term} -F "\t" '$12 ~ search { exit 1 }' ${SensitiveStrain}.output &> /dev/null
	  status=$?
	  if [ "$status" == 1 ]; then #sensitive strain but called as resistant
	     FalsePositive=$((FalsePositive+1))
		 
	     echo -e "$SensitiveStrain,False Positive for $Antibiotic" >> Pipeline_False_Positives.${Antibiotic}.txt
	  else 
	     TrueNegative=$((TrueNegative+1)) #Sensitive strain not called as resistant or intermediate
         echo -e "$SensitiveStrain,True Negative for $Antibiotic" >> Pipeline_True_Negatives.${Antibiotic}.txt		 
	  fi
	fi
done < "$Antibiotic".S.txt


if [ -s "$Antibiotic".R.txt ]; then
  while read ResistantStrain; do
    if [ ! -s ${ResistantStrain}.output  ]; then
	  echo -e "Cannot find AbR output for the resistant strain ${ResistantStrain}"
      echo -e "Please check and re-run"
	  echo "${ResistantStrain}" >> missing_isolates.txt
      #sleep 5
	else
    # echo -e "Testing resistant strain $ResistantStrain\n"
      if [ -s tmp.txt ]; then 
	    rm tmp.txt
	  fi	 
	  #grep -w ${search_term} ${ResistantStrain}.output &> /dev/null
	  awk -v search=${search_term} -F "\t" '$12 ~ search { exit 1 }' ${ResistantStrain}.output &> /dev/null
	  status=$? 
	  if [ "$status" == 1 ]; then  #Resistant strain called as resistant
	   TruePositive=$((TruePositive+1))
	   echo -e "$ResistantStrain,True Positive for $Antibiotic" >> Pipeline_True_Positives.${Antibiotic}.txt
	  else
	   FalseNegative=$((FalseNegative+1)) #Resistant strain not called as resistant	
	   echo -e "$ResistantStrain,False Negative for $Antibiotic" >> Pipeline_False_Negatives.${Antibiotic}.txt
      fi	  
	fi 

  done < "$Antibiotic".R.txt
fi





echo -e "Summary results for $Antibiotic"
echo -e "Number of false negatives($Antibiotic) = $FalseNegative"
echo -e "Number of false positive($Antibiotic) = $FalsePositive"
echo -e "Number of true positives($Antibiotic) = $TruePositive"
echo -e "Number of true negatives($Antibiotic) = $TrueNegative"
#echo -e "Number of true intermediates($Antibiotic) = $TrueIntermediate"
#echo -e "Number of false intermediates($Antibiotic) = $FalseIntermediate"
#echo -e "            Resistant strains identified as intermediate = $FalseIntermediate_res_FN" 
#echo -e "            Sensitive strains identified as intermediate = $FalseIntermediate_sens_FP" 
#echo -e "            Intermediate strains identified as resistant = $FalseIntermediate_int_FP" 
#echo -e "            Intermediate strains identified as sensitive = $FalseIntermediate_int_FN" 
TotalStrains=$(echo "$FalseNegative + $FalsePositive + $TruePositive + $TrueNegative" | bc)
echo -e "Total number of strains($Antibiotic) = $TotalStrains\n"
exit 0
