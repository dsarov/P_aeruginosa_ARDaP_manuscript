antibiotic=$1

#barrio-tofino buhl belkum cabot cdc khaledi kos ramanathan schi sherrard tsang wardell


if [ $antibiotic = ALL ]; then
for i in barrio-tofino buhl belkum cabot cdc khaledi kos ramanathan schi sherrard tsang wardell; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  for f in *.R.txt; do
    qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh ${f/.R.txt/}" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
    #qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
    #depend="-W depend=afterok${qsub_cat_ids}"
    #qsub "$depend" -v command="echo $i >> ~/Summary_${f/.R.txt/}.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh ${f/.R.txt/} >> ~/Summary_${f/.R.txt/}.txt" ~/bin/Header.pbs
  done
done

fi

if [ $antibiotic = TOB ]; then
for i in barrio-tofino cabot khaledi tsang cdc wardell; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
  if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi	
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi


if [ $antibiotic = AMK ]; then
for i in barrio-tofino buhl belkum sherrard cabot tsang cdc ramanathan kos; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
  if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = PIP ]; then
for i in buhl belkum; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = TZP ]; then
for i in tsang barrio-tofino cabot cdc ramanathan buhl; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = CAZ ]; then
for i in buhl barrio-tofino ramanathan tsang cabot khaledi cdc; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = FEP ]; then
for i in buhl barrio-tofino belkum cdc cabot sherrard; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = CIP ]; then
for i in barrio-tofino ramanathan tsang buhl cdc cabot khaledi wardell sherrard; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = CST ]; then
for i in barrio-tofino belkum ramanathan buhl cdc kos cabot sherrard; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = PMB ]; then
for i in belkum sherrard; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = IPM ]; then
for i in cdc cabot barrio-tofino wardell; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt;  ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi

if [ $antibiotic = MEM ]; then
for i in barrio-tofino tsang ramanathan buhl cdc kos cabot khaledi wardell belkum sherrard; do
  cd /home/dsarovich/analyses/pa_ardap/"$i"
  pwd
  rm -r ./tmp
  rm ./Outputs/CARD/*backup*
  rm ./Outputs/AbR_reports/*backup*
  rm ./Header.*
  rm -r ./Outputs/AbR_reports/summary*"$antibiotic"
  rm qsub_id.txt
  for f in *_1_*; do 
    qsub_id=$(qsub -v command="~/bin/re_analyse_ardap.sh ${f}" ~/bin/Header.pbs); echo "$qsub_id" >> qsub_id.txt; 
  done

  qsub_cat_ids=$(cat qsub_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n')
  depend="-W depend=afterok${qsub_cat_ids}"
  cd /home/dsarovich/analyses/pa_ardap/"$i"/Outputs/AbR_reports/
  rm qsub_id2.txt
  qsub_id=$(qsub "$depend" -v command="~/bin/ardap_summary_v2.5_command_line.sh $antibiotic" ~/bin/Header.pbs); echo "$qsub_id" > qsub_id2.txt
  qsub_cat_ids=$(cat qsub_id2.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n' )
  depend="-W depend=afterok${qsub_cat_ids}"
    if [ -s ~/Summary_$antibiotic.txt ]; then 
    rm ~/Summary_$antibiotic.txt
  fi
  qsub "$depend" -v command="echo $i >> ~/Summary_$antibiotic.txt; ~/bin/ardap_summary_v2.5_just_numbers.sh $antibiotic >> ~/Summary_$antibiotic.txt" ~/bin/Header.pbs
  
done

fi
#for i in barrio-tofino belkum buhl cabot cdc khaledi kos ramanathan schi sherrard tsang wardell; do