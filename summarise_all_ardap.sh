

if [ -s ~/ARDaP_full_summary.txt ]; then
rm ~/ARDaP_full_summary.txt
fi


cd /home/dsarovich/analyses/pa_ardap/barrio-tofino/Outputs/AbR_reports
for i in AMK CAZ CIP CST FEP IPM TZP MEM TOB; 
do
echo -e "Summary for Barrio-toffino $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/belkum/Outputs/AbR_reports
for i in AMK MEM CST FEP PIP PMB; 
do
echo -e "Summary for Belkum $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/buhl/Outputs/AbR_reports
for i in AMK MEM CST CIP FEP TZP CAZ PIP; 
do
echo -e "Summary for Buhl $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done

cd /home/dsarovich/analyses/pa_ardap/ramanathan/Outputs/AbR_reports
for i in AMK MEM CST CIP FEP TZP CAZ IPM; 
do
echo -e "Summary for ramanathan $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done

cd /home/dsarovich/analyses/pa_ardap/cabot/Outputs/AbR_reports
for i in AMK MEM CAZ CIP CST FEP TZP IPM TOB; 
do
echo -e "Summary for Cabot $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/cdc/Outputs/AbR_reports
for i in AMK FEP CAZ CZT MEM CIP CST IPM TZP TOB; 
do
echo -e "Summary for cdc $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/khaledi/Outputs/AbR_reports
for i in CAZ MEM CIP TOB; 
do
echo -e "Summary for khaledi $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/kos/Outputs/AbR_reports
for i in MEM AMK CST; 
do
echo -e "Summary for kos $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done

cd /home/dsarovich/analyses/pa_ardap/schi/Outputs/AbR_reports
for i in AMK TOB PIP TZP CAZ FEP CIP CST PMB IPM MEM; 
do
echo -e "Summary for schi $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done

cd /home/dsarovich/analyses/pa_ardap/sherrard/Outputs/AbR_reports
for i in AMK MEM CST PMB CIP FEP; 
do
echo -e "Summary for sherrard $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/tsang/Outputs/AbR_reports
for i in AMK MEM CIP TZP CAZ TOB; 
do
echo -e "Summary for Tsang $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done


cd /home/dsarovich/analyses/pa_ardap/wardell/Outputs/AbR_reports
for i in CIP MEM IPM TOB; 
do
echo -e "Summary for wardell $i\n" >> ~/ARDaP_full_summary.txt
~/bin/ardap_summary_v2.5_just_numbers.sh $i >> ~/ARDaP_full_summary.txt
done
