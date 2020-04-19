#MutantHunter Runcode

cd /MutantHunter/Analysis_Directory/MutantHunter-master/MutantHunter_Code

chmod 777 MutantHunter.sh

./MutantHunter.sh \
-n WTonemillion \
-g /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-f /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome.fa \
-r single \
-s 10 \
-p /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/ploidy_n1.txt \
-d /MutantHunter/Analysis_Directory/My_Analysis \
-o /MutantHunter/Analysis_Directory/MutantHunter_Output_for_SUP22_vs_WT_one_million_reads_each


