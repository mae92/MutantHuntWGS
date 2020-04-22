#MutantHunter Runcode

cd /MutantHunter/Analysis_Directory/MutantHunter-master/MutantHunter_Code

chmod 777 MutantHunter.sh

MutantHunter.sh \
-n WTonemillion \
-g /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-f /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/genome.fa \
-r single \
-s 10 \
-p /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/ploidy_n1.txt \
-d /MutantHunter/Analysis_Directory \
-o /MutantHunter/Analysis_Directory/MutantHunter_Output_for_SUP22_vs_WT_one_million_reads_each_2


