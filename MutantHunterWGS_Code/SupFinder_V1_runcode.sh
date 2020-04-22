#Set working directory
cd /MutantHunter/Analysis_Directory/MutantHunter-master/MutantHunter_Code

#modify script to be executable
chmod 777 /MutantHunter/Analysis_Directory/MutantHunter-master/MutantHunter_Code/Supfinder_V1.1.sh

#Test run on paired end data
./Supfinder_V1.1.sh \
-m /MutantHunter/Analysis_Directory/My_Analysis/Mutant_FASTQ/mut_toy \
-n /MutantHunter/Analysis_Directory/My_Analysis/WT_FASTQ/wt_toy \
-g /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-f /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-r paired \
-s 100 \
-o /MutantHunter/Analysis_Directory/My_Analysis/MutantHunter_paired_end_test

#Test run on single end data
./Supfinder_V1.1.sh \
-m /MutantHunter/Analysis_Directory/My_Analysis/Mutant_FASTQ/mut_toy \
-n /MutantHunter/Analysis_Directory/My_Analysis/WT_FASTQ/wt_toy \
-g /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-f /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-r single \
-s 100 \
-o MutantHunter_single_end_test


#Test run on single end data
./Supfinder_V1.1.gz.sh \
-m /MutantHunter/Analysis_Directory/SUPPRESS22 \
-n /MutantHunter/Analysis_Directory/Analysis_Directory/WT \
-g /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-f /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
-r single \
-s 100 \
-o /MutantHunter/Analysis_Directory/MutantHunter_Sup22_vs_Wt


