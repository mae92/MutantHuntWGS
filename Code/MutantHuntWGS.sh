#!/bin/bash

#MutantHunter script by Mitch Ellison


#Enable user input

while getopts n:g:f:r:s:p:d:o:a: option
do
case "${option}"
in
n) WILD_TYPE=${OPTARG};;
g) GENOME=${OPTARG};;
f) GENOME_FASTA=${OPTARG};;
r) READ_TYPE=${OPTARG};;
s) SCORE=${OPTARG};;
p) PLOIDY_FILE=${OPTARG};;
d) ANALYSIS_DIRECTORY=${OPTARG};;
o) OUTPUT_FILE=${OPTARG};;
a) ALIGNMENT_AND_CALLING=${OPTARG};;
esac
done

# echo "$FASTQ"
# echo "$GENOME"
# echo "$GENOME_FASTA"
# echo "$READ_TYPE"
# echo "$SCORE"
# echo "$OUTPUT_FILE"
# echo "$ANALYSIS_DIRECTORY"


#Remove Old Output directory
rm -rf "$OUTPUT_FILE"

#Make Output directory
mkdir "$OUTPUT_FILE"


#make some directories
mkdir "$OUTPUT_FILE"/BAM
mkdir "$OUTPUT_FILE"/Alignment_Stats
mkdir "$OUTPUT_FILE"/BCF
mkdir "$OUTPUT_FILE"/VCF
mkdir "$OUTPUT_FILE"/VCFtools_Output
mkdir "$OUTPUT_FILE"/SNPeff_Output
mkdir "$OUTPUT_FILE"/SIFT_Output


if ["$ALIGNMENT_AND_CALLING" = "YES" ]
then

	echo -e "\n\n" Running MutantHuntWGS 

	###################################
	###     MODULE 1: Alignment     ###
	###################################

	echo -e "\n\n" Aligning Reads to the Genome "\n"

	#Align all datasets in a single module up-front as Jen Walker suggested

	#Apply the -r option to signify if the input is paired-end or single-end end reads

	if [ "$READ_TYPE" = "paired" ]
	then

	#Alignment of all FASTQ Files if -r paired is set

	#Run Bowtie 2 and convert output to BAM with Samtools

	##Loop to process all files starts here
	##map FASTQ files using for-loop


	for FASTQ_FILE in `ls "$ANALYSIS_DIRECTORY"/FASTQ/*_R1.fastq.gz`
	do

		##Align reads

		#Using this command to get only the file name and lose the path

		NAME_PREFIX=`echo "$FASTQ_FILE" | awk -F "/" '{print $(NF)}' | awk -F "_" '{print  $1}'`

		bowtie2 --no-unal \
			-q \
			-k 2 \
			-p 16 \
			-x "$GENOME" \
			-1 "$ANALYSIS_DIRECTORY"/FASTQ/"$NAME_PREFIX"_R1.fastq.gz -2 "$ANALYSIS_DIRECTORY"/FASTQ/"$NAME_PREFIX"_R2.fastq.gz \
			2> "$OUTPUT_FILE"/Alignment_Stats/"$NAME_PREFIX"_bowtie2.txt \
			| samtools view - -bS -q 30 > "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_unsorted.bam

		##Sort BAM with Samtools

		samtools sort "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_unsorted.bam -o "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam

		##Index BAM with Samtools

		samtools index "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam

		echo  -e Reads from "$ANALYSIS_DIRECTORY"/FASTQ/"$NAME_PREFIX"*.fastq.gz Mapped "\n"and stored in "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam"\n"

		done


	elif [ "$READ_TYPE" = "single" ]
	then


		for FASTQ_FILE in `ls "$ANALYSIS_DIRECTORY"/FASTQ/*.fastq.gz`
		do

		#Using this command to get only the file name and lose the path

		NAME_PREFIX=`echo "$FASTQ_FILE" | awk -F "/" '{print $(NF)}' | awk -F "." '{print  $1}'`

		##Align reads

		bowtie2 --no-unal \
			-q \
			-k 2 \
			-p 16 \
			-x "$GENOME" \
			-U "$ANALYSIS_DIRECTORY"/FASTQ/"$NAME_PREFIX".fastq.gz \
			2> "$OUTPUT_FILE"/Alignment_Stats/"$NAME_PREFIX"_bowtie2.txt \
			| samtools view - -bS -q 30 > "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_unsorted.bam

		##Sort BAM with Samtools

		samtools sort "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_unsorted.bam -o "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam

		##Index BAM with Samtools

		samtools index "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam

		echo -e Reads from "$ANALYSIS_DIRECTORY"/FASTQ/"$NAME_PREFIX"*.fastq.gz mapped"\n"and stored in "$OUTPUT_FILE"/BAM/"$NAME_PREFIX"_sorted.bam"\n"

		done

	else

		echo "Failed to Identify FASTQ Files"

	fi


	###################################
	###  MODULE 2: Variant Calling  ###
	###################################

	echo -e "\n\n" Calling Variants "\n"

	{

	for BAM_FILE in `ls "$OUTPUT_FILE"/BAM/*_sorted.bam`
	do


		#Using this command to get only the file name and lose the path
		NAME_PREFIX=`echo "$BAM_FILE" | awk -F "/" '{print $(NF)}' | awk -F "_" '{print  $1}'`

		#The first step is to use the SAMtools mpileup command to calculate the genotype likelihoods supported by the aligned reads in our sample
	
		{
		samtools mpileup -g -f "$GENOME_FASTA" "$BAM_FILE" -o "$OUTPUT_FILE"/BCF/"$NAME_PREFIX"_variants.bcf
		} &> /dev/null
	
		#-g: directs SAMtools to output genotype likelihoods in the binary call format (BCF). This is a compressed binary format.
		#-f: directs SAMtools to use the specified reference genome. A reference genome must be specified.

		#The second step is to use bcftools:
		#The bcftools call command uses the genotype likelihoods generated from the previous step to call SNPs and indels, and outputs the all identified variants in the variant call format (VFC), the file format created for the 1000 Genomes Project, and now widely used to represent genomic variants.

		echo "$BAM_FILE" | awk '{print $1 "\t" "M"}' > "$OUTPUT_FILE"/BCF/sample_file.txt

		bcftools call -c -v --samples-file "$OUTPUT_FILE"/BCF/sample_file.txt --ploidy-file "$PLOIDY_FILE" "$OUTPUT_FILE"/BCF/"$NAME_PREFIX"_variants.bcf > "$OUTPUT_FILE"/VCF/"$NAME_PREFIX"_variants.vcf

		echo -e Variants have been called and stored in "\n""$OUTPUT_FILE"/VCF/"$NAME_PREFIX"_variants.vcf"\n"

	done

fi

#############################################
###  MODULE 3: Filter and Score Variants  ###
#############################################

echo -e "\n\n" Comparing all strains to the "$WILD_TYPE" wild-type strain "\n"

#WILD_TYPE = the wild type file to compare all other files too

for MUTANT_VCF in `ls "$OUTPUT_FILE"/VCF/*_variants.vcf`
do

	#Using this command to get only the file name and lose the path
	MUTANT=`echo "$MUTANT_VCF" | awk -F "/" '{print $(NF)}' | awk -F "." '{print  $1}'`


	#Compare variations in the two files

	vcftools --vcf "$OUTPUT_FILE"/VCF/"$MUTANT".vcf --diff "$OUTPUT_FILE"/VCF/"$WILD_TYPE"_variants.vcf --diff-site --out "$OUTPUT_FILE"/VCFtools_Output/"$MUTANT"_vs_"$WILD_TYPE"_VCFtools_output.txt

	##These next two steps are to produce VCF files from the output for downstream applications

	##Here I am using awk to remove the third column so that I can make files for use with the grep command below

	awk -F "\t" '$2 != $3' "$OUTPUT_FILE"/VCFtools_Output/"$MUTANT"_vs_"$WILD_TYPE"_VCFtools_output.txt.diff.sites_in_files | awk -F "\t" '$3 == "." {print $1 "\t" $2}' > "$OUTPUT_FILE"/VCFtools_Output/"$MUTANT"_vs_"$WILD_TYPE"_for_grep.txt


	#VCF file construction

	##Grep variants of intrest from Mutant VCF

	grep -f "$OUTPUT_FILE"/VCFtools_Output/"$MUTANT"_vs_"$WILD_TYPE"_for_grep.txt "$OUTPUT_FILE"/VCF/"$MUTANT".vcf > "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered.tmp


	#Get VCF header information so we can add it to the new file allowing it to be viewed in IGV
	grep "^#" "$OUTPUT_FILE"/VCF/"$MUTANT".vcf > "$OUTPUT_FILE"/VCF/VCF_header.tmp

	#Score variants based on user input
	awk ' $6 >= '$SCORE' {print}' "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered.tmp > "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered_and_scored.tmp

	#Add header to filtered variants to generate a proper VCF file
	cat "$OUTPUT_FILE"/VCF/VCF_header.tmp "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered.tmp > "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered.vcf

	#Add header to filtered and scored variants to generate a proper VCF file
	cat "$OUTPUT_FILE"/VCF/VCF_header.tmp "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered_and_scored.tmp > "$OUTPUT_FILE"/VCF/"$MUTANT"_filtered_and_scored.vcf


done

echo -e Comparisons complete


#Remove unwanted directories
rm -rf "$OUTPUT_FILE"/BCF
rm -rf "$OUTPUT_FILE"/VCFtools_Output

#Remove unwanted files
rm "$OUTPUT_FILE"/BAM/*_unsorted.bam
rm "$OUTPUT_FILE"/VCF/*.tmp
rm `ls "$OUTPUT_FILE"/VCF/"$WILD_TYPE"_variants_filtered*.vcf`




###########################################################
###   MODULE 4: Determine Variant Effects with SNPeff   ###
###########################################################

echo -e "\n\n" Running SNPeff  "\n"


cd "$OUTPUT_FILE"/SNPeff_Output

for VCF_FILE in `ls "$OUTPUT_FILE"/VCF/*_filtered_and_scored.vcf`
do

	#Using this command to get only the file name and lose the path
	VCF_NAME=`echo "$VCF_FILE" | awk -F "/" '{print $(NF)}' | awk -F "_" '{print  $1}'`
	
	#Make directory
	mkdir "$OUTPUT_FILE"/SNPeff_Output/"$VCF_NAME"
	
	#set as working directory
	cd "$OUTPUT_FILE"/SNPeff_Output/"$VCF_NAME"

	#Use SNPeff to find out some more information about these SNPs like the amino acid change in the resulting protein
	{
	java -Xmx4G -jar /Main/snpEff/snpEff.jar ann -v Saccharomyces_cerevisiae "$VCF_FILE" > "$OUTPUT_FILE"/SNPeff_Output/"$VCF_NAME"/SNPeff_Annotations.vcf
	} &> /dev/null
done


#################################################################
###   MODULE 5: Predict Variant Impact on Protein with SIFT   ###
#################################################################

echo -e "\n\n" Running SIFT  "\n"

for VCF_FILE in `ls "$OUTPUT_FILE"/VCF/*_filtered_and_scored.vcf`
do

	#Using this command to get only the file name and lose the path
	VCF_NAME=`echo "$VCF_FILE" | awk -F "/" '{print $(NF)}' | awk -F "_" '{print  $1}'`

	##Run SIFT
	{
	java -Xmx4G -jar /Main/SIFT4G_Annotator.jar -c -i "$VCF_FILE"  -d /Main/EF4.74 -r "$OUTPUT_FILE"/SIFT_Output/"$VCF_NAME"
	} &> /dev/null
done

echo -e "\n\n" Run Completed "\n\n\n\n"

exit
