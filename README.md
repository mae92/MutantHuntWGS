# MutantHuntWGS


## A bioinformatics pipeline for identification and characterization of mutations in *Saccharomyces cerevisiae*. MutantHuntWGS compares data, input in FASTQ format, from a mutant strain to a wild-type strain to identify high confidence sequence variants present only in the mutant. This pipeline was designed to be as user friendly as possible.


![picture alt](https://github.com/mae92/MutantHuntWGS/blob/master/Figure_1-01.jpg)


# Setup and Run

1. Follow the instructions at https://docs.docker.com/get-docker/ to download and install Docker on your computer.

2. Create a directory (a folder) in a directory of your choice (on your Desktop is fine) named `./AnalysisDirectory`. Within `./AnalysisDirectory` create a file directory named `./FASTQ` and place all FASTQ files into it (full path: `./PATH_TO_DESKTOP/Analysis_Directory/FASTQ`). 

    1. Ensure that FASTQ files adhere to the naming convention described below. **THIS IS REALLY IMPORTANT.**
    
        1. Single end sequencing FASTQ file should be named: FILENAME.fastq.
        
        2. Paired end sequencing FASTQ files should be named: FILENAME_R1.fastq and FILENAME_R2.fastq.
        
        3. I suggest making copies of your FASTQ files rather then renaming the originals in case a mistake is made during the renaming process.
        
        4. For this naming example "FILENAME" will be used as the input for the -n option below. 
        
        5. "FILENAME" should not have any spaces or punctuation, not even underscores.

3. Open the Terminal and download and run the Docker container for MutantHuntWGS by copying and pasting the following command into the terminal: 

```
docker run -it -v /PATH_TO_DESKTOP/Analysis_Directory:/Main/Analysis_Directory mellison/mutant_hunt_wgs:version1
```

4. Run MutantHuntWGS by running the code below.
```
MutantHuntWGS.sh \
    -n FILENAME \
    -g /Main/MutantHuntWGS/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
    -f /Main/MutantHuntWGS/S_cerevisiae_Bowtie2_Index_and_FASTA/genome.fa \
    -r single \
    -s 100 \
    -p /Main/MutantHuntWGS/S_cerevisiae_Bowtie2_Index_and_FASTA/ploidy_n1.txt \
    -d /Main/Analysis_Directory/FASTQ \
    -o /Main/Analysis_Directory/NAME_YOUR_OUTPUT_FOLDER \
    -a YES
```

## Because the files and directory structure were set up ahead of time (on the Desktop and during the docker build), and this all runs out of the Docker container, *the file paths in the above commands will all stay the same*, but some of the options may change depending upon your needs. 



# MutantHuntWGS Options

## All options are required

### -n 

The -n option takes the prefix of the FASTQ file name for the wild-type strain. For the example of FILENAME.fastq or FILENAME_R1.fastq this prefix would simply be "FILENAME".

### -g

The -g option takes the file PATH to the bowtie index files and the file prefix (genome). **Use exactly what is shown above for this command.**

### -f

The -f option takes the file PATH and file name of the genome FASTA file (genome.fa) **Use exactly what is shown above for this command.**

### -r

The -r option specifies whether the input data contains paired-end or single-end reads and can take values of "paired" or "single".

### -s

The -s option takes a score cutoff for the variant scores. This score is calculated by the following formula: -10 * log10(P) where P is the probability that the variant call (ALT) in the VCF file is wrong. 

So a score of:

10 means a P of 0.1 and a 10% chance the ALT is wrong

20 means a P of 0.01 and a 1% chance the ALT is wrong

30 means a P of 0.001 and a 0.1% chance the ALT is wrong

40 means a P of 0.0001 and a 0.01% chance the ALT is wrong

50 means a P of 0.00001 and a 0.001% chance the ALT is wrong

100 means a P of 0.000000001 and a 0.0000001% chance the ALT is wrong


### -p

The -p option takes the file PATH and file name of the ploidy file (genome.fa) **Use exactly what is shown above for this command.**

Ploidy files are available for haploid (ploidy_n1.txt) and diploid (ploidy_n2.txt) and are in the following format:

chr start   end sex ploidy

ploidy_n1.txt:
```
chrI	1	230218	M	1
chrII	1	813184	M	1
chrIII	1	316620	M	1
```
ploidy_n2.txt
```
chrI	1	230218	M	2
chrII	1	813184	M	2
chrIII	1	316620	M	2
```
so if you wanted to account for a ploidy of 2 on chromosome II in an otherwise haploid yeast strain you should be able to edit the file like this:
```
chrI	1	230218	M	1
chrII	1	813184	M	2
chrIII	1	316620	M	1
```
A sex of M or male was arbitrarily chosen and the MutantHuntWGS program is expecting that so it cannot be changed without editing MutantHuntWGS.sh.

### -d

Directory containing your FASTQ files. If you set things up in the way that the instructions outline above this should stay the same as the example: `/PATH_TO_DESKTOP/Analysis_Directory/FASTQ`. **Use exactly what is shown above for this command.**

### -o

This allows you to specify a folder for your data to output to. This should be structured like the example `/Main/Analysis_Directory/NAME_YOUR_OUTPUT_FOLDER` except you will come up with a descriptive name to replace the `NAME_YOUR_OUTPUT_FOLDER` part of the file PATH.

### -a

This allows you to turn on and off the alignment and calling step. So if you have already aligned reads and called variants and all that you want to do is reanalyze with a different score cuttoff then you can set this to "NO", but if you are starting from FASTQ files that have not gone throught this process yet you set this to "YES"
