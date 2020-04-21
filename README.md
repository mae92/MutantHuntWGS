# MutantHunter


#### A bioinformatics pipeline for identification and characterization of mutations in Saccharomyces cerevisiae from FASTQ files from a wild-type strain and one or more mutant strains.


# Simple Setup / Quick Start

1. Follow the instructions at https://docs.docker.com/get-docker/ to download and install Docker on your computer.

2. Create a directory (a folder) in a directory of your choice (on your Desktop is fine) named `/AnalysisDirectory`. Within `/AnalysisDirectory` create a file directory named `/FASTQ` and place all FASTQ files into it (full path: `/PATH_TO_DESKTOP/Analysis_Directory/FASTQ`). 

    1. Ensure that FASTQ files adhere to the naming convention described below. THIS IS REALLY IMPORTANT.
    
        1. Single end sequencing FASTQ file should be named: FILENAME.fastq.
        
        2. Paried end sequencing FASTQ files should be named: FILENAME_R1.fastq and FILENAME_R2.fastq.
        
        3. I suggest making copies of your FASTQ files rather then renaming the originals in case a mistake is made during the renaming process.
        
        4. For this nameing example "FILENAME" will be used as the input for the -n option below.

3. Open the Terminal and download and run the Docker container for MutantHunter by copying and pasteing the following command into the terminal: `docker run -it -v /PATH_TO_DESKTOP/Analysis_Directory:/MutantHunter/Analysis_Directory mellison/mutant_hunter:version1'

4. Run MutantHunter by running the code below.
```
./MutantHunter.sh \
    -n FILENAME \
    -g /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
    -f /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/genome.fa \
    -r single \
    -s 10 \
    -p /MutantHunter/MutantHunter/S_cerevisiae_Bowtie2_Index_and_FASTA/ploidy_n1.txt \
    -d /MutantHunter/Analysis_Directory \
    -o /MutantHunter/Analysis_Directory/NAME_YOUR_OUTPUT_FOLDER
```

#### Because the files and directory structure were set up ahead of time (on the Desktop and during the docker build), and all run out of the Docker container, the file paths in the above commands will all stay the same, but some of the options may change depending upon your needs. 



# MutantHunter Options

## All options are required

#### -n 

The -n option takes the prefix of the FASTQ file name for the wild-type strain. For the example of FILENAME.fastq or FILENAME_R1.fastq this prefix would simply be "FILENAME".

#### -g

The -g option takes the file PATH to the bowtie index files and the file prefix (genome). Use exactly what is shown above for this command.

#### -f

The -f option takes the file PATH and file name of the genome FASTA file (genome.fa) Use exactly what is shown above for this command.

#### -r

The -r option specifies whether the input data contains paired-end or single-end reads and can take values of "paired" or "single".

#### -s

The -s option takes a score cuttoff for the variant scores. 

#### -p


#### -d


#### -o









# Alternative Setup


## Setting up Directories for your Analysis and Organizing your Input Files

1. Place this folder in a directory of your choice. For these instructions I am going to refer this directory as `/AnalysisDirectory`

2. Download and unzip MutantHunter code from this website

    1. Run the following command in the terminal `wget https://github.com/mae92/MutantHunter/archive/master.zip` (by copying and pasting it into the terminal window and pressing enter) 
    
    2. Alternatively you can download from the website by performing the following steps
    
        1. On this GitHub page click on the green button labeled "Clone or Download"
        
            ![picture alt](https://github.com/mae92/MutantHunter/blob/master/images/image1.png)
            
        2. Click on the "Download ZIP" button in the dropdown menu
            
            ![picture alt](https://github.com/mae92/MutantHunter/blob/master/images/image2.png)
            
3. Place the `/MutantHunter-master` folder into `/AnalysisDirectory` (full path: `/Analysis_Directory/MutantHunter-master`)

4. Create folders in `/Analysis_Directory` to place your data into. FASTQ file directory `/FASTQ` and place all FASTQ files into it (See instructions above)


## Installing Docker

Follow the instructions at https://docs.docker.com/get-docker/ to download and install Docker on your computer

More infromation about how to work with Docker Containers can be found here: https://docs.docker.com


## Building and Running the Container Locally

1. Open the terminal

2. Set working directory to /Analysis_Directory/MutantHunter-master/DockerFile by typing `cd /Analysis_Directory/MutantHunter-master/DockerFile`

3. Run the following command (by copying and pasting it into the terminal window and pressing enter) to build the Docker container on your computer: `docker build -t mutant_hunter - < Dockerfile-MutantHunter` If you are promted for a password (by sudo) just type in your computer password.

4. Run this command to get info on the container you just built: `docker image ls`

5. Run this command to run the container `sudo docker run -it -v /Users/mitchellellison/Desktop/Analysis_Directory:/MutantHunter/Analysis_Directory mutant_hunter` 


## Running Mutant Hunter

1. From within the Docker conatiner navigate to the directory containing the MutantHunter software by running the following command `cd /MutantHunter/Analysis_Directory/MutantHunter-master/MutantHunter_Code`

2. Use the following command to make sure MutantHunter.sh is executable `chmod 777 MutantHunter.sh`

3. Run MutantHunter by running the code below.

```
MutantHunter.sh \
    -n PREFIX_FOR_YOUR_WT_FILE \
    -g /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome \
    -f /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/genome.fa \
    -r single \
    -s 10 \
    -p /MutantHunter/Analysis_Directory/MutantHunter-master/S_cerevisiae_Bowtie2_Index_and_FASTA/ploidy_n1.txt \
    -d /MutantHunter/Analysis_Directory/My_Analysis \
    -o /MutantHunter/Analysis_Directory/NAME_YOUR_OUTPUT_FOLDER
```

#### Because the files and directory structure were set up ahead of time and all run out of the Docker container the file paths in the above commands will all stay the same but some of the options may change depending upon your needs. 





