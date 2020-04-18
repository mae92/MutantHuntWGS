# MutantHunter

A bioinformatics pipeline for identification and characterization of mutations in Saccharomyces cerevisiae from FASTQ files from a wild-type strain and one or more mutant strains.


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

4. Create a folder in `/Analysis_Directory` named descriptively for you specific analysis. For these instructions I am going to refer this directory as `/My_Analysis` the full path to which would be `/Analysis_Directory/My_Analysis`

5. Create folders in `/Analysis_Directory/My_Analysis` to place your data into.

    1. Wild-type FASTQ file `/WT_FASTQ` (full path: `/Analysis_Directory/My_Analysis/WT_FASTQ`)
    
    2. Mutant(s) FASTQ files `Mutant_FASTQ` (full path: `/Analysis_Directory/My_Analysis/Mutant_FASTQ`)
    
    3. Ensure that FASTQ files adhere to the following naming convention
    
        1. Single end sequencing FASTQ file should be named: FILENAME.fastq
        
        2. Paried end sequencing FASTQ files should be named: FILENAME_R1.fastq and FILENAME_R2.fastq
        
        3. I suggest making copies of your FASTQ files rather then renaming the originals in case a mistake is made during the renaming process


## Installing Docker, Building and Running the Container


1. Follow the instructions at https://docs.docker.com/get-docker/ to download and install Docker on your computer

2. Open the terminal

3. Set working directory to /Analysis_Directory/MutantHunter-master/DockerFile by typing `cd /Analysis_Directory/MutantHunter-master/DockerFile`

4. Run the following command (by copying and pasting it into the terminal window and pressing enter) to build the Docker container on your computer: `docker build -t mutant_hunter - < Dockerfile-MutantHunter` If you are promted for a password (by sudo) just type in your computer password.

5. Run this command to get info on the container you just built: `docker image ls`

6. Run this command to run the container `sudo docker run -it -v /Users/mitchellellison/Desktop/Analysis_Directory:/MutantHunter/Analysis_Directory mutant_hunter` 

More infromation about how to work with Docker Containers can be found here: https://docs.docker.com


## Running Mutant Hunter


