#!/usr/bin/env bash

## This specifies the project and what cluster to run on. Don't change it for this course!
#SBATCH -A C3SE2024-2-16 -p vera
## This specifies the number of nodes to use. Don't change it.
#SBATCH -n 1
#SBATCH -c 8
## This specifies the time to allocate for the tasks in this script (format hours:mins:seconds)
## Calculate a reasonable amount of time you think this will take, and multiply by 1.5-ish so you have some margin
#SBATCH -t 5:00:00
## This is the name of the job in the scheduler ('fastqc')
#SBATCH -J fastqc
## Here you can add your e-mail address and you _may_ get an e-mail when the job fails or is done 
#SBATCH --mail-user=caesarn@student.chalmers.se
#SBATCH --mail-type=ALL
## Set the names for the error and output files. It can be smart to set a path to these to your project directory, which you can do by adding that path right after the '=' sign
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out


###
#
# Title: STAR_alignment_full_data.sh
# Date: 2024.02.28
# Author: Caesar Nelson inspired by Johan Bengtsson-Palme and Vi Varga
#
# Description: 
# This script will run STAR alignment on the entire choosen set of RNA-seq data
# 
###


### Set parameters
# This sets up your project folder as the working directory. Change the XX to your group number/name
# CORRECT THIS FILE PATH AS YOU NEED
WORKDIR=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/;

# files used
# location of the container
CONTAINER_LOC=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/containers/Bjork_container.sif;

# temp files directory variable, this specifies a name of a directory on the compute node that will be used to store data temporarily
WORKING_TMP=$TMPDIR/JOB_TMP;


### Load modules, the module purge just makes sure nothing unwanted is loaded on the compute node
module purge
#module load MODULE_NAME/module.version ...;


# Copy relevant files to $TMPDIR
# create a temporary directory to store output files
mkdir $WORKING_TMP;
cd $WORKING_TMP;
# This is a good place to copy files to the working directory on the compute node, if needed for the script. It is not needed in this example, so this is commented out.
#cp FILE_TO_COPY $WORKING_TMP

### 

## First identify a list of files. This command will list all FORWARD read files in the project data directory, remove the last part ('_1.fastq.gz') and store the list in $FILE_LIST


apptainer exec $CONTAINER_LOC bash /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/Sbatch_scripts/script_STAR_alignment.sh

### Copy relevant files back, this is good practice but will actually not do anything for this specific script
cp $WORKING_TMP/* $WORKDIR;