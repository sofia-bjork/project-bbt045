#!/usr/bin/env bash

## This specifies the project and what cluster to run on. Don't change it for this course!
#SBATCH -A C3SE2024-2-16 -p vera
## This specifies the number of nodes to use. Don't change it.
#SBATCH -n 1
## This specifies the time to allocate for the tasks in this script (format hours:mins:seconds)
## Calculate a reasonable amount of time you think this will take, and multiply by 1.5-ish so you have some margin
#SBATCH -t 20:00:00
## This is the name of the job in the scheduler ('fastqc')
#SBATCH -J Bowtie2_Bjork
## Here you can add your e-mail address and you _may_ get an e-mail when the job fails or is done 
#SBATCH --mail-user=bjorkso@student.chalmers.se
#SBATCH --mail-type=ALL
## Set the names for the error and output files. It can be smart to set a path to these to your project directory, which you can do by adding that path right after the '=' sign
#SBATCH --error=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/alignment_Bowtie2/job.%J.err
#SBATCH --output=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/alignment_Bowtie2/job.%J.out


###
#
# Title: Bowtie2.sh
# Date: 2025.03.09
# Author: Sofia Björk, from Canvas (Johan Bengtsson-Palme heavily inspired by Vi Varga)
#
# Description: 
# This script will run Bowtie2 using a pre-made index.
# 
# Usage: 
# sbatch Bowtie2.sh
#
###


### Set parameters
# This sets up your project folder as the working directory. Change the XX to your group number/name
# CORRECT THIS FILE PATH AS YOU NEED
WORKDIR=/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork;

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


### Running Bowtie2

## First identify a list of files. This command will list all FORWARD read files in the project data directory, remove the last part ('_1.fastq.gz') and store the list in $FILE_LIST
FILE_LIST=`ls /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq/*_1.fastq.gz;`

for R1 in $FILE_LIST; do
  
  # Extract the base name
  BN=$(basename "$R1")
  
  # Remove the _1.fastq.gz suffix to get the sample prefix 
  SAMPLE=${BN%%_1.fastq.gz}
  
  # Construct the path to the mate 2 file 
  R2="$WORKDIR/data/trimmed_RNA-seq/${SAMPLE}_2.fastq.gz"
  
  echo "Aligning $SAMPLE with Bowtie2..."
  
  # Run Bowtie2 in local alignment mode for the paired-end reads
  apptainer exec $CONTAINER_LOC bowtie2 --local \
    -x $WORKDIR/data/Bowtie_index/hg38_cdna_index \
    -1 "$R1" -2 "$R2" \
    -S $WORKDIR/data/alignment_Bowtie2/${SAMPLE}.sam \
    --threads 4
  
  # Convert the SAM to a sorted BAM and create an index using samtools
  apptainer exec $CONTAINER_LOC samtools view -bS $WORKDIR/data/alignment_Bowtie2/${SAMPLE}.sam | apptainer exec $CONTAINER_LOC samtools sort -o $WORKDIR/data/alignment_Bowtie2/${SAMPLE}.bam
  
  apptainer exec $CONTAINER_LOC samtools index $WORKDIR/data/alignment_Bowtie2/${SAMPLE}.bam
  
  # Remove the SAM file to save disk space
  rm $WORKDIR/data/alignment_Bowtie2/${SAMPLE}.sam

  
  echo "Done aligning $SAMPLE"
  echo "-----------------------------------------"
done


### Copy relevant files back, this is good practice but will actually not do anything for this specific script
cp $WORKING_TMP/* $WORKDIR;
