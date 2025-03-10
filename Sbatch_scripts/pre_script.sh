#!/usr/bin/env bash
cd /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/rawdata_RNA-seq
mkdir -p fastqc_rawdata
fastqc -o fastqc_rawdata *.fastq.gz 
#multiqc for fastqc reports from rawdata
cd /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/rawdata_RNA-seq/fastqc_rawdata
multiqc . 
#Trimming off adaptors
cd /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/rawdata_RNA-seq
mkdir -p ../trimmed_RNA-seq  # Create output directory if it doesn't exist

for file in *_1.fastq.gz; do
  # Extract the base filename (remove _1.fastq.gz)
  base="${file%%_1.fastq.gz}"
  
  # Define corresponding R2 file
  file_R2="${base}_2.fastq.gz"

  # Run Cutadapt for paired-end reads
  cutadapt \
    -a AGATCGGAAGAG -A AGATCGGAAGAG \
    --minimum-length 25 \
    -o /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq/"${base}_trimmed_1.fastq.gz" \
    -p /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq/"${base}_trimmed_2.fastq.gz" \
    "$file" "$file_R2"
done

#fastqc for trimmed data
cd /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq
mkdir -p fastqc_trimmed
fastqc -o fastqc_trimmed *.fastq.gz 
#multiqc for fastqc reports from trimmed data
cd /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq/fastqc_trimmed
multiqc . 