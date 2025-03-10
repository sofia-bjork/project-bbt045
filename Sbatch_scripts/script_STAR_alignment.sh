#!/usr/bin/env bash

# Define paths
GENOME_DIR="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_index"
FASTQ_DIR="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/trimmed_RNA-seq/"
OUTPUT_DIR="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/alignment_STAR"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through FASTQ files
for file1 in ${FASTQ_DIR}/*_1.fastq.gz; do
    # Extract sample name
    SAMPLE=$(basename "$file1" _1.fastq.gz)
    file2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"

    # Run STAR alignment
    STAR --runThreadN 8 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$file1" "$file2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "${OUTPUT_DIR}/${SAMPLE}_" \
         --sjdbGTFfile /cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/human_genome/Homo_sapiens.GRCh38.113.gtf \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM GeneCounts
done

# Set variables
gtf_file="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/human_genome/Homo_sapiens.GRCh38.113.gtf"
output_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts"
aligned_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/alignment_STAR/"

# Create output directory if it does not exist
mkdir -p "$output_dir"


# Loop through BAM files
for bam_file in ${aligned_dir}/*sortedByCoord.out.bam; do
    # Extract sample name
    SAMPLE=$(basename "$bam_file" _Aligned.sortedByCoord.out.bam)
    
    # Define output files
    sorted_bam_file="$output_dir/${SAMPLE}_sorted.bam"
    output_file="$output_dir/${SAMPLE}_htseq_counts.txt"
    log_file="$output_dir/${SAMPLE}_htseq-count.log"

    # Sort BAM file by name for htseq-count
    samtools sort -n -o "$sorted_bam_file" "$bam_file"

    # Check if samtools sort was successful
    if [ $? -ne 0 ]; then
        echo "Error: samtools sort failed for $bam_file" >&2
        continue
    fi

    # Run htseq-count
    htseq-count -f bam -r name -s no -t exon -i gene_id "$sorted_bam_file" "$gtf_file" > "$output_file" 2> "$log_file"

    # Check if htseq-count ran successfully
    if [ $? -eq 0 ]; then
        echo "htseq-count completed successfully. Output saved in $output_file"
    else
        echo "htseq-count failed for $bam_file. Check log file at $log_file" >&2
        echo "htseq-count log:"
        cat "$log_file"
    fi

done

count_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts_TEST"
sub_dir_count="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts_TEST/sub_dir_count"

mkdir -p "$sub_dir_count"

# Move BAM and log files to the sub-directory
for file in "$count_dir"/*_sorted.bam "$count_dir"/*.log; do
    if [ -e "$file" ]; then  # Check if the file exists
        mv "$file" "$sub_dir_count/"
    fi
done

echo "Files moved successfully to $sub_dir_count"

# Define input directory containing the count files
input_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts_TEST"

# Define output file
output_file="$input_dir/merged_counts.txt"

# Get list of count files
count_files=("$input_dir"/*_trimmed_htseq_counts.txt)

# Check if there are count files
if [ ${#count_files[@]} -eq 0 ]; then
    echo "Error: No count files found in $input_dir" >&2
    exit 1
fi

# Use the first file as a reference to check gene consistency
first_file="${count_files[0]}"
awk '{print $1}' "$first_file" | sort > /tmp/gene_list_reference.txt

# Merge files based on first column
awk -v ref_file="/tmp/gene_list_reference.txt" '
    BEGIN { OFS="\t"; while ((getline < ref_file) > 0) ref_genes[$1] = 1; close(ref_file) }
    FNR==1 { files++ }  # Track number of files
    {
        gene_id = $1
        count = $2
        if (!(gene_id in ref_genes)) {
            print "Error: Gene " gene_id " is missing from one or more files!" > "/dev/stderr"
            exit 1
        }
        counts[gene_id, FILENAME] = count
        genes[gene_id]  # Store unique gene IDs
    }
    END {
        # Print header
        printf "Gene_ID"
        for (i in ARGV) {
            if (ARGV[i] ~ /_trimmed_htseq_counts.txt$/) {
                sample_name = ARGV[i]
                sub(".*/", "", sample_name)  # Remove path
                sub("_trimmed_htseq_counts.txt", "", sample_name)  # Extract sample ID
                printf "\t%s", sample_name
            }
        }
        print ""

        # Print merged counts
        for (gene in genes) {
            split(gene, arr, SUBSEP)
            printf "%s", arr[1]  # Print gene ID
            for (i in ARGV) {
                if (ARGV[i] ~ /_trimmed_htseq_counts.txt$/) {
                    file = ARGV[i]
                    printf "\t%s", (counts[arr[1], file] ? counts[arr[1], file] : 0)
                }
            }
            print ""
        }
    }
' "${count_files[@]}" > "$output_file"

echo "Merged file saved as: $output_file"

echo "STAR_alignment_full_alignment is completed succesfully"