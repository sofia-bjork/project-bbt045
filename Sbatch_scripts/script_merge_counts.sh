count_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts"
sub_dir_count="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts/sub_dir_count"

mkdir -p "$sub_dir_count"

# Move BAM and log files to the sub-directory
for file in "$count_dir"/*_sorted.bam "$count_dir"/*.log; do
    if [ -e "$file" ]; then  # Check if the file exists
        mv "$file" "$sub_dir_count/"
    fi
done

echo "Files moved successfully to $sub_dir_count"

# Define input directory containing the count files
input_dir="/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_Bjork/data/STAR_gene_counts"

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

echo "Merged count file is done"