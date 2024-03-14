#!/bin/sh
echo "Running analysis..."
echo "Starting barrnap"

# Defining variables from widget
input_file="$1"
reference_path="$2"
methods="$3"
threshold="$4"

# Defining working directory
working_directory="$PWD"

# Do txt with names of references bacteria from database
find "$reference_path" -type f -name "*.fasta" -exec readlink -f {} \; > "$reference_path"/genomes_list.txt

file_name="genomes_list.txt"
ref="$reference_path/$file_name"

### Barrnap
barrnap -i "$input_file" -kingdom bac -evalue 1e-6 -lencutoff 0.8 -reject 0.25 -outseq "rRNA.fasta" -threads 6
echo "Barrnap finished"

# Is rRNA empty?
if [ -s rRNA.fasta ]; then
  echo "Known rRNAs have been identified!"
else
  echo "No known rRNA has been identified!"
fi

### Cd-hit
echo "$methods"
if [ "$methods" = "cd_hit" ] || [ "$methods" = "all" ]; then
  echo "Starting Cd-hit"

  # Create multifasta
  cat "$reference_path"/*.fasta > multifasta.fasta

  cd-hit-est-2d -i rRNA.fasta -i2 multifasta.fasta -o cd-hit_output -T 6 -M 0 -c 0.95

  # Analysis of cd-hit results
  echo "Analysis of cd-hit"
  python3 cd_hit_analysis.py "$ref"
fi

if [ "$methods" = "cd_hit" ]; then
  # Print results of analysis
  if [ -s "cluster_similar_bacteria.txt" ]; then
    echo "The bacteria were assigned to the same clusters as the bacteria: "
    cat "cluster_similar_bacteria.txt"
    echo "Cd-hit finished!"
    #elapsed_time=$(echo "$end_time - $start_time" | bc)
  else
    echo "The unknown bacterium was not included in any of the same clusters as the reference bacterium!"
    echo "Cd-hit finished!"
    #elapsed_time=$(echo "$end_time - $start_time" | bc)
  fi
fi

## BLAST

if [ "$methods" = "blast" ] || [ "$methods" = "all" ]; then
  echo "Starting BLAST"
  makeblastdb -in multifasta.fasta -dbtype nucl -out blast_database
  start_time=$(date +%s.%N)
  blastn -query "$input_file" -db blast_database -out "$working_directory"/blast.txt -outfmt "6 qseqid sseqid pident length evalue bitscore" -num_threads 6
  end_time=$(date +%s.%N)

  # Analysis of blast results
  python3 blast_analysis.py "$threshold" "$ref"
fi


if [ "$methods" = "blast" ]; then
  if [ -s "blast_similar_final.txt" ]; then
    echo "A significant match with these bacteria was found for an unknown bacterium: "
    cat "blast_similar_final.txt"
    echo "BLAST finished"
    #elapsed_time=$(echo "$end_time - $start_time" | bc)
    #echo "Čas trvání: $elapsed_time sekund"
  else
    echo "No bacterium is identical enough to an unknown bacterium!"
    echo "BLAST finished"
    #elapsed_time=$(echo "$end_time - $start_time" | bc)
    #echo "Čas trvání: $elapsed_time sekund"
  fi
fi


### Final analysis

cd "$working_directory"

if [ "$methods" = "all" ]; then
  # delete duplicates
  sort -u blast_similar_final.txt -o blast_similar_final.txt
  sort blast_similar_final.txt > blast_similar_final_sorted.txt
  sort cluster_similar_bacteria.txt > cluster_similar_bacteria_sorted.txt

  comm -12 blast_similar_final_sorted.txt | comm -12 - cluster_similar_bacteria_sorted.txt > final.txt
  cat blast_similar_final.txt cluster_similar_bacteria.txt | sort -u > all_possible_bacteria.txt

  rm blast_similar_final_sorted.txt
  rm cluster_similar_bacteria_sorted.txt

elif [ "$methods" = "cd_hit" ]; then
  sort cluster_similar_bacteria.txt > all_possible_bacteria.txt
  sort cluster_similar_bacteria.txt > final.txt

elif [ "$methods" = "blast" ]; then
  # delete duplicates
  sort -u blast_similar_final.txt -o blast_similar_final.txt
  sort blast_similar_final.txt > all_possible_bacteria.txt
  sort blast_similar_final.txt > final.txt

fi

echo "Analysis completed!"
