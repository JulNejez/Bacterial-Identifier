import sys
ref_path = sys.argv[1]
input_file = 'cd-hit_output.clstr'

def cd_hit_analysis(ref_path, input_file):
  with open(input_file, 'r') as cdhit_file:
      cluster_lines = cdhit_file.readlines()

  # Create list for bacteria names
  bacteria_names = []

  for line in cluster_lines:
    if line.startswith(">Cluster"):
      pass

    # Find out the designation of the bacteria
    else:
      parts = line.split("::")
      parts = parts[-1]
      parts = parts.split("_scaf")
      parts = parts[0]

      # Save to  the list
      bacteria_names.append(parts)

      # Removing duplicates
      bacteria_names = list(set(bacteria_names))

  # Open the file with the names of the reference genomes
  with open(ref_path) as genome_file:
    genome_lines = genome_file.readlines()

  # Create a list for full names of bacterias
  names_full = []

  for line in genome_lines:
    for name in bacteria_names:
      if name in line:
        bacterie = line.split(".")
        bacterie = bacterie[0]
        bacterie = bacterie.split("/")
        bacterie =bacterie[-1]
        names_full.append(bacterie)

  # Save to txt
  file_name = "cluster_similar_bacteria.txt"

  with open(file_name, "w") as file:
    for item in names_full:
      file.write(f"{item}\n")

cd_hit_analysis(ref_path, input_file)
