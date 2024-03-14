import sys
threshold = int(sys.argv[1])
ref = sys.argv[2]

def blast_analysis(threshold, ref):
    bacteria_dict = {}
    scaffold_actual = None

    with open("blast.txt", "r") as file:
        for line in file:
            columns = line.strip().split("\t")
            scaffold = columns[0].split("_scaffold_")[-1]

            if scaffold != scaffold_actual:
                names = []
                scaffold_actual = scaffold

            bacteria_name = columns[1].split("_scaffold")[0]
            similarity = float(columns[2])

            if bacteria_name not in names and similarity >= threshold:
                names.append(bacteria_name)
                if bacteria_name in bacteria_dict:
                    bacteria_dict[bacteria_name] += 1
                else:
                    bacteria_dict[bacteria_name] = 1

    # Set threshold
    threshold = round(int(scaffold_actual) * 0.8)

    # Sort bacteria
    sorted_bacteria = sorted(bacteria_dict.items(), key=lambda x: x[1], reverse=True)
    bacterial_names = []

    # Save bacteria with similarity > 98 which are in more than 50% of scaffolds
    with open("blast_analysis_output.txt", "w") as output_file:
        for bacterium, count in sorted_bacteria:
            if count > threshold:
                bacterial_names.append(bacterium)
                count = round((int(count)/int(scaffold_actual))*100, 2)
                output_file.write(f"{bacterium} {count}\n")

    #Otevri soubor s názvy referencnich genomu
    with open(ref) as genome_file:
        genome_lines = genome_file.readlines()

    #Create list for full names of bacteria
    names_full = []

    for line in genome_lines:
        for name in bacterial_names:
            if name in line:
                bacterie = line.split(".")
                bacterie = bacterie[0]
                bacterie = bacterie.split("/")
                bacterie =bacterie[-1]
                names_full.append(bacterie)

    # Ulozeni do textoveho souboru
    file_name = "blast_similar_final.txt"

    with open(file_name, "w") as file:
        for item in names_full:
            file.write(f"{item}\n")

blast_analysis(threshold, ref)
