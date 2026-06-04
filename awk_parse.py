pdb_file = "pdb_files/1CRN.pdb"

print("AWK STYLE PDB PARSING")
print("---------------------")
print("Reading only HELIX and SHEET records from PDB file\n")

with open(pdb_file, "r") as file:
    for line in file:
        if line.startswith("HELIX"):
            print("HELIX RECORD FOUND:")
            print(line.strip())

        elif line.startswith("SHEET"):
            print("SHEET RECORD FOUND:")
            print(line.strip())

print("\nAWK-style parsing completed.")