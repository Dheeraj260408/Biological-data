from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

pdb_file = "pdb_files/1CRN.pdb"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("Protein", pdb_file)

helix_residues = set()
sheet_residues = set()

# Read HELIX and SHEET records from PDB file
with open(pdb_file, "r") as file:
    for line in file:
        if line.startswith("HELIX"):
            start_chain = line[19].strip()
            start_res = int(line[21:25].strip())
            end_res = int(line[33:37].strip())

            for i in range(start_res, end_res + 1):
                helix_residues.add((start_chain, i))

        elif line.startswith("SHEET"):
            start_chain = line[21].strip()
            start_res = int(line[22:26].strip())
            end_res = int(line[33:37].strip())

            for i in range(start_res, end_res + 1):
                sheet_residues.add((start_chain, i))


alpha_count = 0
beta_count = 0
coil_count = 0
total_count = 0

print("\nSECONDARY STRUCTURE PREDICTION USING PDB")
print("----------------------------------------")

for model in structure:
    for chain in model:
        print("\nChain:", chain.id)
        print("Residue No\tAmino Acid\tStructure")
        print("------------------------------------------")

        for residue in chain:
            if is_aa(residue):
                res_name = residue.get_resname()
                res_id = residue.get_id()[1]
                chain_id = chain.id

                if (chain_id, res_id) in helix_residues:
                    structure_type = "Alpha Helix"
                    alpha_count += 1
                elif (chain_id, res_id) in sheet_residues:
                    structure_type = "Beta Sheet"
                    beta_count += 1
                else:
                    structure_type = "Coil / Loop"
                    coil_count += 1

                total_count += 1
                print(f"{res_id}\t\t{res_name}\t\t{structure_type}")


def print_bar(label, count, total):
    percent = (count / total) * 100 if total > 0 else 0
    bar_length = int(percent / 2)
    bar = "█" * bar_length
    print(f"{label:<15} | {bar:<50} {count} residues ({percent:.2f}%)")


print("\nTERMINAL GRAPH")
print("----------------------------------------")
print_bar("Alpha Helix", alpha_count, total_count)
print_bar("Beta Sheet", beta_count, total_count)
print_bar("Coil / Loop", coil_count, total_count)

print("\nSUMMARY")
print("----------------------------------------")
print("Total Residues :", total_count)
print("Alpha Helix    :", alpha_count)
print("Beta Sheet     :", beta_count)
print("Coil / Loop    :", coil_count)