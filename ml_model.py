from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report

pdb_file = "pdb_files/1CRN.pdb"

amino_acid_codes = {
    "ALA": 1, "ARG": 2, "ASN": 3, "ASP": 4, "CYS": 5,
    "GLN": 6, "GLU": 7, "GLY": 8, "HIS": 9, "ILE": 10,
    "LEU": 11, "LYS": 12, "MET": 13, "PHE": 14, "PRO": 15,
    "SER": 16, "THR": 17, "TRP": 18, "TYR": 19, "VAL": 20
}

helix_residues = set()
sheet_residues = set()

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


parser = PDBParser(QUIET=True)
structure = parser.get_structure("Protein", pdb_file)

X = []
y = []

for model in structure:
    for chain in model:
        for residue in chain:
            if is_aa(residue):
                res_name = residue.get_resname()
                res_id = residue.get_id()[1]
                chain_id = chain.id

                if res_name in amino_acid_codes:
                    X.append([amino_acid_codes[res_name]])

                    if (chain_id, res_id) in helix_residues:
                        y.append("Alpha Helix")
                    elif (chain_id, res_id) in sheet_residues:
                        y.append("Beta Sheet")
                    else:
                        y.append("Coil / Loop")


print("MACHINE LEARNING MODEL")
print("----------------------")
print("Dataset created from PDB file")
print("Total samples:", len(X))

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.3, random_state=42
)

model = RandomForestClassifier(random_state=42)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)

print("\nModel training completed.")
print("\nClassification Report:")
print(classification_report(y_test, y_pred, zero_division=0))

print("\nExample Prediction")
print("------------------")
print("Amino acid: CYS")
prediction = model.predict([[5]])
print("Predicted secondary structure:", prediction[0])