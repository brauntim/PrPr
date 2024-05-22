import requests
from rdkit import Chem
from rdkit.Chem import Descriptors

def inchi_key_to_inchi(inchi_key):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{inchi_key}/inchi"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        return None

def get_formula_and_mol_weight(inchi_key):
    # Konvertiere InChIKey zu InChI
    inchi = inchi_key_to_inchi(inchi_key)
    if not inchi:
        return None, None

    # Erzeuge ein RDKit-Mol-Objekt aus dem InChI
    mol = Chem.MolFromInchi(inchi)
    if not mol:
        return None, None

    # Extrahiere die Summenformel
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

    # Berechne das Molekulargewicht
    mol_weight = Descriptors.MolWt(mol)

    return formula, mol_weight

# Beispiel-InChIKey
inchi_key = "TXOFSCODFRHERQ-UHFFFAOYSA-N"
formula, mol_weight = get_formula_and_mol_weight(inchi_key)

print("Summenformel:", formula)
print("Molekulargewicht:", mol_weight)
