from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog, Crippen
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import os
import gzip
import warnings

warnings.filterwarnings("ignore", message="to-Python converter for class boost::shared_ptr")

# ---------- CONFIG ----------
INPUT_FILE = "data set.csv"        #your file with any format supported by workflow
OUTPUT_CSV = "ranked_results.csv"
OUTPUT_SDF = "filtered_results.sdf"
SCAFFOLD_LIMIT = 3   #can be changed 
MAX_MOLS = 10000000   #can be changed
PRINT_EVERY = 5000 
# ----------------------------


def load_molecules(path, max_mols=None):
    ext = os.path.splitext(path)[-1].lower()
    molecules = []

    if ext == ".smi":
        with open(path, "r", encoding="utf-8", errors="ignore") as handle:
            for i, line in enumerate(handle):
                if max_mols is not None and i >= max_mols:
                    break
                if not line.strip():
                    continue

                smiles = line.split()[0]
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    molecules.append((smiles, mol))

    elif ext == ".csv":
        df = pd.read_csv(path)
        if "smiles" not in df.columns:
            raise ValueError("CSV file must contain a column named 'smiles'.")

        for i, smiles in enumerate(df["smiles"]):
            if max_mols is not None and i >= max_mols:
                break
            if not isinstance(smiles, str) or not smiles.strip():
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                molecules.append((smiles, mol))

    elif ext == ".sdf":
        supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=True)
        for i, mol in enumerate(supplier):
            if max_mols is not None and i >= max_mols:
                break
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                molecules.append((smiles, mol))

            if PRINT_EVERY and (i + 1) % PRINT_EVERY == 0:
                print(f"Read {i + 1} molecules")

    elif ext == ".gz" and path.endswith(".sdf.gz"):
        with gzip.open(path, "rb") as handle:
            supplier = Chem.ForwardSDMolSupplier(handle, sanitize=True, removeHs=True)
            for i, mol in enumerate(supplier):
                if max_mols is not None and i >= max_mols:
                    break
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    molecules.append((smiles, mol))

                if PRINT_EVERY and (i + 1) % PRINT_EVERY == 0:
                    print(f"Read {i + 1} molecules")

    else:
        raise ValueError("Supported input formats are: .smi, .csv, .sdf, .sdf.gz")

    return molecules


METALS = {
    3, 4, 11, 12, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
    31, 37, 38, 39, 40, 47, 48, 49, 50, 55, 56
}


def keep_organic_single(mol):
    smiles = Chem.MolToSmiles(mol)

    if "." in smiles:
        return None

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in METALS:
            return None

    return mol


def make_catalog():
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
    return FilterCatalog.FilterCatalog(params)


CATALOG = make_catalog()


def has_alerts(mol):
    return CATALOG.HasMatch(mol)


def lipinski(mol):
    return (
        Descriptors.MolWt(mol) <= 500
        and Descriptors.MolLogP(mol) <= 5
        and rdMolDescriptors.CalcNumHBD(mol) <= 5
        and rdMolDescriptors.CalcNumHBA(mol) <= 10
    )


def veber(mol):
    return (
        rdMolDescriptors.CalcNumRotatableBonds(mol) <= 10
        and rdMolDescriptors.CalcTPSA(mol) <= 140
    )


def ghose(mol):
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    atoms = mol.GetNumAtoms()

    try:
        mr = Crippen.MolMR(mol)
    except Exception:
        mr = -1.0

    return 160 <= mw <= 480 and -0.4 <= logp <= 5.6 and 20 <= atoms <= 70 and 40 <= mr <= 130


def egan(mol):
    return Descriptors.MolLogP(mol) <= 5.88 and rdMolDescriptors.CalcTPSA(mol) <= 131


def muegge(mol):
    return (
        200 <= Descriptors.MolWt(mol) <= 600
        and -2 <= Descriptors.MolLogP(mol) <= 5
        and rdMolDescriptors.CalcNumHBA(mol) <= 10
        and rdMolDescriptors.CalcNumHBD(mol) <= 5
        and rdMolDescriptors.CalcNumRotatableBonds(mol) <= 15
        and rdMolDescriptors.CalcTPSA(mol) <= 150
    )


def bioavailability(mol):
    return 0.55 if (lipinski(mol) and veber(mol)) else 0.11


def murcko(mol):
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold is None or scaffold.GetNumAtoms() == 0:
            return None
        return Chem.MolToSmiles(scaffold)
    except Exception:
        return None


def analyze_library(path, out_csv=OUTPUT_CSV, out_sdf=OUTPUT_SDF, scaffold_limit=SCAFFOLD_LIMIT, max_mols=MAX_MOLS):
    molecules = load_molecules(path, max_mols=max_mols)
    print(f"Loaded {len(molecules)} molecules from {path}")

    results = []
    passed_records = []
    scaffold_counts = {}

    for idx, (smiles, mol) in enumerate(molecules, start=1):
        mol = keep_organic_single(mol)
        if mol is None:
            continue

        if has_alerts(mol):
            continue

        if not lipinski(mol) or not veber(mol):
            continue

        scaffold = murcko(mol)
        if scaffold:
            current_count = scaffold_counts.get(scaffold, 0)
            if current_count >= scaffold_limit:
                continue
            scaffold_counts[scaffold] = current_count + 1

        molwt = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        num_atoms = mol.GetNumAtoms()
        num_hba = rdMolDescriptors.CalcNumHBA(mol)
        num_hbd = rdMolDescriptors.CalcNumHBD(mol)
        num_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        molar_refractivity = Crippen.MolMR(mol)

        ghose_pass = ghose(mol)
        egan_pass = egan(mol)
        muegge_pass = muegge(mol)
        bio = bioavailability(mol)
        final_score = round(bio + 0.1 * (int(ghose_pass) + int(egan_pass) + int(muegge_pass)), 4)

        record = {
            "SMILES": smiles,
            "MolWt": molwt,
            "MolLogP": logp,
            "NumAtoms": num_atoms,
            "NumHBA": num_hba,
            "NumHBD": num_hbd,
            "NumRotatableBonds": num_rot,
            "TPSA": tpsa,
            "MolarRefractivity": molar_refractivity,
            "Lipinski": True,
            "Veber": True,
            "Ghose": ghose_pass,
            "Egan": egan_pass,
            "Muegge": muegge_pass,
            "Bioavailability": bio,
            "Final_Score": final_score,
            "Scaffold": scaffold,
        }
        results.append(record)

        mol_copy = Chem.Mol(mol)
        mol_copy.SetProp("SMILES", smiles)
        mol_copy.SetProp("MolWt", str(round(molwt, 4)))
        mol_copy.SetProp("MolLogP", str(round(logp, 4)))
        mol_copy.SetProp("NumAtoms", str(num_atoms))
        mol_copy.SetProp("NumHBA", str(num_hba))
        mol_copy.SetProp("NumHBD", str(num_hbd))
        mol_copy.SetProp("NumRotatableBonds", str(num_rot))
        mol_copy.SetProp("TPSA", str(round(tpsa, 4)))
        mol_copy.SetProp("MolarRefractivity", str(round(molar_refractivity, 4)))
        mol_copy.SetProp("Lipinski", "True")
        mol_copy.SetProp("Veber", "True")
        mol_copy.SetProp("Ghose", str(ghose_pass))
        mol_copy.SetProp("Egan", str(egan_pass))
        mol_copy.SetProp("Muegge", str(muegge_pass))
        mol_copy.SetProp("Bioavailability", str(bio))
        mol_copy.SetProp("Final_Score", str(final_score))
        if scaffold is not None:
            mol_copy.SetProp("Scaffold", scaffold)

        passed_records.append((final_score, smiles, mol_copy))

        if PRINT_EVERY and idx % PRINT_EVERY == 0:
            print(f"Processed {idx} of {len(molecules)} molecules")

    if not results:
        print("No molecules passed the filtering steps.")
        pd.DataFrame().to_csv(out_csv, index=False)

        writer = Chem.SDWriter(out_sdf)
        writer.close()
        return pd.DataFrame()

    df = pd.DataFrame(results).sort_values("Final_Score", ascending=False).reset_index(drop=True)
    df.to_csv(out_csv, index=False)

    passed_records.sort(key=lambda item: item[0], reverse=True)

    writer = Chem.SDWriter(out_sdf)
    for _, _, mol in passed_records:
        writer.write(mol)
    writer.close()

    print(f"Saved {len(df)} molecules to {out_csv}")
    print(f"Saved {len(df)} molecules to {out_sdf}")

    return df


if __name__ == "__main__":
    df = analyze_library(INPUT_FILE)
    print(df.head())
