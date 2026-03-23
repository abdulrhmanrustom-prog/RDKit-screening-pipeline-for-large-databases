# RDKit-screening-pipeline-for-large-databases
An RDKit-based pipeline for large-scale molecular library screening, including structural cleaning, PAINS/Brenk filtering, drug-likeness evaluation, scaffold diversity control, and compound ranking.
# Workflow
<img width="2813" height="2892" alt="Multi-Panel Vertical Timeline (Layout 2x2) (3)" src="https://github.com/user-attachments/assets/e301be31-1980-468f-bfc1-d52424c78b8e" />

### Requirements

Install the required Python packages:

```bash
pip install rdkit pandas
```

---

### Run the pipeline

Execute the script:

```bash
python pipeline.py
```

---

## Input

The pipeline supports multiple molecular data formats commonly used in cheminformatics workflows, including:

* **CSV files** containing a column named `smiles`
* **SMILES files** (`.smi`)
* **Structure Data Files** (`.sdf`)
* **Compressed SDF files** (`.sdf.gz`)

The input dataset should consist of valid chemical structures represented as SMILES strings or molecular objects. Invalid or improperly formatted entries will be automatically ignored during processing.

---

## Output

* CSV file containing:
Molecular descriptors (e.g., molecular weight, LogP, TPSA, HBA/HBD)
drug-likeness evaluation flags (Lipinski, Veber, Ghose, Egan, Muegge)
final ranking score
murcko scaffold information
* SDF file containing:
The cleaned and filtered molecular structures
ready for downstream applications such as molecular docking or further cheminformatics analysis

---

## Notes

* The maximum number of processed molecules is configurable within the script.
* Scaffold diversity is controlled using a user-defined per-scaffold limit.
* The ranking score is derived from predicted bioavailability combined with additional drug-likeness rules.
* Only molecules passing structural cleaning, medicinal chemistry filters, and primary drug-likeness criteria are retained.

---

## Citation

If you use this pipeline in your research, please cite it as:


