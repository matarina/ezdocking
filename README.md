# Molecular Docking Pipeline

Automated protein–ligand docking using AutoDock Vina with binding-pocket prediction (PRANK/P2Rank) and basic interaction analysis.

## Environment Setup
1) Create and activate the Python environment (includes AutoDock Vina):
   ```bash
   conda env create -f environment.yml
   conda activate easydocker
   # or: pixi shell -e easydocker
   ```
2) Install P2Rank/prank manually (download from the official releases), ensure Java is available, and add its `bin` directory to your `PATH`. Confirm it works:
   ```bash
   command -v prank && prank --version
   ```
3) Install geostd yourself and place it at `./geostd` (the code uses this path by default):
   ```bash
   git clone https://github.com/phenix-project/geostd.git geostd
   ```
4) Ensure AutoDock Vina is on your `PATH` (if not using the conda package shown above):
   ```bash
   vina --version
   ```

## Single Docking
```bash
python molecular_docking.py --pdb <PDB_ID> --smiles "<SMILES>" [--prefix <outdir>] [--chain A] [--exhaustiveness 128]
# or
python molecular_docking.py --local_pdb /path/to/protein.pdb --smiles "<SMILES>" [--prefix <outdir>] [--chain A]
```
- `--pdb` downloads from RCSB; `--local_pdb` uses a local PDB (mutually exclusive).
- `--prefix` is optional; by default outputs go under `runs/<protein>_<timestamp>/`.
- Outputs include `<pdb>_ligand_vina_best.pdbqt` plus interaction CSV/HTML files in the run folder.
Example:
```bash
python molecular_docking.py --pdb 1a0q --smiles "CC(C)NC1=NC=NC2=C1N=CN2" --chain A --exhaustiveness 128
```

## Batch Docking
Prepare a text file with `ID <tab/space> SMILES`, one per line (comments start with `#`).
```bash
python batch_docking.py --pdb <PDB_ID> --batch_file batch_smiles.txt --chain A
# or
python batch_docking.py --local_pdb /path/to/protein.pdb --batch_file batch_smiles.txt --chain A
```
Each ligand is docked into its own `<ID>/` folder under `runs/` (or your chosen `--prefix` root) with the same outputs as single docking.
Example:
```bash
cat > batch_smiles.txt <<'EOF'
lig1 CC(C)NC1=NC=NC2=C1N=CN2
lig2 CCOC1=CC=CC=C1
EOF
python batch_docking.py --pdb 1a0q --batch_file batch_smiles.txt --chain A
```
