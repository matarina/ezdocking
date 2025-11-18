# Molecular Docking Pipeline

Automated protein–ligand docking using AutoDock Vina with binding-pocket prediction (PRANK/P2Rank) and basic interaction analysis.

## Environment
- Create the env: `conda env create -f environment.yml` (or `pixi shell -e easydocker` with `pixi.toml`).
- Install helper tools with the provided script (downloads P2Rank/prank and geostd into the repo):
  ```bash
  bash scripts/setup_tools.sh
  export PATH="$PWD/p2rank_2.5.1/bin:$PATH"   # add prank to PATH
  ```
- AutoDock Vina binaries must be available in your environment.

## Single Docking
```bash
python molecular_docking.py --pdb <PDB_ID> --smiles "<SMILES>" [--prefix <outdir>] [--chain A] [--exhaustiveness 128]
# or
python molecular_docking.py --local_pdb /path/to/protein.pdb --smiles "<SMILES>" [--prefix <outdir>] [--chain A]
```
- `--pdb` downloads from RCSB; `--local_pdb` uses a local PDB (mutually exclusive).
- `--prefix` is optional; by default outputs go under `runs/<protein>_<timestamp>/`.
- Outputs include `<pdb>_ligand_vina_best.pdbqt` plus interaction CSV/HTML files in the run folder.

## Batch Docking
Prepare a text file with `ID <tab/space> SMILES`, one per line (comments start with `#`).
```bash
python batch_docking.py --pdb <PDB_ID> --batch_file batch_smiles.txt --chain A
# or
python batch_docking.py --local_pdb /path/to/protein.pdb --batch_file batch_smiles.txt --chain A
```
Each ligand is docked into its own `<ID>/` folder under `runs/` (or your chosen `--prefix` root) with the same outputs as single docking.
