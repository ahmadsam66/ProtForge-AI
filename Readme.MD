# ProtForge-AI

**ProtForge-AI is a Python-based protein sequence analysis tool that calculates physicochemical properties, estimates exploratory developability metrics, retrieves predicted 3D structures, and generates comprehensive PDF reports.**
<img src="ProtForge-AI Protein Analysis Workflow(3).png" alt="ProtForge-AI Workflow" />

## Overview

ProtForge-AI is a command-line tool designed for preliminary analysis of protein sequences from FASTA files. It calculates physicochemical properties, including molecular weight, isoelectric point (pI), instability index, GRAVY, approximate net charge, and molar extinction coefficients. It also evaluates amino acid composition, identifies candidate regions using simple flexibility and polarity rules, and generates Kyte–Doolittle hydrophobicity plots.

For sequences of up to 400 amino acids, the tool attempts to retrieve a predicted 3D structure through the public ESMFold API. When successful, the predicted structure and available secondary-structure annotations are included in the PDF report. Additional exploratory metrics include aliphatic index, Boman index, heuristic solubility and aggregation scores, protease cleavage-site counts, selected sequence-liability motifs, and simplified MHC-I binding and melting-temperature estimates. Optional SMILES generation is supported through RDKit. **These heuristic estimates are intended for preliminary research and must not be interpreted as experimentally validated or clinically reliable predictions.**

## Features

- Protein sequence analysis from FASTA files.
- Physicochemical property calculations using Biopython.
- Exploratory developability metrics.
- Candidate region identification using sequence-based heuristics.
- Predicted 3D structure retrieval through ESMFold.
- C-alpha trace visualization.
- Kyte–Doolittle hydrophobicity profiling.
- Automated PDF report generation for each accepted FASTA record.
- Optional SMILES generation using RDKit.

## Requirements

- Python 3.9 or newer recommended.
- Internet connection for ESMFold structure retrieval.
- Required packages:
  - `requests`
  - `numpy`
  - `matplotlib`
  - `biopython`
  - `reportlab`
- Optional package: `rdkit`

## Installation on Linux

### 1. Clone the repository

Replace the example repository URL with your actual GitHub repository URL.

```bash
git clone https://github.com/ahmadsam66/ProtForge-AI.git
cd ProtForge-AI
```

Alternatively, download `ProtForge-AI.py` directly and navigate to its directory.

### 2. Create a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

### 3. Install dependencies

```bash
python -m pip install requests numpy matplotlib biopython reportlab
```

### 4. Install RDKit (optional)

```bash
python -m pip install rdkit
```

If pip installation is unsuccessful, Conda is an alternative:

```bash
conda create -n protforge -c conda-forge python=3.11 rdkit requests numpy matplotlib biopython reportlab
conda activate protforge
```

Verify the core dependencies:

```bash
python -c "import requests, numpy, matplotlib, Bio, reportlab; print('Dependencies OK')"
```

## Usage

General syntax:

```bash
python ProtForge-AI.py --input /path/to/sequences.fasta
```

The short option is also supported:

```bash
python ProtForge-AI.py -i /path/to/sequences.fasta
```

### Example 1: Analyze a local FASTA file

```bash
python ProtForge-AI.py -i protein.fasta
```

The input file must exist in the current directory. Reports are written to `reports/`.

### Example 2: Use an absolute input path

```bash
python ProtForge-AI.py \
  --input /home/$USER/data/candidate_proteins.fasta
```

Replace the example path with the actual location of your FASTA file.

### Example 3: Analyze multiple protein sequences

```bash
python ProtForge-AI.py \
  -i /home/$USER/data/protein_candidates.fasta
```

A FASTA file may contain multiple records. The script processes each record separately and attempts to generate an individual report for each accepted sequence.

### Example 4: Run from a dedicated project directory

```bash
mkdir -p ~/projects/protforge
cd ~/projects/protforge

cp /path/to/ProtForge-AI.py .

python3 -m venv .venv
source .venv/bin/activate

python -m pip install requests numpy matplotlib biopython reportlab

python ProtForge-AI.py \
  -i /home/$USER/data/protein_candidates.fasta

find reports -maxdepth 1 -type f -print
```

Install RDKit in the same environment if you also want SMILES files.

## Input

The main input is a FASTA-formatted text file containing one or more protein sequences.

Example:

```fasta
>Protein_A
MKWVTFISLLFLFSSAYS
>Protein_B
MALWMRLLPLLALLALWGPGPG
```

Input requirements:

- Each sequence record must begin with a header line starting with `>`.
- Sequences should use the standard one-letter codes for the 20 canonical amino acids.
- Sequences are converted to uppercase during processing.
- Records containing non-standard amino acid characters are skipped by the current validation logic.
- The FASTA record ID is used to name output files.
- The required command-line argument is `-i` or `--input`.

## Output

The script creates a `reports/` directory relative to the current working directory.

| Output | Description |
|---|---|
| `<ID>_report.pdf` | PDF report containing sequence properties, heuristic metrics, candidate regions, and available plots and structural information |
| `<ID>.smi` | Optional SMILES file generated when RDKit is installed and conversion succeeds |
| Terminal messages | Processing status, skipped sequences, structure retrieval status, and output information |

The PDF report may include:

1. Physicochemical properties.
2. Simplified ADMET and drug-likeness proxies.
3. Candidate regions for further investigation.
4. Available secondary-structure information.
5. Binding-related sequence descriptors.
6. Exploratory developability metrics.
7. Predicted 3D structure visualization, if retrieval succeeds.
8. Hydrophobicity profile, if the sequence is long enough.

**Note:** The current script does not provide a command-line option for changing the output directory. The `reports/` directory is created relative to the directory from which the script is executed.

## Potential Applications

- Preliminary screening of candidate protein and peptide sequences.
- Comparative analysis of physicochemical properties.
- Exploratory investigation of hydrophobicity and polar regions.
- Early-stage protein and peptide research.
- Automated generation of sequence-analysis reports.
- Research documentation and preliminary candidate prioritization.

## Scientific Limitations

ProtForge-AI is an exploratory analysis tool and does not replace specialized prediction software or experimental validation.

- **ADMET estimates:** LogP, solubility, blood–brain barrier penetration, and CYP450 risk are estimated using simplified amino acid-based heuristics. They are not validated pharmacological predictions for proteins or peptides.
- **MHC-I binding:** The simplified motif-based method is not a substitute for validated tools such as NetMHCpan.
- **Melting temperature:** The reported Tm is a rough heuristic estimate, not a validated thermodynamic prediction.
- **Solubility and aggregation:** These scores are approximate sequence-based indicators and may not correlate with experimental measurements.
- **Candidate binding regions:** Flexibility and hydrophilicity alone do not establish solvent accessibility, epitope identity, binding affinity, or biological activity.
- **Secondary structure:** The script parses available `HELIX` and `SHEET` records from PDB text. If these annotations are unavailable, it uses a simplified Biopython fallback.
- **3D structure:** Structure retrieval requires internet access and depends on the availability of the public ESMFold service. Predicted structures are not experimentally determined structures.
- **Structure visualization:** The plot represents the C-alpha trace rather than a complete molecular structure.
- **SMILES generation:** RDKit output may not fully represent post-translational modifications, disulfide connectivity, terminal modifications, or cyclic peptide chemistry.
- **Input validation:** Empty or very short sequences may cause errors or incomplete calculations.
- **Data privacy:** Sequences submitted for structure prediction are sent to a public external service. Do not submit confidential or unpublished sequences without appropriate authorization.

## Troubleshooting

### Missing Python modules

Install dependencies in the same environment used to run the script:

```bash
python -m pip install requests numpy matplotlib biopython reportlab
```

### SMILES files are missing

RDKit is optional. Check its availability:

```bash
python -c "from rdkit import Chem; print('RDKit available')"
```

### No 3D structure appears in the report

Check internet connectivity, ESMFold service availability, sequence length, and the terminal messages. The script only attempts structure retrieval for sequences up to 400 amino acids.

### Reports are not in the expected location

The output directory is relative to the current working directory:

```bash
pwd
find reports -maxdepth 1 -type f -print
```

## License and Citation

Specify an appropriate license before distributing this project publicly. This README does not assign a license to the source code.

For scientific publications, document the software version, relevant parameters, data sources, and methodological limitations.

Author :
Seyed Ahmad Mousavi
ahmad.moousavi@gmail.com

---

**Project:** ProtForge-AI  
**Language:** Python  
**Input format:** FASTA  
**Primary output:** PDF reports  
**Optional output:** SMILES files  
