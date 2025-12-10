# PzDE-HMM: A Simple Tool for Identifying Plasticizer-Degrading Enzyme Genes


PzDE-HMM is a lightweight Python wrapper for HMMER that detects plasticizer-degrading enzyme (PzDE) genes in protein sequences using curated Hidden Markov Models (HMMs). It combines `hmmsearch` with filtering and annotation to enable rapid identification of key enzymes (e.g., hydrolases, esterases, lipases) involved in plasticizer degradation.

This tool is designed for metagenomic studies, especially in environments like landfill leachate, soil, and wastewater, where microbial plastic degradation potential is of interest.


---
Dependence
python: >=2.7
HMMER3: >=3.3.2


Usage
python PzDE-HMM.py -i input.faa -o output_prefix [options]

Required Arguments
Option	Description
-i, --input	Input protein FASTA file (e.g., ORFs.faa)
-o, --output	Output prefix (e.g., results)
Optional Arguments
Option	Default	Description
-db, --hmm-db	data/PzDE-HMM.hmm	Path to HMM database
--evalue	1e-5	E-value threshold
--min-score	None	Minimum bit score (no filter by default)
--min-modelcov	None	Minimum model coverage (e.g., 0.7)
--min-seqcov	None	Minimum sequence coverage (e.g., 0.7)
-n, --nproc	4	Number of CPUs
--ko-map	data/hmm_label-KO.txt	KO annotation file
--symbol-map	data/hmm_label-symbol.txt	Gene symbol file


Project Structure
 
PzDE-HMM/
├── PzDE-HMM.py                 # Main script
├── test.fasta                  # Test input (small subset)
├── data/
│   ├── PzDE-HMM.hmm           # HMM database
│   ├── hmm_label-KO.txt         # KO annotations
│   └── hmm_label-symbol.txt     # Gene symbol annotations
└── README.md


Run with built-in test data
python PzDE-HMM.py -i test.fasta -o test_result --evalue 1e-9 --min-modelcov 0.35 --min-score 30 -n 8
Expected output:
test_result.domtblout: Raw HMMER domain table
test_result.filtered.csv: Filtered and annotated results (CSV)


