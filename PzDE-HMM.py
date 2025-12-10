#!/usr/bin/env python
# -*- coding: utf-8 -*-
# PzDE-HMM: Simple HMM-based detector for plasticizer-degrading enzyme genes
# Usage: python PzDE-HMM.py -i input.faa -o output_prefix
# Compatible with Python 2.7 and 3.5+

__author__ = ("Ju Feng",
              "Kang Xiaoxi")
__version__ = 'v1.0'
__date__ = '2025.12.10'

import os
import sys
import subprocess
import argparse

def main():
    parser = argparse.ArgumentParser(description="PzDE-HMM: Detect plasticizer-degrading enzyme genes using HMMER")
    parser.add_argument('-i', '--input', required=True, help='Input protein FASTA file (e.g., ORFs.faa)')
    parser.add_argument('-o', '--output', required=True, help='Output prefix (e.g., results)')
    parser.add_argument('-db', '--hmm-db', default='data/PzDE-HMM.hmm', help='Path to HMM database [default: data/PzDE-HMM.hmm]')
    parser.add_argument('--evalue', type=float, default=1e-5, help='E-value threshold [default: 1e-5]')
    parser.add_argument('--min-score', type=float, default=None, help='Minimum bit score (e.g., 50)')
    parser.add_argument('--min-modelcov', type=float, default=None, help='Minimum model coverage (e.g., 0.7)')
    parser.add_argument('--min-seqcov', type=float, default=None, help='Minimum sequence coverage (e.g., 0.7)')
    parser.add_argument('-n', '--nproc', type=int, default=4, help='Number of CPUs [default: 4]')
    parser.add_argument('--ko-map', default='data/hmm_label-KO.txt', help='KO mapping file')
    parser.add_argument('--symbol-map', default='data/hmm_label-symbol.txt', help='Gene symbol mapping file')

    args = parser.parse_args()

    input_faa = args.input
    output_prefix = args.output
    hmm_db = args.hmm_db
    evalue = args.evalue
    min_score = args.min_score
    min_modelcov = args.min_modelcov
    min_seqcov = args.min_seqcov
    cpus = args.nproc
    ko_map_file = args.ko_map
    symbol_map_file = args.symbol_map

    # Step 1: Check input files
    for f, name in [(input_faa, "Input FASTA"), (hmm_db, "HMM database")]:
        if not os.path.exists(f):
            print("Error: %s not found: %s" % (name, f))
            sys.exit(1)

    # Step 2: Run hmmsearch
    domtblout = "%s.domtblout" % output_prefix
    print("Running hmmsearch on %s..." % input_faa)
    cmd = [
        "hmmsearch",
        "--domtblout", domtblout,
        "-E", str(evalue),
        "--cpu", str(cpus),
        hmm_db,
        input_faa
    ]
    with open(os.devnull, 'w') as devnull:
        try:
            subprocess.check_call(cmd, stdout=devnull, stderr=devnull)
        except subprocess.CalledProcessError:
            print("Error: hmmsearch failed. Is HMMER installed and in your PATH?")
            sys.exit(1)

    # Step 3: Parse and filter results
    print("Parsing and filtering HMM hits...")
    try:
        with open(domtblout, "r") as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    except IOError:
        print("Error: failed to read %s" % domtblout)
        sys.exit(1)

    filtered = []
    for line in lines:
        parts = line.split(None, 22)
        if len(parts) < 22:
            continue

        try:
            target = parts[0]
            query = parts[3]
            tlen = int(parts[2])
            qlen = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
        except (ValueError, IndexError):
            continue

        modelcov = float(hmm_to - hmm_from + 1) / qlen
        seqcov = float(ali_to - ali_from + 1) / tlen

        if full_score < min_score:
            continue
        if min_modelcov is not None and modelcov < min_modelcov:
            continue
        if min_seqcov is not None and seqcov < min_seqcov:
            continue

        filtered.append({
            'target': target,
            'query': query,
            'full_evalue': full_evalue,
            'full_score': full_score,
            'modelcov': round(modelcov, 4),
            'seqcov': round(seqcov, 4)
        })

    # Step 4: Load annotations
    ko_map = {}
    symbol_map = {}

    if os.path.exists(ko_map_file):
        try:
            with open(ko_map_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        parts = line.split(None, 1)
                        if len(parts) == 2:
                            ko_map[parts[0]] = parts[1]
        except IOError:
            print("Warning: could not read KO map file: %s" % ko_map_file)

    if os.path.exists(symbol_map_file):
        try:
            with open(symbol_map_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        parts = line.split(None, 1)
                        if len(parts) == 2:
                            symbol_map[parts[0]] = parts[1]
        except IOError:
            print("Warning: could not read symbol map file: %s" % symbol_map_file)

    # Step 5: Write final result
    output_file = "%s.filtered.csv" % output_prefix
    try:
        with open(output_file, "w") as f:
            f.write("target,query,full_evalue,full_score,modelcov,seqcov,KO,symbol\n")
            for hit in filtered:
                query = hit['query']
                ko = ko_map.get(query, "NA")
                symbol = symbol_map.get(query, "NA")
                f.write("%s,%s,%.2e,%.2f,%.4f,%.4f,%s,%s\n" % (
                    hit['target'], hit['query'],
                    hit['full_evalue'], hit['full_score'],
                    hit['modelcov'], hit['seqcov'],
                    ko, symbol))
        print("Done! Results saved to %s" % output_file)
    except IOError:
        print("Error: could not write output file: %s" % output_file)
        sys.exit(1)

if __name__ == "__main__":
    main()