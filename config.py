import sqlite3
import os
from rich.console import Console

console = Console()

PDB_DATABASE_DIR = "/app/mmcif"
BLAST_DATABASE_DIR = "/app/data/EPI_SEQ.fasta"
SQLITE_DATABASE_DIR = "/app/data/epitopedia.sqlite3"
OUTPUT_DIR = "/app/output"
DICE_DIR = "/app/output/dices"
TMALIGN_DIR = "/app/output/TMalign_results"

# PDB_DATABASE_DIR = "/Volumes/ssd/mmCIF"
# BLAST_DATABASE_DIR = "/Users/christianbalbin/bioinformatics/data/EPI_SEQ.fasta"
# SQLITE_DATABASE_DIR = "/Users/christianbalbin/bioinformatics/data/epitopedia.sqlite3"
# OUTPUT_DIR = "/Users/christianbalbin/bioinformatics/output"
# DICE_DIR = "/Users/christianbalbin/bioinformatics/epitope-molecular-mimicry-pipeline/dices"
# TMALIGN_DIR = "/Users/christianbalbin/bioinformatics/epitope-molecular-mimicry-pipeline/TMalign_results"


os.makedirs(DICE_DIR, exist_ok=True)
os.makedirs(TMALIGN_DIR, exist_ok=True)


con = sqlite3.connect(SQLITE_DATABASE_DIR)


protein_3to1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "?": "?",
}


residue_max_acc = {
    # Miller max acc: Miller et al. 1987 https://doi.org/10.1016/0022-2836(87)90038-6
    # Wilke: Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635
    # Sander: Sander & Rost 1994 https://doi.org/10.1002/prot.340200303
    "Miller": {
        "A": 113.0,
        "R": 241.0,
        "N": 158.0,
        "D": 151.0,
        "C": 140.0,
        "Q": 189.0,
        "E": 183.0,
        "G": 85.0,
        "H": 194.0,
        "I": 182.0,
        "L": 180.0,
        "K": 211.0,
        "M": 204.0,
        "F": 218.0,
        "P": 143.0,
        "S": 122.0,
        "T": 146.0,
        "W": 259.0,
        "Y": 229.0,
        "V": 160.0,
    },
    "Wilke": {
        "A": 129.0,
        "R": 274.0,
        "N": 195.0,
        "D": 193.0,
        "C": 167.0,
        "Q": 225.0,
        "E": 223.0,
        "G": 104.0,
        "H": 224.0,
        "I": 197.0,
        "L": 201.0,
        "K": 236.0,
        "M": 224.0,
        "F": 240.0,
        "P": 159.0,
        "S": 155.0,
        "T": 172.0,
        "W": 285.0,
        "Y": 263.0,
        "V": 174.0,
    },
    "Sander": {
        "A": 106.0,
        "R": 248.0,
        "N": 157.0,
        "D": 163.0,
        "C": 135.0,
        "Q": 198.0,
        "E": 194.0,
        "G": 84.0,
        "H": 184.0,
        "I": 169.0,
        "L": 164.0,
        "K": 205.0,
        "M": 188.0,
        "F": 197.0,
        "P": 136.0,
        "S": 130.0,
        "T": 142.0,
        "W": 227.0,
        "Y": 222.0,
        "V": 142.0,
    },
}
