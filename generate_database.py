#!/usr/bin/python3.9

# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.


### Steps to generate database used by pipeline written as one long monolithic procedural python script
import sqlite3
import subprocess
from dataclasses import dataclass

import os
import glob
from Bio.PDB.DSSP import make_dssp_dict
from multiprocessing import Pool
from rich.progress import track
from rich.console import Console

console = Console()
# convert mysql into sqlite3 db


with console.status("[bold green]Converting IEDB to sqlite3..."):
    subprocess.run("mysql2sqlite /app/iedb | sqlite3 /app/data/iedb_public.sqlite3", shell=True)

console.log("IEDB converted to sqlite3")


# Data model for each representative unique epitope
@dataclass
class Epitope:
    epitope_id: int
    description: str
    linear_peptide_seq: str
    linear_peptide_modified_seq: str
    linear_peptide_modification: str
    source_antigen_accession: str
    starting_position: int
    ending_position: int
    database: str
    name: str
    organism_id: int
    organism_name: str
    sequence: str
    internal_acc: int = None


## Create flatten and filter iedb_public.sqlite3 into indexed local DB
with console.status("[bold green]Flattening and filtering database IEDB..."):
    # Establish connection to db
    con = sqlite3.connect("/app/data/iedb_public.sqlite3")
    cur = con.cursor()
    # Create Joined tables with required fields
    cur.execute(
        f"""
CREATE TABLE epitope_data AS
SELECT DISTINCT epitope.epitope_id, epitope.description, epitope.linear_peptide_seq, epitope.linear_peptide_modified_seq, epitope.linear_peptide_modification, epitope_object.source_antigen_accession, object.starting_position, object.ending_position, source.database, source.name, source.organism_id, source.organism_name, source.sequence
FROM epitope
INNER JOIN epitope_object ON  epitope.epitope_id = epitope_object.epitope_id
INNER JOIN object  ON epitope_object.object_id = object.object_id
INNER JOIN source ON epitope_object.source_antigen_accession = source.accession
WHERE epitope.epitope_id IS NOT NULL
AND epitope.linear_peptide_seq IS NOT NULL
AND epitope_object.source_antigen_accession IS NOT NULL
AND object.starting_position IS NOT NULL
AND object.ending_position IS NOT NULL
AND source.database IS NOT NULL
AND source.name IS NOT NULL
AND source.organism_id IS NOT NULL
AND source.organism_name IS NOT NULL
ORDER BY epitope.epitope_id DESC;
"""
    )

    con.commit()
    cur.execute("CREATE INDEX idx_epitope_id ON epitope_data(epitope_id);")

    con.commit()

    cur.execute(
        f"""
CREATE TABLE epitope_assays AS
SELECT epitope_id, bcell.as_char_value, tcell.as_char_value, mhc_bind.as_char_value, mhc_elution.as_char_value FROM epitope_object
LEFT JOIN curated_epitope ON epitope_object.object_id = curated_epitope.e_object_id
LEFT JOIN bcell ON curated_epitope.curated_epitope_id = bcell.curated_epitope_id
LEFT JOIN tcell ON curated_epitope.curated_epitope_id = tcell.curated_epitope_id
LEFT JOIN mhc_bind ON curated_epitope.curated_epitope_id = mhc_bind.curated_epitope_id
LEFT JOIN mhc_elution ON curated_epitope.curated_epitope_id = mhc_elution.curated_epitope_id
"""
    )

    cur.execute("CREATE INDEX idx_epitope_id_assay ON epitope_assays(epitope_id);")

    con.commit()
    con.close()

console.log("IEDB flattened and filtered")
# closing and reopening connection to DB to ensure all changes have been comitted
con = sqlite3.connect("/app/data/iedb_public.sqlite3")
cur = con.cursor()


# function for determining if epitope_id has at least 1 positive assay
positive_strings = ["Positive", "Positive-Low", "Positive-High", "Positive-Intermediate"]


def hasPositiveAssay(epitope_id):
    for row in cur.execute(
        f'SELECT "as_char_value", "as_char_value:1", "as_char_value:2", "as_char_value:3" FROM epitope_assays WHERE epitope_id = {epitope_id}'
    ):
        if (
            (row[0] in positive_strings)
            or (row[1] in positive_strings)
            or (row[2] in positive_strings)
            or (row[3] in positive_strings)
        ):
            return True
    return False


# Extract and store epitope ids, filter out those with only negative assays
# TODO: add possitve assau filter back
epitope_ids = []
cur.execute("SELECT DISTINCT epitope_id FROM epitope_data")
rows = cur.fetchall()

for row in track(rows, description="Extracting epitopes with positive assays..."):

    if row[0] == None:
        continue
    if hasPositiveAssay(row[0]):
        epitope_ids.append(row[0])
epitope_ids = set(epitope_ids)
console.log("Epitopes with positive assays extracted")


# Each epitope ID may have multiple rows due to both reduntant references and taxa. Keeping only the first reference for each taxa
epitopes = []
for epitope_id in track(epitope_ids, description="Populating epitope data structure..."):
    epitope_version = 1
    taxid_set = []
    for row in cur.execute(f"SELECT * FROM epitope_data WHERE epitope_id = {epitope_id}"):
        organism_id = int(row[10])
        if organism_id in taxid_set:
            continue
        else:
            taxid_set.append(int(row[10]))
            epitopes.append(
                Epitope(
                    f"{row[0]}.{epitope_version}",
                    row[1],
                    row[2],
                    row[3],
                    row[4],
                    row[5],
                    row[6],
                    row[7],
                    row[8],
                    row[9],
                    row[10],
                    row[11],
                    row[12],
                )
            )
            epitope_version += 1

console.log("Epitope data structures populated")

# assigning internal accesion number to each unique epitope source sequence
source_internal_acession = {}
counter = 0
for epitope in track(epitopes, description="Assigning internal accesion numbers..."):
    if epitope.sequence in source_internal_acession:
        epitope.internal_acc = source_internal_acession[epitope.sequence]
    else:
        source_internal_acession[epitope.sequence] = counter
        epitope.internal_acc = counter

        counter += 1
console.log("Internal accession numbers assigned")
# Write out data to csv file
with open("/app/data/IEDB_FILT.tsv", "w") as handle:
    handle.write(
        "epitope_id\tdescription\tlinear_peptide_seq\tlinear_peptide_modified_seq\tlinear_peptide_modification\tsource_antigen_accession\tstarting_position\tending_position\tdatabase\tname\torganism_id\torganism_name\tsequence\tinternal_sequence_acc\n"
    )
    for epitope in epitopes:
        handle.write(
            f"{epitope.epitope_id}\t{epitope.description}\t{epitope.linear_peptide_seq}\t{epitope.linear_peptide_modified_seq}\t{epitope.linear_peptide_modification}\t{epitope.source_antigen_accession}\t{epitope.starting_position}\t{epitope.ending_position}\t{epitope.database}\t{epitope.name}\t{epitope.organism_id}\t{epitope.organism_name}\t{epitope.sequence}\t{epitope.internal_acc}\n"
        )


# Write out fasta file
with open("/app/data/EPI_SEQ.fasta", "w") as handle:
    for epitope in epitopes:
        handle.write(
            f">{epitope.epitope_id} {epitope.database}:{epitope.source_antigen_accession}[{epitope.starting_position}:{epitope.ending_position}] {epitope.name} {epitope.organism_name}\n"
        )
        handle.write(f"{epitope.linear_peptide_seq}\n")

# Write out taxid map
with open("/app/data/EPI_SEQ_taxid_map.txt", "w") as handle:
    for epitope in epitopes:
        handle.write(f"{epitope.epitope_id} {epitope.organism_id}\n")

con.close()

with console.status("[bold green]Writing IEDB_FILT to DB..."):

    sql = f""".separator "\t"
.import /app/data/IEDB_FILT.tsv IEDB_FILT
CREATE INDEX idx_epitope_id ON IEDB_FILT(epitope_id);"""

    subprocess.run(["sqlite3", "/app/data/epitopedia.sqlite3"], input=sql.encode())

console.log("IEDB_FILT written to DB")
## generate blast db from epitope linear peptide sequences with a tax map
with console.status("[bold green]Generating EPI_SEQ..."):

    subprocess.run(
        [
            "makeblastdb",
            "-in",
            "/app/data/EPI_SEQ.fasta",
            "-parse_seqids",
            "-taxid_map",
            "/app/data/EPI_SEQ_taxid_map.txt",
            "-title",
            "EPI_SEQ",
            "-dbtype",
            "prot",
        ]
    )

console.log("EPI_SEQ generated")

# extract sequence data from mmcif files
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


paths = glob.glob("/app/mmcif/*/*.cif")

with open("/app/data/mmCIF_seqs.csv", "w") as f_cif_seq_dat:
    f_cif_seq_dat.write("pdb_id,pdb_id_asym_id,seqres,seqsolv,seqnums,asym_id,icode\n")
    for path in track(paths, description="Extracting sequences from mmCIF files for EPI_PDB generation"):
        name = os.path.basename(path).split(".")[0]
        in_data = False
        chain = []
        asym_id = []
        seqres = []
        seqsolv = []
        icode = []
        auth_seq_nums = []
        with open(path) as cif_file:
            for line in cif_file:
                if "_pdbx_poly_seq_scheme.asym_id" in line:
                    in_data = True

                if in_data:
                    if line.startswith("#"):
                        break
                    elif line.startswith("_"):
                        continue

                    line = line.split()
                    chain.append(line[9])
                    asym_id.append(line[0])
                    seqres.append(line[3])
                    icode.append(line[10])
                    auth_seq_nums.append(line[6])
                    seqsolv.append(line[7])

        if len(chain) == 0:
            continue

        cur_chain_id = chain[0]
        cur_asym_id = []
        cur_seqres = []
        cur_seqsolv = []
        cur_icode = []
        cur_auth_seq_nums = []

        for i, chain_ident in enumerate(chain):
            if chain_ident == cur_chain_id:
                cur_asym_id.append(asym_id[i])
                cur_seqres.append(seqres[i])
                cur_seqsolv.append(seqsolv[i])
                cur_icode.append(icode[i])
                cur_auth_seq_nums.append(auth_seq_nums[i])
            else:
                try:
                    cur_asym_id = " ".join(cur_asym_id)
                    cur_seqres = "".join([protein_3to1[res] for res in cur_seqres])
                    cur_seqsolv = "".join([protein_3to1[res] for res in cur_seqsolv])
                    cur_auth_seq_nums = " ".join(cur_auth_seq_nums)
                    cur_icode = " ".join(cur_icode)
                except KeyError:
                    cur_chain_id = chain_ident
                    cur_asym_id = []
                    cur_seqres = []
                    cur_seqsolv = []
                    cur_icode = []
                    cur_auth_seq_nums = []
                    cur_asym_id.append(asym_id[i])
                    cur_seqres.append(seqres[i])
                    cur_seqsolv.append(seqsolv[i])
                    cur_icode.append(icode[i])
                    cur_auth_seq_nums.append(auth_seq_nums[i])
                    continue

                f_cif_seq_dat.write(
                    f"{name}_{cur_chain_id},{name}_{cur_asym_id[0]},{cur_seqres},{cur_seqsolv},{cur_auth_seq_nums},{cur_asym_id},{cur_icode}\n"
                )

                cur_chain_id = chain_ident
                cur_asym_id = []
                cur_seqres = []
                cur_seqsolv = []
                cur_icode = []
                cur_auth_seq_nums = []
                cur_asym_id.append(asym_id[i])
                cur_seqres.append(seqres[i])
                cur_seqsolv.append(seqsolv[i])
                cur_icode.append(icode[i])
                cur_auth_seq_nums.append(auth_seq_nums[i])

        try:
            cur_asym_id = " ".join(cur_asym_id)
            cur_seqres = "".join([protein_3to1[res] for res in cur_seqres])
            cur_seqsolv = "".join([protein_3to1[res] for res in cur_seqsolv])
            cur_auth_seq_nums = " ".join(cur_auth_seq_nums)
            cur_icode = " ".join(cur_icode)
        except KeyError:
            continue

        f_cif_seq_dat.write(
            f"{name}_{cur_chain_id},{name}_{cur_asym_id[0]},{cur_seqres},{cur_seqsolv},{cur_auth_seq_nums},{cur_asym_id},{cur_icode}\n"
        )


console.log("Sequences from mmCIF files extracted for EPI_PDB generation")
with console.status("[bold green]Writing mmCIF_seqs to DB..."):
    sql = f""".separator ","
.import /app/data/mmCIF_seqs.csv mmCIF_seqs
CREATE INDEX idx_pdb_id ON mmCIF_seqs(pdb_id);
CREATE INDEX idx_pdb_id_asym_id ON mmCIF_seqs(pdb_id_asym_id);"""

    subprocess.run(["sqlite3", "/app/data/epitopedia.sqlite3"], input=sql.encode())

console.log("mmCIF_seqs written to DB")

con = sqlite3.connect("/app/data/epitopedia.sqlite3")
cur = con.cursor()


console.log("Generating EPI_PDB...")
# generate fasta file of PDB sequences for MMseqs
with open("/app/data/mmCIF_seqs.fa", "w") as handle:
    for row in cur.execute("SELECT DISTINCT pdb_id, seqres FROM mmCIF_seqs"):
        handle.write(f">{row[0]}\n{row[1]}\n")

# generate fasta of epitope source sequence
with open("/app/data/IEDB_source_seqs.fa", "w") as handle:
    for row in cur.execute("SELECT DISTINCT source_antigen_accession, sequence FROM IEDB_FILT"):
        handle.write(f">{row[0]}\n{row[1]}\n")
con.close()
# run mmseqs iedb against pdb
subprocess.run(["mmseqs", "createdb", "/app/data/IEDB_source_seqs.fa", "/app/data/IEDB_source_seqs_DB"])
subprocess.run(["mmseqs", "createdb", "/app/data/mmCIF_seqs.fa", "/app/data/mmCIF_seqs_DB"])
os.makedirs("/app/data/tmp", exist_ok=True)

subprocess.run(["mmseqs", "createindex", "/app/data/IEDB_source_seqs_DB", "/app/data/tmp"])
subprocess.run(["mmseqs", "createindex", "/app/data/mmCIF_seqs_DB", "/app/data/tmp"])

subprocess.run(
    [
        "mmseqs",
        "search",
        "/app/data/IEDB_source_seqs_DB",
        "/app/data/mmCIF_seqs_DB",
        "/app/data/EPI_PDB_DB",
        "/app/data/tmp",
        "-a",
        "-s",
        "7.5",
    ]
)

subprocess.run(
    [
        "mmseqs",
        "convertalis",
        "/app/data/IEDB_source_seqs_DB",
        "/app/data/mmCIF_seqs_DB",
        "/app/data/EPI_PDB_DB",
        "/app/data/EPI_PDB_DB.m8",
        "--format-output",
        "query,target,qcov,pident,evalue,cigar,qaln,taln",
    ]
)
console.log("EPI_PDB generated")

with open("/app/data/EPI_PDB_DB.m8") as input_handle, open(
    "/app/data/EPI_PDB_DB_filt.m8", "w"
) as output_handle, console.status("[bold green]Filtering EPI-PDB..."):
    output_handle.write("query\ttarget\tqcov\tpident\tevalue\tcigar\tqaln\ttaln\n")
    for line in input_handle:
        split = line.split("\t")
        qcov = split[2]
        pident = split[3]

        if float(pident) >= 90.0 and float(qcov) >= 0.20:
            output_handle.write(line)

console.log("EPI_PDB filtered")
with console.status("[bold green]Writing EPI-PDB to DB..."):
    sql = f""".separator "\t"
.import /app/data/EPI_PDB_DB_filt.m8 EPI_PDB
CREATE INDEX idx_query ON EPI_PDB(query);
CREATE INDEX idx_target ON EPI_PDB(target);"""

    subprocess.run(["sqlite3", "/app/data/epitopedia.sqlite3"], input=sql.encode())
console.log("EPI_PDB written to DB")
os.makedirs("/app/data/dssp", exist_ok=True)


def runDSSP(path):
    basename = os.path.basename(path).split(".")[0]

    subprocess.run(f"mkdssp {path} > /app/data/dssp/{basename}.dssp", shell=True, capture_output=True)


def parseDSSP(path):

    basename = os.path.basename(path).split(".")[0]

    chains_acc = []

    try:

        dssp_dict = make_dssp_dict(path)
        chains = {key[0] for key in dssp_dict[0].keys()}
    except:
        return chains_acc

    for chain in chains:

        try:

            cur.execute(
                f'SELECT seqnums, icode, asym_id, pdb_id FROM mmCIF_seqs WHERE pdb_id_asym_id = "{basename}_{chain}"'
            )
            row = cur.fetchone()
            if row == None:
                continue
            seq_nums = row[0].split(" ")
            icodes = row[1].split(" ")
            asym_ids = row[2].split(" ")
            pdb_id = row[3]

            chain_acc = []
            for seq_num, icode, asym_id in zip(seq_nums, icodes, asym_ids):
                if seq_num != "?":
                    try:
                        if icode != ".":

                            chain_acc.append(str(dssp_dict[0][(asym_id, (" ", int(seq_num), icode))][2]))
                        else:
                            chain_acc.append(str(dssp_dict[0][(asym_id, (" ", int(seq_num), " "))][2]))
                    except KeyError:
                        chain_acc.append("?")  # this should fix issue with weird residues
                else:
                    chain_acc.append("?")

            chain_acc = " ".join(chain_acc)

        except:
            continue

        chains_acc.append(f"{pdb_id},{chain_acc}\n")

    return chains_acc


if __name__ == "__main__":
    paths = glob.glob("/app/mmcif/*/*.cif")
    with Pool(12) as p:
        data = list(
            track(
                p.imap(runDSSP, paths, chunksize=1024),
                total=len(paths),
                description="Running DSSP on PDB (mmCIF) files...",
            )
        )

    console.log("DSSP ran on PDB (mmCIF) files")

    con = sqlite3.connect("/app/data/epitopedia.sqlite3")
    cur = con.cursor()

    dssp_files = glob.glob("/app/data/dssp/*")
    with Pool(12) as p:
        data = list(
            track(
                p.imap(parseDSSP, dssp_files, chunksize=1024),
                total=len(dssp_files),
                description="Parsing DSSP files...",
            )
        )
        data = [chains for chains in data if chains]
        with open("/app/data/PDB_DSSP.csv", "w") as output_handle:
            output_handle.write(f"pdb_id,acc\n")
            for chains in data:
                for chain in chains:
                    output_handle.write(chain)
    con.close()
    console.log("Parsed DSSP files")
    with console.status("[bold green]Writing PDB_DSSP to DB..."):
        sql = f""".separator ","
.import /app/data/PDB_DSSP.csv PDB_DSSP
CREATE INDEX idx_acc_pdb_id ON PDB_DSSP(pdb_id);"""

    subprocess.run(["sqlite3", "/app/data/epitopedia.sqlite3"], input=sql.encode())
    console.log("PDB_DSSP written to DB")
