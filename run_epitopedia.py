#!/usr/bin/python3.9

# docker run -p 8888:8888 -v /Users/christianbalbin/bioinformatics/data:/app/data -v /Volumes/ssd/mmCIF:/app/mmcif -v /Users/christianbalbin/bioinformatics/output:/app/output -it mimicrypipeline bash

# docker run -p 5000:5000 -v /Users/christianbalbin/bioinformatics/data:/app/data -v /Volumes/ssd/mmCIF:/app/mmcif -v /Users/christianbalbin/bioinformatics/output:/app/output -it mimicrypipeline bash
# move tmalign and dices to data folder

import tempfile


import config
from MMCIFMeta import PDBSeqs
from blastparser import BLASTParser
from subepi_finder import hit_to_pdb
import sqlite3
from reduce import reduce_results
from dataclasses import dataclass
import json
from dataclasses_json import dataclass_json
import argparse
from rich.progress import track
from rich.console import Console
from rich import print
from flask import Flask, render_template


console = Console()


parser = argparse.ArgumentParser(description="Epitopedia")
parser.add_argument("PDB_IDS", type=str, nargs="+", help="List of PDB_IDS formatted as PDBID_CHAIN. e.g. 6xr8_A")
parser.add_argument("--span", type=int, default=5, help="Minimum span length to consider a mimic")
parser.add_argument("--rasa", type=float, default=0.20, help="Relative accessible surface area cutoff")
parser.add_argument("--rasa_span", type=int, default=3, help="Minimum span length for surface accessibility filter")
parser.add_argument("--taxid_filter", type=str, help="Filter all taxaony at or below the level described by this taxid")

args = parser.parse_args()


@dataclass_json
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
    internal_source_seq_acc: int


def parseHit(hit, pdb_seq, query_pdb_base, query_pdb_chain, pdb_input_str):
    # looping through each hit in the inital blast search of the query structure sequence against the short IEDB epitope sequences

    con = sqlite3.connect(config.SQLITE_DATABASE_DIR)
    cur = con.cursor()
    motifs = []
    motif_pdb_ranges = []
    motifs_acc = []

    for index, (start, end) in enumerate(hit.match_ranges):
        # convert to 0 based indexing, check to make sure that motif is at least 5 residues long
        if end - (start - 1) >= args.span:
            # extract the motif seq from the pdb structure used as a query against iedb
            query_motif = pdb_seq.seqsolv[start - 1 : end]
            # make sure there are not any disordered residues in the motif seq
            if "?" not in query_motif:
                # add it to the list of motifs to search for representative structures using the source epitope sequence against pdb later on, also add the corresponding seq res nums
                motifs.append(query_motif)
                motif_pdb_ranges.append(pdb_seq.seqnums[start - 1 : end])
                motifs_acc.append(pdb_seq.binaryrasa[start - 1 : end])

    # get epitope data on this hit against the iedb db

    cur.execute(f"SELECT * FROM IEDB_FILT WHERE epitope_id={hit.subject_accession}")
    ret = cur.fetchone()
    epitope = Epitope(*ret)

    # get mmseqs results for this hits source sequence against pdb, looking for representative structures

    cur.execute(
        f'SELECT query, target, qcov, pident, evalue, seqres, seqsolv, seqnums FROM EPI_PDB LEFT JOIN mmCIF_seqs ON EPI_PDB.target = mmCIF_seqs.pdb_id WHERE EPI_PDB.query = "{epitope.source_antigen_accession}";'
    )
    mmseq_db_rows = cur.fetchall()

    pdb_hits_motifs = []
    # for each of the sub motifs in this  upper level hit
    for motif, motif_pdb_range, motif_acc in zip(motifs, motif_pdb_ranges, motifs_acc):
        # map it to the representative pdb strcutrues

        pdb_hits = hit_to_pdb(
            motif,
            motif_pdb_range,
            motif_acc,
            query_pdb_base,
            query_pdb_chain,
            mmseq_db_rows,
            epitope,
            hit,
            pdb_input_str,
        )
        if pdb_hits:
            # if data is returned ( succefully mapped to representative pdb structures), sort it by TMscore
            pdb_hits_motifs.append(sorted(pdb_hits, key=lambda x: x.TMalign_TMscore, reverse=True))

    if pdb_hits_motifs:

        return {
            "blasthit": hit.to_dict(),
            "hitepitopedata": epitope.to_dict(),
            "pdb_motif_hits": [[motif_hit.to_dict() for motif_hit in motif] for motif in pdb_hits_motifs],
        }


if __name__ == "__main__":
    con = sqlite3.connect(config.SQLITE_DATABASE_DIR)

    PDB_INPUTS = args.PDB_IDS
    pdb_input_str = "_".join(PDB_INPUTS)
    data_m = []

    for PDB_INPUT in PDB_INPUTS:
        query_pdb_base = PDB_INPUT.split("_")[0].lower()
        query_pdb_chain = PDB_INPUT.split("_")[1]

        pdb_path = f"{config.PDB_DATABASE_DIR}/{query_pdb_base[1:3]}/{query_pdb_base}.cif"

        with console.status("[bold green]Extracting query protein..."):
            query_pdb_seq = PDBSeqs(query_pdb_base, query_pdb_chain, compute_acc=True)
            fp = tempfile.NamedTemporaryFile(mode="w", delete=False)
            file_path = fp.name
            fp.write(query_pdb_seq.seqres)
            fp.close()
        console.log("Query protein extracted")

        with console.status("[bold green]BlASTing query protein against EPI-SEQ..."):
            if args.taxid_filter:
                bp = BLASTParser(
                    file_path,
                    PDB_INPUT,
                    config.BLAST_DATABASE_DIR,
                    acc_seq=query_pdb_seq.rasa,
                    taxids=args.taxid_filter,
                )
            else:
                bp = BLASTParser(file_path, PDB_INPUT, config.BLAST_DATABASE_DIR, acc_seq=query_pdb_seq.rasa)

        hits = bp.gethits()
        console.log("Query protein BLASTed against EPI-SEQ")

        console.log(f"Number of unfiltered hits for {PDB_INPUT}: {len(hits)}")
        hits.tocsv(f"{config.OUTPUT_DIR}/Blast_episeq_{pdb_input_str}.tsv")

        with console.status("[bold green]Filtering hits by span length..."):
            hits = hits.filterbymatchlen(args.span)
            hits.tocsv(f"{config.OUTPUT_DIR}/Blast_length_{pdb_input_str}.tsv")
            console.log("Hits filtered by span length")
            console.log(f"Number of hits remaining after span length filter {PDB_INPUT}: {len(hits)}")

        with console.status("[bold green]Filtering hits by surface accessibility..."):
            hits = hits.filterbyacc(args.rasa_span, args.rasa)
            hits.tocsv(f"{config.OUTPUT_DIR}/Blast_length_acc_{pdb_input_str}.tsv")
            console.log("Hits filtered by surface accessibility")
            console.log(f"Number of hits remaining after surface accessibility filter {PDB_INPUT}: {len(hits)}")

        data = [
            parseHit(
                hit,
                pdb_seq=query_pdb_seq,
                query_pdb_base=query_pdb_base,
                query_pdb_chain=query_pdb_chain,
                pdb_input_str=pdb_input_str,
            )
            for hit in track(hits, description="Searching and aligning SeqBMM structural representatives")
        ]
        data = [datum for datum in data if datum]
        data_m.append(data)

    data = [data for data in data_m if data]
    with open(f"{config.OUTPUT_DIR}/EPISEQ_struct_fragment_{pdb_input_str}.json", "w") as output_handle:
        json.dump(data, output_handle)

    config.con.close()

    # console.log(f"Number of hits remaining after surface accessibility filter {PDB_INPUT}: {len(hits)}")
    # console.log(f"Number of unique SeqBMM's identified {PDB_INPUT}: {len(hits)}")

    reduce_results(f"{config.OUTPUT_DIR}/EPISEQ_struct_fragment_{pdb_input_str}.json")

    with open(f"{config.OUTPUT_DIR}/EPISEQ_struct_fragment_{pdb_input_str}_best_per_source_seq.json") as input_handle:
        data = json.load(input_handle)
    app = Flask(__name__)

    @app.get("/")
    def main():
        return render_template("index.html", data=data.items())

    app.run(host="0.0.0.0")

    print("[bold green]View results in browser at http://0.0.0.0:5000[/bold green]")
