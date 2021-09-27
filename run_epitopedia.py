#!/usr/bin/python3.9

# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.


import argparse
import json
import sqlite3
import tempfile
from functools import partial
from multiprocessing import Pool

from flask import Flask, render_template
from jinja2 import Environment, FileSystemLoader
from rich import print
from rich.progress import track

from epitopedia.app import config
from epitopedia.app.blastparser import BLASTParser
from epitopedia.app.config import console
from epitopedia.app.hitparser import parseHit
from epitopedia.app.MMCIFSeqs import MMCIFSeqs
from epitopedia.app.reduce import reduce_results
from epitopedia.utils.utils import remove_previous_files

parser = argparse.ArgumentParser(description="Epitopedia")
parser.add_argument("PDB_IDS", type=str, nargs="+", help="List of PDB_IDS formatted as PDBID_CHAIN. e.g. 6xr8_A")
parser.add_argument("--span", type=int, default=5, help="Minimum span length to consider a mimic")
parser.add_argument("--rasa", type=float, default=0.20, help="Relative accessible surface area cutoff")
parser.add_argument("--rasa_span", type=int, default=3, help="Minimum span length for surface accessibility filter")
parser.add_argument("--taxid_filter", type=str, help="Filter all taxaony at or below the level described by this taxid")
parser.add_argument("--view", type=str, help="View results from a previous run")
parser.add_argument("--use-afdb", action="store_true", help="Include AFDB in database generation")


args = parser.parse_args()


if __name__ == "__main__":

    if not args.view:
        con = sqlite3.connect(config.SQLITE_DATABASE_DIR)

        PDB_INPUTS = args.PDB_IDS
        pdb_input_str = "_".join(PDB_INPUTS)
        remove_previous_files(config, pdb_input_str)
        data_m = []

        for PDB_INPUT in PDB_INPUTS:
            query_pdb_base = PDB_INPUT.split("_")[0].lower()
            query_pdb_chain = PDB_INPUT.split("_")[1]

            pdb_path = f"{config.PDB_DATABASE_DIR}/{query_pdb_base[1:3]}/{query_pdb_base}.cif"

            with console.status("[bold green]Extracting query protein..."):
                query_pdb_seq = MMCIFSeqs(query_pdb_base, query_pdb_chain, compute_acc=True)
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
                        pdb_seqnums=query_pdb_seq.seqnums,
                        pdb_seqsolv=query_pdb_seq.seqsolv,
                    )
                else:
                    bp = BLASTParser(
                        file_path,
                        PDB_INPUT,
                        config.BLAST_DATABASE_DIR,
                        acc_seq=query_pdb_seq.rasa,
                        pdb_seqnums=query_pdb_seq.seqnums,
                        pdb_seqsolv=query_pdb_seq.seqsolv,
                    )

            hits = bp.gethits()
            console.log("Query protein BLASTed against EPI-SEQ")

            console.log(f"Number of unfiltered hits for {PDB_INPUT}: {len(hits)}")
            hits.tocsv(f"{config.OUTPUT_DIR}/EPI_SEQ_hits_{pdb_input_str}.tsv")

            with console.status("[bold green]Filtering hits by span length..."):
                hits = hits.filterbymatchlen(args.span)
                hits.tocsv(f"{config.OUTPUT_DIR}/EPI_SEQ_span_filt_hits_{pdb_input_str}.tsv")
                console.log("Hits filtered by span length")
                console.log(f"Number of hits remaining after span length filter {PDB_INPUT}: {len(hits)}")

            with console.status("[bold green]Filtering hits by surface accessibility..."):
                hits = hits.filterbyacc(args.rasa_span, args.rasa)
                hits.tocsv(f"{config.OUTPUT_DIR}/EPI_SEQ_span_filt_acc_hits_{pdb_input_str}.tsv")
                console.log("Hits filtered by surface accessibility")
                console.log(f"Number of hits remaining after surface accessibility filter {PDB_INPUT}: {len(hits)}")

            parseHit = partial(
                parseHit,
                span=args.span,
                pdb_seq=query_pdb_seq,
                query_pdb_base=query_pdb_base,
                query_pdb_chain=query_pdb_chain,
                pdb_input_str=pdb_input_str,
            )
            with Pool(12) as p:
                data = list(
                    track(
                        p.imap(parseHit, hits, chunksize=8),
                        total=len(hits),
                        description="Searching and aligning SeqBMM structural representatives",
                    )
                )

            data = [datum for datum in data if datum]
            data_m.append(data)

        data = [data for data in data_m if data]
        with open(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}.json", "w") as output_handle:
            json.dump(data, output_handle)

        config.con.close()

        # console.log(f"Number of hits remaining after surface accessibility filter {PDB_INPUT}: {len(hits)}")
        # console.log(f"Number of unique SeqBMM's identified {PDB_INPUT}: {len(hits)}")

        reduce_results(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}.json")

        with open(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}_best.json") as input_handle:
            data = json.load(input_handle)
        app = Flask(__name__)

        with open(f"{config.OUTPUT_DIR}/Epitopedia_{pdb_input_str}_output.html", "w") as handle:

            env = Environment(loader=FileSystemLoader("templates"))
            template = env.get_template("index.html")
            output_from_parsed_template = template.render(data=data.items())
            handle.write(output_from_parsed_template)

        @app.get("/")
        def main():
            return render_template("index.html", data=data.items())

        print("[bold green]View results in browser at http://0.0.0.0:5000[/bold green]")

        app.run(host="0.0.0.0")
    else:
        with open(args.view) as input_handle:
            data = json.load(input_handle)

        app = Flask(__name__)

        @app.get("/")
        def main():
            return render_template("index.html", data=data.items())

        print("[bold green]View results in browser at http://0.0.0.0:5000[/bold green]")

        app.run(host="0.0.0.0")
