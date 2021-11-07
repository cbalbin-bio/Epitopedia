#!/usr/bin/python3.10

# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import os
import json
import sqlite3
import tempfile
from functools import partial
import subprocess
from multiprocessing import Pool

from rich.progress import track

from epitopedia.app import config
from epitopedia.app.blastparser import BLASTParser
from epitopedia.app.config import console
from epitopedia.app.hitparser import parseHit
from epitopedia.app.MMCIFSeqs import MMCIFSeqs
from epitopedia.app.reduce import reduce_results
from epitopedia.utils.utils import remove_previous_files

# from epitopedia.viz.serve import write_html, serve_html
from epitopedia.app.args import args


def main():

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

            parseHitP = partial(
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
                        p.imap(parseHitP, hits),
                        total=len(hits),
                        description="Searching and aligning SeqBMM structural representatives",
                    )
                )

            # data = [parseHit(hit) for hit in track(hits)]

            data = [datum for datum in data if datum]
            data_m.append(data)

        data = [data for data in data_m if data]
        with open(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}.json", "w") as output_handle:
            json.dump({"parameters": vars(args), "results": data}, output_handle)

        config.con.close()

        # console.log(f"Number of hits remaining after surface accessibility filter {PDB_INPUT}: {len(hits)}")
        # console.log(f"Number of unique SeqBMM's identified {PDB_INPUT}: {len(hits)}")

        reduce_results(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}.json")

        with open(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}_best.json") as input_handle:
            data = json.load(input_handle)

        # write_html(f"{config.OUTPUT_DIR}/Epitopedia_{pdb_input_str}_output.html", data)
        print(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}_best.json")
        if not args.headless:

            os.environ["EPITOPEDIA_DATA_DIR"] = f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}_best.json"
            print(f"[bold green]View results in browser at http://0.0.0.0:{args.port}[/bold green]")
            subprocess.run(["flask", "run", "--host=0.0.0.0", f"--port={args.port}"])

    else:
        os.environ["EPITOPEDIA_DATA_DIR"] = args.view
        print(f"[bold green]View results in browser at http://0.0.0.0:{args.port}[/bold green]")
        subprocess.run(["flask", "run", "--host=0.0.0.0", f"--port={args.port}"])


# if __name__ == "epitopedia.run_epitopedia":
#     main()
# elif __name__ == "__main__":
#     main()
