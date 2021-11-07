#!/usr/bin/python3.10

# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import argparse
import subprocess


import argparse

parser = argparse.ArgumentParser(description="Epitopedia")
parser.add_argument("output_dir", type=str, help="Directory to write output files")
parser.add_argument("mmcif_dir", type=str, help="Path to directory containing mmcif files")
parser.add_argument(
    "data_dir",
    type=str,
    help="Directory containing epitopedia.sqlite3 and EPI-SEQ.fasta",
)
parser.add_argument(
    "--afdb-dir",
    type=str,
    help="Optional path to directory containing afdb mmcif files",
)
parser.add_argument(
    "--PDB-IDS",
    type=str,
    nargs="+",
    help="List of PDB_IDS formatted as PDBID_CHAIN. e.g. 6xr8_A",
)
parser.add_argument("--span", type=int, default=5, help="Minimum span length to consider a mimic")
parser.add_argument("--rasa", type=float, default=0.20, help="Relative accessible surface area cutoff")
parser.add_argument(
    "--rasa-span",
    type=int,
    default=3,
    help="Minimum span length for surface accessibility filter",
)
parser.add_argument(
    "--taxid-filter",
    type=str,
    nargs="+",
    help="Filter all taxaony at or below the level described by this taxid",
)
parser.add_argument("--rmsd", type=float, help="Max RMSD to consider match a structural mimic")
parser.add_argument("--view", type=str, help="View results from a previous run")

parser.add_argument("--port", type=int, default=5000, help="Port used by webserver")
# parser.add_argument("--use-afdb", action="store_true", help="Include AFDB in database generation")
parser.add_argument(
    "--gplddt",
    type=float,
    default=0.70,
    help="Minimum global pLDDT score a structure predicted by alphafold must have to be considered",
)
parser.add_argument(
    "--lplddt",
    type=float,
    default=0.90,
    help="Minimum average local pLDDT score a region predicted by alphafold must have to be considered",
)
parser.add_argument("--headless", action="store_true", help="Do not start up webserver")


args = parser.parse_args()

args_dict = vars(args)
mounts = [
    "-v",
    f"{args.output_dir}:/app/output",
    "-v",
    f"{args.mmcif_dir}:/app/mmcif",
    "-v",
    f"{args.data_dir}:/app/data",
]

base = ["docker", "run", "--rm", "-it"]

commands = ["cbalbin/epitopedia:beta", "run_epitopedia"]

if args.afdb_dir:
    mounts += ["-v", f"{args.afdb_dir}:/app/afdb"]
    commands += ["--use-afdb"]


if args.PDB_IDS:
    commands += ["--PDB-IDS", "".join(args.PDB_IDS)]

if args.span:
    commands += ["--span", args.span]

if args.rasa:
    commands += ["--rasa", args.rasa]

if args.rasa_span:
    commands += ["--rasa-span", args.rasa_span]

if args.taxid_filter:
    commands += ["--taxid-filter", " ".join(list(args.taxid_filter))]

if args.rmsd:
    commands += ["--rmsd", args.rmsd]


if args.view:
    commands += ["--view", args.view]

if args.gplddt:
    commands += ["--gplddt", args.gplddt]

if args.lplddt:
    commands += ["--lplddt", args.lplddt]

if args.headless:
    commands += ["--headless"]
else:
    base += ["-p", f"{args.port}:{args.port}"]


subprocess.run(base + mounts + [str(command) for command in commands])
