import argparse
import subprocess

parser = argparse.ArgumentParser(description="Epitopedia")

parser.add_argument("output_dir", type=str, help="Directory to write output files")
parser.add_argument("mmcif_dir", type=str, help="Path to directory containing mmcif files")
parser.add_argument("data_dir", type=str, help="Directory containing epitopedia.sqlite3 and EPI-SEQ.fasta")
parser.add_argument("PDB_IDS", type=str, nargs="+", help="List of PDB_IDS formatted as PDBID_CHAIN. e.g. 6xr8_A")
parser.add_argument("--taxid-filter", type=str, help="Filter all taxaony at or below the level described by this taxid")


parser.add_argument("--afdb-dir", type=str, help="Optional path to directory containing afdb mmcif files")
parser.add_argument("--headless", action="store_true", help="Do not start up webserver")
args = parser.parse_args()

mounts = [
    "-v",
    f"{args.output_dir}:/app/output",
    "-v",
    f"{args.mmcif_dir}:/app/mmcif",
    "-v",
    f"{args.data_dir}:/app/data",
]
commands = ["epitopedia:afdb", "run_epitopedia"] + args.PDB_IDS


if args.afdb_dir:
    mounts += ["-v", f"{args.afdb_dir}:/app/afdb"]
    commands += ["--use-afdb"]

if args.taxid_filter:
    commands += ["--taxid-filter", args.taxid_filter]

if args.headless:
    commands += ["--headless"]

print(mounts + commands)
subprocess.run(["docker", "run", "--rm", "-it"] + mounts + commands)
