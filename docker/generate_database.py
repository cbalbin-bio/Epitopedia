#!/usr/bin/python3.10

# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import argparse
import subprocess
from sqlite3.dbapi2 import complete_statement

parser = argparse.ArgumentParser(description="Epitopedia")

parser.add_argument("iedb_path", type=str, help="Path to iedb_public.sql database")
parser.add_argument("mmcif_dir", type=str, help="Path to directory containing mmcif files")
parser.add_argument("data_dir", type=str, help="Output path for database generation")


parser.add_argument("--afdb-dir", type=str, help="Optional path to directory containing afdb mmcif files")
args = parser.parse_args()

mounts = ["-v", f"{args.iedb_path}:/app/iedb", "-v", f"{args.mmcif_dir}:/app/mmcif", "-v", f"{args.data_dir}:/app/data"]
commands = ["cbalbin/epitopedia:beta", "generate_database"]
if args.afdb_dir:
    mounts += ["-v", f"{args.afdb_dir}:/app/afdb"]
    commands += ["--use-afdb"]


subprocess.run(["docker", "run", "--rm", "-it"] + mounts + commands)
