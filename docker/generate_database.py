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
commands = ["epitopedia:afdb", "generate_database.py"]
if args.afdb_dir:
    mounts += ["-v", f"{args.afdb_dir}:/app/afdb"]
    commands += ["--use-afdb"]


subprocess.run(["docker", "run", "--rm", "-it"] + mounts + commands)

# mounts = [types.Mount("/app/iedb", args.iedb_path, type='bind'),types.Mount("/app/mmcif", args.mmcif_dir, type='bind'),types.Mount("/app/data", args.data_dir, type='bind')]
# commands = ["epitopedia/data/generate_database.py"]
# if args.afdb_dir:
#     mounts += [types.Mount("/app/afdb", args.afdb_dir, type='bind')]
#     commands += ["--use-afdb"]


# client = docker.from_env()

# container = client.containers.run(
#     image="epitopedia:afdb",
#     command=commands,
#     remove=True,
#     detach=False,
#     mounts=mounts,
#     stdin_open=True,

# )


# signal.signal(signal.SIGINT,lambda unused_sig, unused_frame: container.kill())

# for line in container.logs(stream=True):
#     print(line.strip().decode('utf-8'))
