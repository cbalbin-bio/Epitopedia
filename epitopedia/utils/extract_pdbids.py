# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import sqlite3
import sys

con = sqlite3.connect(sys.argv[1])
cur = con.cursor()
pdbids = []

for row in cur.execute("SELECT DISTINCT target FROM EPI_PDB"):
    pdbids.append(row[0].split("_")[0].lower())


pdbids = set(pdbids)
with open("pdb_id_list.txt", "w") as handle:
    for pdbid in pdbids:

        handle.write(f"""*/{pdbid}.cif.gz\n""")
con.close()
