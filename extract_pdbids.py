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

        # handle.write(f"rsync.rcsb.org::ftp_data/structures/divided/mmCIF/{pdbid[1:3]}/{pdbid}.cif.gz\n")

        handle.write(f"""*/{pdbid}.cif.gz\n""")
con.close()
