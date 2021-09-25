# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import os
import subprocess
from Bio.PDB.DSSP import make_dssp_dict
from epitopedia.data.config import epitopedia_con

epitopedia_con

def runDSSP(path):
    basename = os.path.basename(path).split(".")[0]

    if not os.path.exists(f"/app/data/dssp/{basename}.dssp"):
        subprocess.run(f"mkdssp {path} > /app/data/dssp/{basename}.dssp", shell=True, capture_output=True)


def parseDSSP(path):

    basename = os.path.basename(path).split(".")[0]

    chains_acc = []

    # try:

    dssp_dict = make_dssp_dict(path)
    chains = {key[0] for key in dssp_dict[0].keys()}
    # except:
    #     return chains_acc

    for chain in chains:

        # try:
        cur = epitopedia_con.cursor()
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

        # except:
        #     continue

        chains_acc.append(f"{pdb_id},{chain_acc}\n")

    return chains_acc