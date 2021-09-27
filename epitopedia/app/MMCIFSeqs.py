# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import sqlite3
from dataclasses import dataclass, field

from epitopedia.app import config


class AccAgree:
    def __init__(self, query, target, acc_char="A"):
        query_len = len(query)
        target_len = len(target)

        assert query_len == target_len, "Cant compute acc agreement for sequences of different length"

        query_acc_count = 0
        target_acc_count = 0
        agree_count = 0
        for query_char, target_char in zip(query, target):
            if query_char == acc_char:
                query_acc_count += 1
            if target_char == acc_char:
                target_acc_count += 1

            if query_char == target_char:
                agree_count += 1

        self.percentAccQuery = query_acc_count / query_len
        self.percentAccTarget = target_acc_count / target_len
        self.percentAccAgree = agree_count / target_len


class MMCIFSeqs:
    def __init__(self, pdbid, chain, compute_acc=False, threshold=0.25):
        con = sqlite3.connect(config.SQLITE_DATABASE_DIR)
        cur = con.cursor()
        cur.execute(f'SELECT seqres, seqsolv, seqnums FROM mmCIF_seqs WHERE pdb_id = "{pdbid.lower()}_{chain}"')
        row = cur.fetchone()

        self.seqres = row[0]
        self.seqsolv = row[1]
        self.seqnums = row[2].split(" ")

        if compute_acc:
            cur.execute(f'SELECT acc FROM PDB_DSSP WHERE pdb_id = "{pdbid.lower()}_{chain}"')
            row = cur.fetchone()
            if row == None:
                self.acc = None
                self.rasa = None
                self.binaryrasa = None

            else:
                self.acc = row[0].split(" ")

                self.rasa = []
                for val, res in zip(self.acc, self.seqsolv):
                    if res == "?":
                        self.rasa.append("?")
                    else:
                        self.rasa.append(float(val) / config.residue_max_acc["Wilke"][res])

                self.binaryrasa = []
                for val in self.rasa:
                    if val == "?":
                        self.binaryrasa.append("?")
                    elif val >= threshold:
                        self.binaryrasa.append("A")
                    else:
                        self.binaryrasa.append("B")

        else:
            self.acc = None
            self.rasa = None
            self.binaryrasa = None
