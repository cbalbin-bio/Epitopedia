# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.


# _pdbx_poly_seq_scheme.asym_id asym_id used in mmcif
# _pdbx_poly_seq_scheme.mon_id sequences studied
# _pdbx_poly_seq_scheme.pdb_mon_id amino acids in structure studied
# _pdbx_poly_seq_scheme.pdb_strand_id  classic strand id

# auth seq num _pdbx_poly_seq_scheme.auth_seq_num


from dataclasses import dataclass, field
import sqlite3
import config


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


class PDBSeqs:
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


# def parseDSSPMMCIF(structure, path):
#     model = structure[0]
#     dssp = DSSP(model, path, file_type="MMCIF", dssp='mkdssp', acc_array='Wilke')
#     return dssp


# @dataclass
# class ChainData():
#     studied_seq : list
#     struc_seq : list
#     seq_nums : list
#     acc_seq : list[float] = field(default_factory=list)


# class MMCIFMeta():

#     def __init__(self, path, structure, compute_acc=False):
#         self.__mmcif_path__ = path
#         self.__mmcif_structure__ = structure
#         mmcif_dict = MMCIF2Dict(path)

#         self.__mon_seq__ = mmcif_dict["_pdbx_poly_seq_scheme.mon_id"]
#         self.__pdb_mon_seq__ = mmcif_dict["_pdbx_poly_seq_scheme.pdb_mon_id"]
#         self.__pdb_strand_id__ = mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"]
#         self.__auth_seq_num__ = mmcif_dict["_pdbx_poly_seq_scheme.auth_seq_num"]


#         self.__splitOnChain__()

#         if compute_acc:
#             self.__add_acc_data__()


#     def __splitOnChain__(self):
#         self.chain_dict = {}

#         split_mon_seq = []
#         split_pdb_mon_seq = []
#         split_seq_num = []

#         curr_chain_id = None
#         for chain, mon, pdb_mon, seq_num in zip(self.__pdb_strand_id__, self.__mon_seq__, self.__pdb_mon_seq__, self.__auth_seq_num__):
#             if curr_chain_id == None:
#                 curr_chain_id = chain
#                 split_mon_seq.append(mon)
#                 split_pdb_mon_seq.append(pdb_mon)
#                 split_seq_num.append(seq_num)
#             elif curr_chain_id == chain:
#                 split_mon_seq.append(mon)
#                 split_pdb_mon_seq.append(pdb_mon)
#                 split_seq_num.append(seq_num)

#             elif curr_chain_id != chain:
#                 self.chain_dict[curr_chain_id] = ChainData(split_mon_seq,split_pdb_mon_seq,split_seq_num )

#                 split_mon_seq = []
#                 split_pdb_mon_seq = []
#                 split_seq_num = []

#                 curr_chain_id = chain

#         self.chain_dict[curr_chain_id] = ChainData(split_mon_seq,split_pdb_mon_seq,split_seq_num)

#     def __add_acc_data__(self):
#         dssp = parseDSSPMMCIF(self.__mmcif_structure__, self.__mmcif_path__)

#         for chain, chain_data in self.chain_dict.items():
#             acc_seq = []
#             for res in chain_data.seq_nums:
#                 if res == "?":
#                     acc_seq.append(res)
#                 else:
#                     acc_seq.append(dssp[(chain, int(res))][3])

#             chain_data.acc_seq = acc_seq
