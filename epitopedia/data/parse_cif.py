# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

from collections import defaultdict
from io import TextIOWrapper
from typing import Union

from gemmi import cif

protein_3to1 = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
    "?": "?",
}


def extract_data(path: str, use_auth: bool = True) -> dict[str, dict[str, list[str]]]:

    # Create defaultdict with default type of another defaultdict with list default
    data = defaultdict(lambda: defaultdict(list))

    doc = cif.read_file(path)
    block = doc.sole_block()
    scheme = block.get_mmcif_category("_pdbx_poly_seq_scheme.", raw=True)

    # Ensure all columns contains same number of rows
    _lenths = [len(value) for value in scheme.values()]
    assert _lenths.count(_lenths[0]) == len(_lenths)

    # Split data on auth chains (NCBI still uses auth (stand_id) instead of asym)
    if use_auth:
        for index, pdb_strand_id in enumerate(scheme["pdb_strand_id"]):
            for key in scheme.keys():
                data[pdb_strand_id][key].append(scheme[key][index])
    else:
        for index, asym_id in enumerate(scheme["asym_id"]):
            for key in scheme.keys():
                data[asym_id][key].append(scheme[key][index])

    return data


def extract_plddt(path: str) -> tuple[list[str], str]:
    doc = cif.read_file(path)
    block = doc.sole_block()
    return list(block.find_values("_ma_qa_metric_local.metric_value")), block.find_value(
        "_ma_qa_metric_global.metric_value"
    )


def write_cif_data_csv(
    data: dict[str, dict[str, list[str]]],
    handle: TextIOWrapper,
    pdb_id: str,
    plddt: Union[bool, tuple[list[str], str]] = False,
):
    if plddt:
        # writing af prediction

        try:
            handle.write(
                f'{pdb_id}_A,{pdb_id}_A,{"".join([protein_3to1[res] for res in data["A"]["mon_id"]])},{"".join([protein_3to1[res] for res in data["A"]["mon_id"]])},{" ".join(data["A"]["auth_seq_num"])},{" ".join(data["A"]["asym_id"])},{" ".join(data["A"]["pdb_ins_code"])},{" ".join(plddt[0])},{plddt[1]},TRUE\n'
            )
        except KeyError:
            return
    else:

        for chain, chain_data in data.items():
            try:
                handle.write(
                    f'{pdb_id}_{chain},{pdb_id}_{chain_data["asym_id"][0]},{"".join([protein_3to1[res] for res in chain_data["mon_id"]])},{"".join([protein_3to1[res] for res in chain_data["pdb_mon_id"]])},{" ".join(chain_data["auth_seq_num"])},{" ".join(chain_data["asym_id"])},{" ".join(chain_data["pdb_ins_code"])},NULL,NULL,FALSE\n'
                )
            except KeyError:
                continue
