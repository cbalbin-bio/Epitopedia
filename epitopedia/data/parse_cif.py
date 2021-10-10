# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

from collections import defaultdict
from io import TextIOWrapper
from posix import times_result
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


def extract_data(
    path: str, use_auth: bool = True, af: bool = False
) -> tuple[dict[str, dict[str, list[str]]], str, str]:

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

    if af:
        species = (
            list(block.find_values("_ma_target_ref_db_details.organism_scientific"))[0]
            .replace("\t", "")
            .replace("\n", "")
            .replace("'", "")
            .replace('"', "")
            .replace(";", "")
        )
        title = (
            list(block.find_values("_entity.pdbx_description"))[0]
            .replace("\t", "")
            .replace("\n", "")
            .replace("'", "")
            .replace('"', "")
            .replace(";", "")
        )
    else:
        species = (
            ", ".join(set(block.find_values("_entity_src_gen.pdbx_gene_src_scientific_name")))
            .replace("\t", "")
            .replace("\n", "")
            .replace("'", "")
            .replace('"', "")
            .replace(";", "")
        )
        title = (
            list(block.find_values("_struct.title"))[0]
            .replace("\t", "")
            .replace("\n", "")
            .replace("'", "")
            .replace('"', "")
            .replace(";", "")
        )

    return data, species, title


def extract_plddt(path: str) -> tuple[list[str], str]:
    doc = cif.read_file(path)
    block = doc.sole_block()
    return list(block.find_values("_ma_qa_metric_local.metric_value")), block.find_value(
        "_ma_qa_metric_global.metric_value"
    )


def write_cif_data_csv(
    data_s_t: tuple[dict[str, dict[str, list[str]]], str, str],
    handle: TextIOWrapper,
    pdb_id: str,
    plddt: Union[bool, tuple[list[str], str]] = False,
):

    data = data_s_t[0]
    species = data_s_t[1]
    title = data_s_t[2]
    if plddt:
        # writing af prediction

        try:
            handle.write(
                f'{pdb_id}_A\t{pdb_id}_A\t{"".join([protein_3to1[res] for res in data["A"]["mon_id"]])}\t{"".join([protein_3to1[res] for res in data["A"]["mon_id"]])}\t{" ".join(data["A"]["auth_seq_num"])}\t{" ".join(data["A"]["asym_id"])}\t{" ".join(data["A"]["pdb_ins_code"])}\t{" ".join(plddt[0])}\t{plddt[1]}\tTRUE\t{title}\t{species}\n'
            )
        except KeyError:
            return
    else:

        for chain, chain_data in data.items():
            try:
                handle.write(
                    f'{pdb_id}_{chain}\t{pdb_id}_{chain_data["asym_id"][0]}\t{"".join([protein_3to1[res] for res in chain_data["mon_id"]])}\t{"".join([protein_3to1[res] for res in chain_data["pdb_mon_id"]])}\t{" ".join(chain_data["auth_seq_num"])}\t{" ".join(chain_data["asym_id"])}\t{" ".join(chain_data["pdb_ins_code"])}\t\t\tFALSE\t{title}\t{species}\n'
                )
            except KeyError:
                continue
