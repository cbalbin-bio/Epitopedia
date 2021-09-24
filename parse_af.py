from collections import defaultdict
from io import TextIOWrapper
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


def extract_data(cif: str) -> dict[str, dict[str, list[str]]]:

    # Create defaultdict with default type of another defaultdict with list default
    data = defaultdict(lambda: defaultdict(list))

    doc = cif.read_file(cif)
    block = doc.sole_block()
    scheme = block.get_mmcif_category("_pdbx_poly_seq_scheme.", raw=True)

    # Ensure all columns contains same number of rows
    _lenths = [len(value) for value in scheme.values()]
    assert _lenths.count(_lenths[0]) == len(_lenths)

    # Split data on auth chains (NCBI still uses auth (stand_id) instead of asym)
    for index, pdb_strand_id in enumerate(scheme["pdb_strand_id"]):
        for key in scheme.keys():
            data[pdb_strand_id][key].append(scheme[key][index])

    return data


def extract_plddt(cif: str) -> list[str]:
    doc = cif.read_file(cif)
    block = doc.sole_block()
    scheme = block.get_mmcif_category("_ma_qa_metric_local.", raw=True)
    return scheme["metric_value"]


def write_cif_data_csv(data: dict[str, dict[str, list[str]]], handle: TextIOWrapper, pdb_id: str):
    for chain, chain_data in data.items():
        try:
            handle.write(
                f'{pdb_id}_{chain},{pdb_id}_{chain_data["asym_id"][0]},{"".join([protein_3to1[res] for res in chain_data["mon_id"]])},{"".join([protein_3to1[res] for res in chain_data["pdb_mon_id"]])},{" ".join(chain_data["auth_seq_num"])},{" ".join(chain_data["asym_id"])},{" ".join(chain_data["pdb_ins_code"])}'
            )
        except KeyError:
            continue


def write_cif_data_csv_af(data: dict[str, dict[str, list[str]]], plddt: list[str], handle: TextIOWrapper, pdb_id: str):
    for chain, chain_data in data.items():
        try:
            handle.write(
                f'{pdb_id}_{chain}, {"".join([protein_3to1[res] for res in chain_data["mon_id"]])},{" ".join(plddt)}'
            )
        except KeyError:
            continue
