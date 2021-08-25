# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import gemmi


def extract(input, chain_name, start, stop, output):
    w_st = gemmi.Structure()
    w_m = w_st.add_model(gemmi.Model("0"), 0)
    w_c = w_m.add_chain("A")

    start = int(start)
    stop = int(stop)

    r_st = gemmi.read_structure(input)
    for chain in r_st[0]:
        if chain.name == chain_name:
            for res in chain:
                if start <= res.seqid.num <= stop:
                    w_c.add_residue(res)

    # w_st.make_mmcif_document().write_file(output)
    w_st.write_minimal_pdb(output)
