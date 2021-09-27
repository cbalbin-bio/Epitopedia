# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import os
import time
from pathlib import Path


def remove_previous_files(config, pdb_inputs):

    try:
        os.remove(f"{config.OUTPUT_DIR}/EPI_SEQ_hits_{pdb_inputs}.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_SEQ_span_filt_hits_{pdb_inputs}.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_SEQ_span_filt_acc_hits_{pdb_inputs}.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_PDB_hits_{pdb_inputs}.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_inputs}.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_inputs}_ranked.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_inputs}_best.tsv")
        os.remove(f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_inputs}_best.json")
    except OSError:
        pass


def obtain_lock(resource):
    locked = False
    while not locked:
        try:
            Path(resource + ".lock").touch(exist_ok=False)
            locked = True
        except FileExistsError:
            time.sleep(0.1)


def release_lock(resource):
    os.remove(resource + ".lock")
