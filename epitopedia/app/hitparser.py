# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import sqlite3
from dataclasses import dataclass

from dataclasses_json import dataclass_json
from epitopedia.app import config
from epitopedia.app.epipdbfinder import hit_to_pdb


@dataclass_json
@dataclass
class Epitope:
    epitope_id: int
    description: str
    linear_peptide_seq: str
    linear_peptide_modified_seq: str
    linear_peptide_modification: str
    source_antigen_accession: str
    starting_position: int
    ending_position: int
    database: str
    name: str
    organism_id: int
    organism_name: str
    sequence: str
    internal_source_seq_acc: int


def parseHit(hit, span, pdb_seq, query_pdb_base, query_pdb_chain, pdb_input_str):
    # looping through each hit in the inital blast search of the query structure sequence against the short IEDB epitope sequences

    con = sqlite3.connect(config.SQLITE_DATABASE_DIR)
    cur = con.cursor()
    motifs = []
    motif_pdb_ranges = []
    motifs_acc = []

    for index, (start, end) in enumerate(hit.match_ranges):
        # convert to 0 based indexing, check to make sure that motif is at least 5 residues long
        if end - (start - 1) >= span:
            # extract the motif seq from the pdb structure used as a query against iedb
            query_motif = pdb_seq.seqsolv[start - 1 : end]
            # make sure there are not any disordered residues in the motif seq
            if "?" not in query_motif:
                # add it to the list of motifs to search for representative structures using the source epitope sequence against pdb later on, also add the corresponding seq res nums
                motifs.append(query_motif)
                motif_pdb_ranges.append(pdb_seq.seqnums[start - 1 : end])
                motifs_acc.append(pdb_seq.binaryrasa[start - 1 : end])

    # get epitope data on this hit against the iedb db

    cur.execute(f"SELECT * FROM IEDB_FILT WHERE epitope_id={hit.subject_accession}")
    ret = cur.fetchone()
    epitope = Epitope(*ret)

    # get mmseqs results for this hits source sequence against pdb, looking for representative structures

    cur.execute(
        f'SELECT query, target, qcov, pident, evalue, seqres, seqsolv, seqnums, lplddt, gplddt, AF, title, species FROM EPI_PDB LEFT JOIN mmCIF_seqs ON EPI_PDB.target = mmCIF_seqs.pdb_id WHERE EPI_PDB.query = "{epitope.source_antigen_accession}";'
    )
    mmseq_db_rows = cur.fetchall()

    pdb_hits_motifs = []
    # for each of the sub motifs in this  upper level hit
    for motif, motif_pdb_range, motif_acc in zip(motifs, motif_pdb_ranges, motifs_acc):
        # map it to the representative pdb strcutrues

        pdb_hits = hit_to_pdb(
            motif,
            motif_pdb_range,
            motif_acc,
            query_pdb_base,
            query_pdb_chain,
            mmseq_db_rows,
            epitope,
            hit,
            pdb_input_str,
        )
        if pdb_hits:
            # if data is returned ( succefully mapped to representative pdb structures), sort it by TMscore
            pdb_hits_motifs.append(sorted(pdb_hits, key=lambda x: x.TMalign_TMscore, reverse=True))

    if pdb_hits_motifs:

        return {
            "blasthit": hit.to_dict(),
            "hitepitopedata": epitope.to_dict(),
            "pdb_motif_hits": [[motif_hit.to_dict() for motif_hit in motif] for motif in pdb_hits_motifs],
        }
