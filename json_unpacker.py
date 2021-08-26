display = "SeqBMM Motif\tSeqBMM Input Struc Res Nums\tSeqBMM Acc\tEPI_PDB Epi Source Acc\tEPI_PDB Rep PDB\tEPI_PDB Qcov\tEPI_PDB Pident\tEPI_PDB Evalue\tmmCIF_SEQ Rep Res\tmmCIF_SEQ Rep Solv\tmmCIF_SEQ Rep Num\tEPI_PDB Rep Res Nums\tEPI_PDB Input Dice Path\tEPI_PDB Rep Dice Path\tEPI_PDB TMalign RMSD\tEPI_PDB TMalign TMscore\tEPI_PDB TMalign PDB\tEPI_PDB Rep Acc\tEPI_PDB Input Motif Perc Acc\tEPI_PDB Rep Motif Perc Acc\tEPI_PDB Perc Acc Agree\tIEDB_FILT Epitope ID\tIEDB_FILT Epitope Seq\tIEDB_FILT Source Seq Acc\tIEDB_FILT Start Pos\tIEDB_FILT Stop Pos\tIEDB_FILT Source DB\tIEDB_FILT Source Title\tIEDB_FILT Source Taxid\tIEDB_FILT Source Org\tIEDB_FILT Source Seq\tIEDB_FILT Iacc\tEPI_SEQ Input Structure\tEPI_SEQ Epitope ID\tEPI_SEQ Input Structure Seq Start Pos\tEPI_SEQ Input Structure Seq Stop Pos\tEPI_SEQ Epitope Start Pos\tEPI_SEQ Epitope End Pos\tEPI_SEQ Aln Input Struc Seq\tEPI_SEQ Aln Epitope Seq\tEPI_SEQ Evalue\tEPI_SEQ Qcov\tEPI_SEQ Pident\tEPI_SEQ Epitope Taxid\tEPI_SEQ Span Ranges\tEPI_SEQ Aln Cigar\tEPI_SEQ Span Lengths\tEPI_SEQ Span Seqs\tPDB_DSSP Input Struc ASA\tmmCIF_SEQ Input Struc Solv Seq\tmmCIF_SEQ Input Struc Res Nums\n"
keys = "{pdbhit.motif_seq}\t{pdbhit.motif_res_nums_query}\t{pdbhit.query_acc_motif}\t{pdbhit.query}\t{pdbhit.target}\t{pdbhit.qcov}\t{pdbhit.pident}\t{pdbhit.evalue}\t{pdbhit.seqres}\t{pdbhit.seqsolv}\t{pdbhit.seqnums}\t{pdbhit.motif_res_nums_target}\t{pdbhit.query_struc_dice_path}\t{pdbhit.target_struc_dice_path}\t{pdbhit.TMalign_RMSD}\t{pdbhit.TMalign_TMscore}\t{pdbhit.TMalign_PDB_file}\t{pdbhit.target_acc_motif}\t{pdbhit.query_perc_acc}\t{pdbhit.target_perc_acc}\t{pdbhit.perc_acc_agree}\t{epitope.epitope_id}\t{epitope.linear_peptide_seq}\t{epitope.source_antigen_accession}\t{epitope.starting_position}\t{epitope.ending_position}\t{epitope.database}\t{epitope.name}\t{epitope.organism_id}\t{epitope.organism_name}\t{epitope.sequence}\t{epitope.internal_source_seq_acc}\t{hit.query_accession}\t{hit.subject_accession}\t{hit.query_start}\t{hit.query_end}\t{hit.subject_start}\t{hit.subject_end}\t{hit.aln_query_seq}\t{hit.aln_subject_seq}\t{hit.evalue}\t{hit.qcovs}\t{hit.pident}\t{hit.staxid}\t{hit.match_ranges}\t{hit.cigar}\t{hit.match_lengths}\t{hit.submatch_seqs}\t{hit.acc_seq}\t{hit.pdb_seqsolv}\t{hit.pdb_seqnums}\n"

display = display.strip().split("\t")
keys = keys.strip().replace("{", "").replace("}", "").split("\t")


key_dict = {
    "hit": "hit",
    "epitope": "epitope",
    "pdbhit": "motif_type_pdb_hit",
}

with open("unpack.txt", "w") as output:
    for key, dis in zip(keys, display):
        key = key.split(".")

        new_key = key_dict[key[0]]

        output.write(f'"{dis}" : {new_key}["{key[1]}"],\n')
