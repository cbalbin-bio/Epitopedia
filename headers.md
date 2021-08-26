## Description of Headers

This document contains the description for the headers used across various output files

### EPI-SEQ Headers

Header | Description
------------ | -------------
EPI_SEQ Input Structure | The structure provided to Epitopedia that this row relates to. Its sequence is extracted and used as a query against EPI_SEQ
EPI_SEQ Epitope ID | The Epitope ID for the subject in the resulting BLAST alignment
EPI_SEQ Input Structure Seq Start Pos | The starting position of the query sequence, extracted from the input structure, in the resulting BLAST Alignment (NOTE: due to inconsistencies in the PDB/mmCIF format, this is not necessarily the starting position in the Structure itself)
EPI_SEQ Input Structure Seq End Pos | Ending position for what is described above
EPI_SEQ Epitope Start Pos | The starting position for the subject (Epitope) in the resulting BLAST alignment
EPI_SEQ Epitope Stop Pos | Ending position for what is described above
EPI_SEQ Aln Input Struc Seq | The query portion of the aligned sequence resulting from the BLAST alignment
EPI_SEQ Aln Epitope Seq | The subject portion of the aligned sequence resulting from the BLAST alignment
EPI_SEQ Evalue | The Evalue for the resulting BLAST alignment
EPI_SEQ Qcov | The query coverage for the resulting BLAST alignment
EPI_SEQ Pident | The percent identity for the resulting BLAST alignment
EPI_SEQ Epitope Taxid | The taxid for the organism of which the epitope is sourced from
EPI_SEQ Aln Cigar | the cigar sequence for the BLAST alignment
EPI_SEQ Span Ranges | List of lists containing the starting and ending position for spans of consecutive matches in the BLAST alignment
EPI_SEQ Span Lengths | List of the lengths of the EPI_SEQ Span Ranges
EPI_SEQ Span Seqs | List of amino acid sequences resulting from these spans
PDB_DSSP Input Struc ASA | List of the ASA values as calculaed from DSSP for the aligned portion of the input structure
mmCIF_SEQ Input Struc Solv Seq | String where "?" describes unsolved residues in the aligned input structure (query) sequence
mmCIF_SEQ Input Struc Res Nums | List mapping the aligned input structure (query) sequence to the associated residue number, if present, in the PDB structure


### EPI-PDB Headers

Header | Description
------------ | -------------
SeqBMM Motif | motif of the SeqBMM found in the EPI_SEQ search step
SeqBMM Input Struc Res Nums | List containing the residue numbers in the input structure that correspond to the reidues in the SeqBMM motif
SeqBMM Acc | List contianing the surface accessibility classification for each residue in the SeqBMM motif; A = accessible, B = Buried
EPI_PDB Epi Source Acc | The epitope soruce sequence accession which is used as a query against EPI_PDB
EPI_PDB Rep PDB | The PDB_ID for the resulting representative strucure (target) from the query against EPI_PDB
EPI_PDB Qcov | The query coverage for the mmSEQ result
EPI_PDB Pident | The percent identity for the mmSEQ result
EPI_PDB Evalue | The Evalue for the mmSEQ result
mmCIF_SEQ Rep Res | The residue sequence of the structural representative
mmCIF_SEQ Rep Solv | String where "?" describes unsolved residues in the structural representative
mmCIF_SEQ Rep Num | List mapping the sequence of the structural representative to the associated residue number, if present, in the PDB structure
EPI_PDB Rep Res Nums | List containing the residue numbers in the structural representative that correspond to the reidues in the EPI_PDB motif
EPI_PDB Input Dice Path | Path to the SeqBMM motif dice from the input structure 
EPI_PDB Rep Dice Path | Path to the SeqBMM motif dice from the representative srructure
EPI_PDB TMalign RMSD | resulting RMSD value for the alignment of the dices described above
EPI_PDB TMalign TMscore | TMscore for above
EPI_PDB TMalign PDB | Path to the PDB file containing the structural superposition of the two dices described above
EPI_PDB Rep Acc | List contianing the surface accessibility classification for each residue in the representative motif; A = accessible, B = Buried
EPI_PDB Input Motif Perc Acc | Percentage of accessible residue in the motif span of the input structure
EPI_PDB Rep Motif Perc Acc | Percentage of accessible residue in the motif span of the representative structure
EPI_PDB Perc Acc Agree | Percentage of surface accessibility agreement between the motif span of the input structure and the motif span of the representative structure
IEDB_FILT Epitope ID | Epitope ID of the epitope which contains the SeqBMM
IEDB_FILT Epitope Seq | Full sequence of the Epitope
IEDB_FILT Source Seq Acc | Acession for the source sequence of which the epitope is a subsequence of
IEDB_FILT Start Pos | Starting position of the epitope in the epitope source sequence
IEDB_FILT Stop Pos | Stoping position for above
IEDB_FILT Source DB | Database from which the source sequence was obtained
IEDB_FILT Source Title | Title for the source sequence protein
IEDB_FILT Source Taxid | Taxid for the organism associated with the source sequence
IEDB_FILT Source Org | Organism name the source sequence is associated with
IEDB_FILT Source Seq | Epitope Source sequence
IEDB_FILT Iacc | Internal acession number

[For remaining headers in this file see EPI-SEQ headers](###EPI-SEQ Headers)


