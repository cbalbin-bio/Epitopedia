## Description of headers

This document contains the description for the headers used across various output files

### EPI-SEQ Headers

Header | Description
------------ | -------------
EPI_SEQ Input Structure | The structure provided to Epitopedia that this row relates to. Its sequence is extracted and used as a query against EPI_SEQ
EPI_SEQ Epitope ID | The Epitope ID for the subject in the resulting BLAST alignment
EPI_SEQ Input Structure Seq Start Pos | The starting position of the query sequence, extracted from the input structure, in the resulting BLAST Alignment (NOTE: due to inconsistencies in the PDB/mmCIF format, this is not necessarily the starting position in the Structure itself)
EPI_SEQ Input Structure Seq End Pos | Ending position for what is described above
EPI_SEQ Epitope Start Pos | The starting position for the subject (Epitope) in the resulting BLAST alignment
EPI_SEQ Epitope End Pos | Ending position for what is described above
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
EPI_SEQ Input Struc ASA | List of the ASA values as calculaed from DSSP for the aligned portion of the input structure
EPI_SEQ Input Struc Solv Seq | String where "?" describes unsolved residues in the aligned input structure (query) sequence
EPI_SEQ Input Struc Res Nums | List mapping the aligned input structure (query) sequence to the associated residue number, if present, in the PDB structure


### EPI-PDB Headers

Header | Description
------------ | -------------