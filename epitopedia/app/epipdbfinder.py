# Copyright (c) 2021 Christian Balbin
# This work is licensed under the terms of the MIT license.
# For a copy, see <https://opensource.org/licenses/MIT>.

import os
import subprocess
from dataclasses import dataclass, field

from dataclasses_json import dataclass_json
from epitopedia.app import config
from epitopedia.app.config import console
from epitopedia.utils.make_dice_gemmi import extract
from epitopedia.utils.utils import obtain_lock, release_lock, wait_unlock
from epitopedia.app.MMCIFSeqs import AccAgree, MMCIFSeqs
from epitopedia.app.args import args


@dataclass_json
@dataclass
class PDBHit:
    query: str
    target: str
    qcov: float
    pident: float
    evalue: float
    seqres: str
    seqsolv: str
    seqnums: str
    lplddt: str
    gplddt: float
    isAF: bool
    pdb_title: str
    pdb_species: str
    motif_seq: str = None
    motif_res_nums_query: list[int] = field(default_factory=list)
    motif_res_nums_target: list[int] = field(default_factory=list)
    motif_lplddt: list[int] = field(default_factory=list)
    avg_motif_lplddt: float = -1
    query_struc_dice_path: str = ""
    target_struc_dice_path: str = ""
    TMalign_RMSD: float = -1.0
    TMalign_TMscore: float = -1.0
    TMalign_PDB_file: str = ""
    TMalign_png_file: str = ""
    query_acc_motif: str = ""
    target_acc_motif: str = ""
    query_perc_acc: float = -1
    target_perc_acc: float = -1
    perc_acc_agree: float = -1

    def __post_init__(self):
        self.qcov = float(self.qcov)
        self.pident = float(self.pident)
        self.evalue = float(self.evalue)

        if self.isAF == "TRUE":
            self.isAF = True
            self.gplddt = float(self.gplddt)
        else:
            self.isAF = False


def hit_to_csv(file_path, motif, query_pdb_res_nums, motif_acc, pdbhit, epitope, hit):
    if not os.path.exists(file_path):
        with open(file_path, "a") as output_handle:
            output_handle.write(
                "SeqBMM Motif\tSeqBMM Input Struc Res Nums\tSeqBMM Acc\tEPI_PDB Epi Source Acc\tEPI_PDB Rep PDB\tEPI_PDB Qcov\tEPI_PDB Pident\tEPI_PDB Evalue\tmmCIF_SEQ Rep Res\tmmCIF_SEQ Rep Solv\tmmCIF_SEQ Rep Num\tEPI_PDB Rep Res Nums\tEPI_PDB Input Dice Path\tEPI_PDB Rep Dice Path\tEPI_PDB TMalign RMSD\tEPI_PDB TMalign TMscore\tEPI_PDB TMalign PDB\tEPI_PDB Rep Acc\tEPI_PDB Input Motif Perc Acc\tEPI_PDB Rep Motif Perc Acc\tEPI_PDB Perc Acc Agree\tIEDB_FILT Epitope ID\tIEDB_FILT Epitope Seq\tIEDB_FILT Source Seq Acc\tIEDB_FILT Start Pos\tIEDB_FILT Stop Pos\tIEDB_FILT Source DB\tIEDB_FILT Source Title\tIEDB_FILT Source Taxid\tIEDB_FILT Source Org\tIEDB_FILT Source Seq\tIEDB_FILT Iacc\tEPI_SEQ Input Structure\tEPI_SEQ Epitope ID\tEPI_SEQ Input Structure Seq Start Pos\tEPI_SEQ Input Structure Seq Stop Pos\tEPI_SEQ Epitope Start Pos\tEPI_SEQ Epitope End Pos\tEPI_SEQ Aln Input Struc Seq\tEPI_SEQ Aln Epitope Seq\tEPI_SEQ Evalue\tEPI_SEQ Qcov\tEPI_SEQ Pident\tEPI_SEQ Epitope Taxid\tEPI_SEQ Span Ranges\tEPI_SEQ Aln Cigar\tEPI_SEQ Span Lengths\tEPI_SEQ Span Seqs\tPDB_DSSP Input Struc ASA\tmmCIF_SEQ Input Struc Solv Seq\tmmCIF_SEQ Input Struc Res Nums\n"
            )

    with open(file_path, "a") as output_handle:
        output_handle.write(
            f"{motif}\t{query_pdb_res_nums}\t{motif_acc}\t{pdbhit.query}\t{pdbhit.target}\t{pdbhit.qcov}\t{pdbhit.pident}\t{pdbhit.evalue}\t{pdbhit.seqres}\t{pdbhit.seqsolv}\t{pdbhit.seqnums}\t{pdbhit.motif_res_nums_target}\t{pdbhit.query_struc_dice_path}\t{pdbhit.target_struc_dice_path}\t{pdbhit.TMalign_RMSD}\t{pdbhit.TMalign_TMscore}\t{pdbhit.TMalign_PDB_file}\t{pdbhit.target_acc_motif}\t{pdbhit.query_perc_acc}\t{pdbhit.target_perc_acc}\t{pdbhit.perc_acc_agree}\t{epitope.epitope_id}\t{epitope.linear_peptide_seq}\t{epitope.source_antigen_accession}\t{epitope.starting_position}\t{epitope.ending_position}\t{epitope.database}\t{epitope.name}\t{epitope.organism_id}\t{epitope.organism_name}\t{epitope.sequence}\t{epitope.internal_source_seq_acc}\t{hit.query_accession}\t{hit.subject_accession}\t{hit.query_start}\t{hit.query_end}\t{hit.subject_start}\t{hit.subject_end}\t{hit.aln_query_seq}\t{hit.aln_subject_seq}\t{hit.evalue}\t{hit.qcovs}\t{hit.pident}\t{hit.staxid}\t{hit.match_ranges}\t{hit.cigar}\t{hit.match_lengths}\t{hit.submatch_seqs}\t{hit.acc_seq}\t{hit.pdb_seqsolv}\t{hit.pdb_seqnums}\n"
        )


def hit_to_pdb(
    motif,
    query_pdb_res_nums,
    motif_acc,
    query_structure_basename,
    query_pdb_chain,
    mmseqs_db_rows,
    epitope,
    hit,
    pdb_input_str,
    min_pident=95.0,
):

    pdbhits = []
    row_count = 0

    # if motif == "TQLPP":
    #     breakpoint()

    for row in mmseqs_db_rows:
        # if row_count >= 10:
        #     break

        pdbhit = PDBHit(*row)
        if pdbhit.pident < min_pident:
            continue

        if pdbhit.isAF and not args.use_afdb:
            continue

        if pdbhit.isAF and pdbhit.gplddt < args.gplddt:
            continue

        hit_to_csv(
            f"{config.OUTPUT_DIR}/EPI_PDB_hits_{pdb_input_str}.tsv",
            motif,
            query_pdb_res_nums,
            motif_acc,
            pdbhit,
            epitope,
            hit,
        )

        idx = pdbhit.seqsolv.find(motif)

        if idx != -1:

            # found a match

            # ensure query structure dice has been generated, else generate
            query_pdb_filename = f"{config.DICE_DIR}/{query_structure_basename}_{query_pdb_chain}_{query_pdb_res_nums[0]}-{query_pdb_res_nums[-1]}.pdb"

            # if file doesnt eixist and not locked
            # obtain lock
            # generate file
            # release lock
            # if file is locked
            # sleep until it exists and isnt locked

            if not os.path.exists(query_pdb_filename) and obtain_lock(query_pdb_filename):
                try:
                    extract(
                        f"{config.PDB_DATABASE_DIR}/{query_structure_basename[1:3]}/{query_structure_basename}.cif",
                        query_pdb_chain,
                        query_pdb_res_nums[0],
                        query_pdb_res_nums[-1],
                        query_pdb_filename,
                    )
                    release_lock(query_pdb_filename)
                except:
                    release_lock(query_pdb_filename)
                    console.log(f"[bold red] {query_structure_basename}.cif not Found")
                    continue
            else:
                wait_unlock(query_pdb_filename)
                # result = subprocess.run([f"pdb_selchain -{query_pdb_chain} {PDB_DATABASE_DIR}/{query_structure_basename[1:3]}/pdb{query_structure_basename}.ent | pdb_selres -{query_pdb_res_nums[0]}:{query_pdb_res_nums[-1]} > {query_pdb_filename}"],shell=True, capture_output=True)

                # if result.stderr:
                #     continue

            pdbhit.motif_res_nums_query = query_pdb_res_nums
            pdbhit.query_struc_dice_path = query_pdb_filename

            # slice out the motif sequence from the pdb strcutre
            pdbhit.motif_seq = pdbhit.seqsolv[idx : idx + len(motif)]

            # make sure they are the same
            assert pdbhit.motif_seq == motif, "Mismatch between motif from blast and motif extracted from PDB sequence"

            # get the residue numbers that correspond to the motif in the pdb structure
            pdbhit.motif_res_nums_target = pdbhit.seqnums.split(" ")[idx : idx + len(motif)]

            if pdbhit.isAF:
                pdbhit.motif_lplddt = pdbhit.lplddt.split(" ")[idx : idx + len(motif)]
                pdbhit.motif_lplddt = [float(val) for val in pdbhit.motif_lplddt]
                pdbhit.avg_motif_lplddt = (sum(pdbhit.motif_lplddt) / len(pdbhit.motif_lplddt)) /100

                if pdbhit.avg_motif_lplddt < args.lplddt:
                    continue

            target_base_pdb_name = pdbhit.target.rsplit("_", 1)[0]
            target_chain_pdb_name = pdbhit.target.rsplit("_", 1)[1]

            # TODO: implement surface acc on target side
            #

            # you can get rid of this step later and just incorporate it into the pdb hit query against the sql db with another left join
            # doing it this way to save dev time
            try:
                target_pdb_seqs = MMCIFSeqs(target_base_pdb_name, target_chain_pdb_name, compute_acc=True)
            except:
                # print(target_base_pdb_name, target_chain_pdb_name)
                continue

            if target_pdb_seqs.binaryrasa:
                target_pdb_acc = target_pdb_seqs.binaryrasa[idx : idx + len(motif)]

                pdbhit.query_acc_motif = motif_acc
                pdbhit.target_acc_motif = target_pdb_acc
                acc_agreement = AccAgree(motif_acc, target_pdb_acc)

                pdbhit.query_perc_acc = acc_agreement.percentAccQuery
                pdbhit.target_perc_acc = acc_agreement.percentAccTarget
                pdbhit.perc_acc_agree = acc_agreement.percentAccAgree

            # ensure target structure dice has been generated, else generate
            target_pdb_filename = f"{config.DICE_DIR}/{target_base_pdb_name}_{target_chain_pdb_name}_{pdbhit.motif_res_nums_target[0]}-{pdbhit.motif_res_nums_target[-1]}.pdb"
            if not os.path.exists(target_pdb_filename) and obtain_lock(target_pdb_filename):
                try:
                    extract(
                        f"{config.AFDB_DIR}/{target_base_pdb_name}.cif"
                        if pdbhit.isAF
                        else f"{config.PDB_DATABASE_DIR}/{target_base_pdb_name[1:3]}/{target_base_pdb_name}.cif",
                        target_chain_pdb_name,
                        pdbhit.motif_res_nums_target[0],
                        pdbhit.motif_res_nums_target[-1],
                        target_pdb_filename,
                    )
                    release_lock(target_pdb_filename)
                except:
                    release_lock(target_pdb_filename)
                    console.log(f"[bold red] {target_base_pdb_name}.cif not Found")
                    continue
            else:
                wait_unlock(target_pdb_filename)

            pdbhit.target_struc_dice_path = target_pdb_filename

            TMalign_prefix = f"{query_structure_basename}_{query_pdb_chain}_{query_pdb_res_nums[0]}-{query_pdb_res_nums[-1]}___{target_base_pdb_name}_{target_chain_pdb_name}_{pdbhit.motif_res_nums_target[0]}-{pdbhit.motif_res_nums_target[-1]}"

            # ensure TM align has been gernerated, else run
            if not os.path.exists(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout") and obtain_lock(
                f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout"
            ):

                # force TM to properly align using  motif sequence alignment
                with open(f"{config.TMALIGN_DIR}/{TMalign_prefix}.fa", "w") as tmfasta:
                    tmfasta.write(
                        f">{query_structure_basename}\n{motif}\n>{target_base_pdb_name}\n{pdbhit.motif_seq}\n"
                    )

                # run tmalign, and dump output to txt file
                result = subprocess.run(
                    [
                        "TMalign",
                        query_pdb_filename,
                        target_pdb_filename,
                        "-I",
                        f"{config.TMALIGN_DIR}/{TMalign_prefix}.fa",
                        "-o",
                        f"{config.TMALIGN_DIR}/{TMalign_prefix}",
                    ],
                    capture_output=True,
                    text=True,
                )

                if result.stderr:
                    release_lock(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout")
                    # print(result.stderr)
                    continue
                with open(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout", "w") as output_handle:
                    output_handle.write(result.stdout)
                release_lock(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout")
            else:
                wait_unlock(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout")

            # read in tmalign output

            with open(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout") as input_handle:

                for line in input_handle:
                    if line.startswith("Aligned length="):
                        pdbhit.TMalign_RMSD = line[30:34]  # line[30:34])
                    if line.startswith("TM-score= "):
                        pdbhit.TMalign_TMscore = line[10:17]  # float(line[10:17])
                        break

            pdbhit.TMalign_PDB_file = f"{config.TMALIGN_DIR}/{TMalign_prefix}_atm.pdb"

            if args.rmsd:
                if float(pdbhit.TMalign_RMSD) >= args.rmsd:
                    continue
            hit_to_csv(
                f"{config.OUTPUT_DIR}/EPI_PDB_fragment_pairs_{pdb_input_str}.tsv",
                motif,
                query_pdb_res_nums,
                motif_acc,
                pdbhit,
                epitope,
                hit,
            )
            pdbhits.append(pdbhit)
            row_count += 1

    return pdbhits
