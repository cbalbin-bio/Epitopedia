from dataclasses import dataclass, field
from make_dice_gemmi import extract
import config

from MMCIFMeta import PDBSeqs, AccAgree
import os
import subprocess
from dataclasses_json import dataclass_json


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
    motif_seq: str = None
    motif_res_nums_query: list[int] = field(default_factory=list)
    motif_res_nums_target: list[int] = field(default_factory=list)
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


def hit_to_csv(file_path, motif, query_pdb_res_nums, motif_acc, pdbhit, epitope, hit):
    if not os.path.exists(file_path):
        with open(file_path, "a") as output_handle:
            output_handle.write(
                "motif\tquery_pdb_res_nums\tmotif_acc\tpdbhit_query\tpdbhit_target\tpdbhit_qcov\tpdbhit_pident\tpdbhit_evalue\tpdbhit_seqres\tpdbhit_seqsolv\tpdbhit_seqnums\tpdbhit_motif_seq\tpdbhit_motif_res_nums_query\tpdbhit_motif_res_nums_target\tpdbhit_query_structure_dice_path\tpdbhit_target_struc_dice_path\tpdbhit_TMalign_RMSD\tpdbhit_TMalign_TMscore\tpdbhit_TMalign_PDB_file\tpdbhit_TMalign_png_file\tpdbhit_query_acc_motif\tpdbhit_target_acc_motif\tpdbhit_query_perc_acc\tpdbhit_target_perc_acc\tpdbhit_perc_acc_agree\tepitope_epitope_id\tepitope_description\tepitope_linear_peptide_seq\tepitope_linear_peptide_modified_seq\tepitope_linear_peptide_modification\tepitope_source_antigen_accession\tepitope_starting_position\tepitope_ending_position\tepitope_database\tepitope_name\tepitope_organism_id\tepitope_organism_name\tepitope_sequence\tepitope_internal_source_seq_acc\thit_query_accession\thit_subject_accession\thit_query_start\thit_query_end\thit_subject_start\thit_subject_end\thit_aln_query_seq\thit_aln_subject_seq\thit_evalue\thit_qcovs\thit_pident\thit_staxid\thit_match_ranges\thit_cigar\thit_match_lengths\thit_submatch_seqs\thit_acc_seq\thit_pdb_seqsolv\thit_pdb_seqnums\n"
            )

    with open(file_path, "a") as output_handle:
        output_handle.write(
            f"{motif}\t{query_pdb_res_nums}\t{motif_acc}\t{pdbhit.query}\t{pdbhit.target}\t{pdbhit.qcov}\t{pdbhit.pident}\t{pdbhit.evalue}\t{pdbhit.seqres}\t{pdbhit.seqsolv}\t{pdbhit.seqnums}\t{pdbhit.motif_seq}\t{pdbhit.motif_res_nums_query}\t{pdbhit.motif_res_nums_target}\t{pdbhit.query_struc_dice_path}\t{pdbhit.target_struc_dice_path}\t{pdbhit.TMalign_RMSD}\t{pdbhit.TMalign_TMscore}\t{pdbhit.TMalign_PDB_file}\t{pdbhit.TMalign_png_file}\t{pdbhit.query_acc_motif}\t{pdbhit.target_acc_motif}\t{pdbhit.query_perc_acc}\t{pdbhit.target_perc_acc}\t{pdbhit.perc_acc_agree}\t{epitope.epitope_id}\t{epitope.description}\t{epitope.linear_peptide_seq}\t{epitope.linear_peptide_modified_seq}\t{epitope.linear_peptide_modification}\t{epitope.source_antigen_accession}\t{epitope.starting_position}\t{epitope.ending_position}\t{epitope.database}\t{epitope.name}\t{epitope.organism_id}\t{epitope.organism_name}\t{epitope.sequence}\t{epitope.internal_source_seq_acc}\t{hit.query_accession}\t{hit.subject_accession}\t{hit.query_start}\t{hit.query_end}\t{hit.subject_start}\t{hit.subject_end}\t{hit.aln_query_seq}\t{hit.aln_subject_seq}\t{hit.evalue}\t{hit.qcovs}\t{hit.pident}\t{hit.staxid}\t{hit.match_ranges}\t{hit.cigar}\t{hit.match_lengths}\t{hit.submatch_seqs}\t{hit.acc_seq}\t{hit.pdb_seqsolv}\t{hit.pdb_seqnums}\n"
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

        hit_to_csv(
            f"{config.OUTPUT_DIR}/EPI_PDB_hits_{pdb_input_str}.tsv",
            motif,
            query_pdb_res_nums,
            motif_acc,
            pdbhit,
            epitope,
            hit,
        )

        # look for motif from epitope blast in solv seq of pdb structure
        idx = pdbhit.seqsolv.find(motif)

        if idx != -1:

            # found a match

            # ensure query structure dice has been generated, else generate
            query_pdb_filename = f"{config.DICE_DIR}/{query_structure_basename}_{query_pdb_chain}_{query_pdb_res_nums[0]}-{query_pdb_res_nums[-1]}.pdb"
            if not os.path.exists(query_pdb_filename):

                extract(
                    f"{config.PDB_DATABASE_DIR}/{query_structure_basename[1:3]}/{query_structure_basename}.cif",
                    query_pdb_chain,
                    query_pdb_res_nums[0],
                    query_pdb_res_nums[-1],
                    query_pdb_filename,
                )

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

            target_base_pdb_name = pdbhit.target.split("_")[0].lower()
            target_chain_pdb_name = pdbhit.target.split("_")[1]

            # TODO: implement surface acc on target side
            #

            # you can get rid of this step later and just incoperate it into the pdb hit query against the sql db with another left join
            # doing it this way to save dev time
            try:
                target_pdb_seqs = PDBSeqs(target_base_pdb_name, target_chain_pdb_name, compute_acc=True)
            except:
                # print(target_base_pdb_name, target_chain_pdb_name)
                continue

            if target_pdb_seqs.binaryrasa:
                target_pdb_acc = target_pdb_seqs.binaryrasa[idx : idx + len(motif)]

                pdbhit.query_acc_motif = motif_acc
                pdbhit.target_acc_motif = target_pdb_acc
                if len(target_pdb_acc) == 0:
                    breakpoint()
                acc_agreement = AccAgree(motif_acc, target_pdb_acc)

                pdbhit.query_perc_acc = acc_agreement.percentAccQuery
                pdbhit.target_perc_acc = acc_agreement.percentAccTarget
                pdbhit.perc_acc_agree = acc_agreement.percentAccAgree

            # ensure target structure dice has been generated, else generate
            target_pdb_filename = f"{config.DICE_DIR}/{target_base_pdb_name}_{target_chain_pdb_name}_{pdbhit.motif_res_nums_target[0]}-{pdbhit.motif_res_nums_target[-1]}.pdb"
            if not os.path.exists(target_pdb_filename):

                extract(
                    f"{config.PDB_DATABASE_DIR}/{target_base_pdb_name[1:3]}/{target_base_pdb_name}.cif",
                    target_chain_pdb_name,
                    pdbhit.motif_res_nums_target[0],
                    pdbhit.motif_res_nums_target[-1],
                    target_pdb_filename,
                )

                # result = subprocess.run([f"pdb_selchain -{target_chain_pdb_name} {PDB_DATABASE_DIR}/{target_base_pdb_name[1:3]}/pdb{target_base_pdb_name}.ent | pdb_selres -{pdbhit.motif_res_nums_target[0]}:{pdbhit.motif_res_nums_target[-1]} > {target_pdb_filename}"],shell=True, capture_output=True)

                # if result.stderr:

                #     continue

            pdbhit.target_struc_dice_path = target_pdb_filename

            TMalign_prefix = f"{query_structure_basename}_{query_pdb_chain}_{query_pdb_res_nums[0]}-{query_pdb_res_nums[-1]}___{target_base_pdb_name}_{target_chain_pdb_name}_{pdbhit.motif_res_nums_target[0]}-{pdbhit.motif_res_nums_target[-1]}"

            # ensure TM align has been gernerated, else run
            if not os.path.exists(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout"):

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
                    # print(result.stderr)
                    continue
                with open(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout", "w") as output_handle:
                    output_handle.write(result.stdout)

            # read in tmalign output

            with open(f"{config.TMALIGN_DIR}/{TMalign_prefix}.stdout") as input_handle:

                for line in input_handle:
                    if line.startswith("Aligned length="):
                        pdbhit.TMalign_RMSD = line[30:34]  # line[30:34])
                    if line.startswith("TM-score= "):
                        pdbhit.TMalign_TMscore = line[10:17]  # float(line[10:17])
                        break

            # if pymol image of tmalign superposiiton has not been generated, generate it
            # if not os.path.exists(f"{TMALIGN_DIR}/{TMalign_prefix}_atm.pdb.png"):
            #     subprocess.run(["render-pymol.sh", f"{TMALIGN_DIR}/{TMalign_prefix}_atm.pdb"], capture_output=True)

            pdbhit.TMalign_PDB_file = f"{config.TMALIGN_DIR}/{TMalign_prefix}_atm.pdb"
            # pdbhit.TMalign_png_file = f"{TMALIGN_DIR}/{TMalign_prefix}_atm.pdb.png"

            ## remeber to remove later
            # if pdbhit.TMalign_TMscore == -1:
            #     breakpoint()
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
