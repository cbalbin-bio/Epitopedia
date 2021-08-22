import json
import os
from collections import defaultdict
import sys
import csv
import config


def reduce_results(path):
    with open(path) as input_handle:
        data = json.load(input_handle)

    basename = os.path.basename(path).rstrip(".json")

    motif_dict = defaultdict(list)

    data = [dat for datum in data for dat in datum if dat]

    for hit in data:

        if hit == None:
            continue
        blasthit = hit["blasthit"]
        epitope = hit["hitepitopedata"]
        motifs = hit["pdb_motif_hits"]

        for pdb_hit in motifs:
            for motif_type in pdb_hit:
                # print(motif_type["motif_seq"])
                if motif_type["TMalign_RMSD"] == " NaN":

                    continue
                # print(motif_type["TMalign_RMSD"])
                motif_dict[motif_type["motif_seq"]].append(
                    {
                        "epitope_id": epitope["epitope_id"],
                        "epitope_database": epitope["database"],
                        "source_antigen_accession": epitope["source_antigen_accession"],
                        "epitope_starting_position": epitope["starting_position"],
                        "epitope_ending_position": epitope["ending_position"],
                        "epitope_name": epitope["name"],
                        "epitope_organism_name": epitope["organism_name"],
                        "epitope_organism_id": epitope["organism_id"],
                        "epitope_internal_source_seq_acc": epitope["internal_source_seq_acc"],
                        "blast1_query_accession": blasthit["query_accession"],
                        "blast1_subject_accession": blasthit["subject_accession"],
                        "blast1_match_ranges": blasthit["match_ranges"],
                        "blast1_match_lengths": blasthit["match_lengths"],
                        "blast1_query_start": blasthit["query_start"],
                        "blast1_aln_query_seq": blasthit["aln_query_seq"],
                        "blast1_query_end": blasthit["query_end"],
                        "blast1_cigar": blasthit["cigar"],
                        "blast1_subject_start": blasthit["subject_start"],
                        "blast1_aln_subject_seq": blasthit["aln_subject_seq"],
                        "blast1_subject_end": blasthit["subject_end"],
                        "mmseqs2_query": motif_type["query"],
                        "mmseqs2_target": motif_type["target"],
                        "TMalign_png_file": os.path.basename(motif_type["TMalign_png_file"]),
                        "TMalign_RMSD": motif_type["TMalign_RMSD"],
                        "TMalign_TMscore": motif_type["TMalign_TMscore"],
                        "mmseqs2_motif_res_nums_target": motif_type["motif_res_nums_target"],
                        "mmseqs2_qcov": motif_type["qcov"],
                        "mmseqs2_pident": motif_type["pident"],
                        "mmseqs2_evalue": motif_type["evalue"],
                        "TMalign_PDB_file": motif_type["TMalign_PDB_file"],
                        "query_acc_motif": motif_type["query_acc_motif"],
                        "target_acc_motif": motif_type["target_acc_motif"],
                        "query_perc_acc": motif_type["query_perc_acc"],
                        "target_perc_acc": motif_type["target_perc_acc"],
                        "perc_acc_agree": motif_type["perc_acc_agree"],
                    }
                )

    for motif, data in motif_dict.items():
        motif_dict[motif] = sorted(data, key=lambda x: x["TMalign_RMSD"])

    motif_dict = {k: v for k, v in sorted(motif_dict.items(), key=lambda item: item[1][0]["TMalign_RMSD"])}

    with open(f"{config.OUTPUT_DIR}/{basename}_all_ranked.tsv", "w") as outhandle:
        outhandle.write("motif\t")
        w = csv.DictWriter(outhandle, list(motif_dict.items())[0][1][0].keys(), delimiter="\t")
        w.writeheader()
        for motif, data in motif_dict.items():

            for dat in data:
                outhandle.write(f"{motif}\t")
                w.writerow(dat)

    for key, instances in motif_dict.items():
        filtered_instacnes = []
        acc_visited = []
        for instance in instances:
            i_acc = int(instance["epitope_internal_source_seq_acc"])
            if i_acc in acc_visited:
                continue
            else:
                acc_visited.append(i_acc)
                filtered_instacnes.append(instance)

        motif_dict[key] = filtered_instacnes

    # all dicts are ordered dicts in 3.7 +, this breaks if using lower version of python.

    with open(f"{config.OUTPUT_DIR}/{basename}_best_per_source_seq.tsv", "w") as outhandle:
        outhandle.write("motif\t")
        w = csv.DictWriter(outhandle, list(motif_dict.items())[0][1][0].keys(), delimiter="\t")
        w.writeheader()
        for motif, data in motif_dict.items():

            for dat in data:
                outhandle.write(f"{motif}\t")
                w.writerow(dat)

    with open(f"{config.OUTPUT_DIR}/{basename}_best_per_source_seq.json", "w") as output_handle:
        json.dump(motif_dict, output_handle)
