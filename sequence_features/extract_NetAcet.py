#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import argparse



def extract_data(in_file, out_file, parent_file, org_file):
    data_folder = "/media/suwisa/CAFA3/"
    
    with gzip.open(data_folder + parent_file, "rb") as f:
        all_org = [line.split("\t")[0] for i, line in enumerate(f) if "\t2759\t" in line and i != 0]
        
    with gzip.open(data_folder + org_file, "rb") as f:
        all_protein = [line.split("\t")[0] for i, line in enumerate(f) if line.strip("\n").split("\t")[1] in all_org and i != 0]

    print len(all_org)
    print len(all_protein)
    
    with open(data_folder + "map_cafa_organism.txt", "rb") as f:
        name_dict = {line.split("\t")[0]: line.split("\t")[1] for i, line in enumerate(f) if i != 0} 
    loc_dict = {'1': 0, '2': 4, '3': 8}
    netacet = []
    pred_dict = dict()
    with gzip.open(data_folder + in_file, "rb") as f:
        value_list = 4*['0.0']
        empty_list = 4*['-']
        for line in f:
            if "no Ala, Gly, Ser or Thr at positions 1-3" in line:
                protein_id = line.split("\t")[0]
            else:
                protein_id = line.split()[0]
            pred_dict.setdefault(protein_id, []).append(line)

    for k, v in pred_dict.iteritems():    
        text_list = value_list*3
        for item in v:
            if "no Ala, Gly, Ser or Thr at positions 1-3" not in item:
                item_l = item.split()
                i = loc_dict[item_l[1]]
                text_list[i:i+4] = item_l[4:8]
        else:
            try:
                if name_dict[k] in all_protein:
                    netacet.append("\t".join([name_dict[k]] + text_list))
                # else:
                #     netacet.append("\t".join([name_dict[k]] + empty_list*3))
            except:
                if k in all_protein:
                    netacet.append("\t".join([k] + text_list))
                # else:
                #     netacet.append("\t".join([k] + empty_list*3))
        
    print len(netacet)
    with gzip.open(data_folder + out_file, "wb") as f:
        f.write("\n".join(netacet) + "\n")


        
extract_data("Swissprot_sequence_NetAcet.tsv.gz", "Swissprot_sequence_NetAcet_extract.tsv.gz", "map_training_organism_parent.tsv.gz", "map_training_organism.tsv.gz")


extract_data("target_all_NetAcet.tsv.gz", "target_all_NetAcet_extract.tsv.gz", "map_cafa3target_organism_parent.tsv.gz", "map_cafa2target_organism.tsv.gz")

