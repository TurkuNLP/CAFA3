#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
from Bio import SeqIO
import argparse
import re
import gzip                
import os


def process_nodes_dmp(tax_folder):
    """
    extract data from file:nodes.dmp    
    create 2 map_tables:
    map_org2org
    map_org2rank
    """
    parent_tax_dict = dict()
    tax_tree_dict = dict()
    hierarchy_dict = dict()
    with gzip.open(tax_folder + 'nodes.dmp.gz', 'rt') as f:
        for line in f:
            tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div_flag, genetic_code_id, inherited_gc_flag, mitochondrial_genetic_code_id, inherited_mgc_flag, genbank_hidden_flag, hidden_subtree_root_flag, comments = line.split("\t|\t")
            parent_tax_dict.setdefault(tax_id, parent_tax_id)    
    for tax_id, parent_tax_id in parent_tax_dict.items():
        tax_tree_dict.setdefault(tax_id, []).append(parent_tax_id)
        while parent_tax_dict[tax_tree_dict[tax_id][-1]] != tax_tree_dict[tax_id][-1]:
            tax_tree_dict[tax_id].append(parent_tax_dict[tax_tree_dict[tax_id][-1]])
    for tax_id, parent_tax_ids in tax_tree_dict.items():
        for level, parent_tax_id in enumerate(parent_tax_ids):
            hierarchy_dict.setdefault(int(tax_id), []).append(int(parent_tax_id))
    return hierarchy_dict


def process_names_dmp(tax_folder):
    """ 
    extract data from file:names.dmp
    name type: scientific name, synonym, acronym, anamorph, misspelling, misnomer, common name 
    """
    map_symb2org = dict()
    with gzip.open(tax_folder + "names.dmp.gz", "rt") as f:
        for line in f:
            tax_id, name_txt, unique_name, name_class = line.split("\t|\t")
            if name_class == 'scientific name\t|\n':
                map_symb2org.setdefault(int(tax_id), name_txt)
    return map_symb2org


def fasta_organism(in_folder):
    org_dict = dict()
    for in_file in glob.glob(in_folder + '*'):
        for record in SeqIO.parse(in_file, 'fasta'):
            desc = record.description            
            fields = re.finditer(r'\s(\w+)=', desc)
            for item in fields:
                if item.group() == ' OS=':
                    offB = item.span()[1]
                    next_item = next(fields)
                    offE = next_item.span()[0]
                    tax_name = desc[offB:offE]                          
                    if next_item.group() == ' OX=':
                        offB = next_item.span()[1]
                        next_item = next(fields)
                        offE = next_item.span()[0]
                        tax_id = desc[offB:offE]                                
                elif item.group() == ' OX=':
                    offB = item.span()[1]
                    next_item = next(fields)
                    offE = next_item.span()[0]
                    tax_name = desc[offB:offE]
            org_dict.setdefault(record.id, [int(tax_id), tax_name])
    return org_dict                    


def map_names(hierarchy_dict, map_symb2org, org_dict, out_folder):
    org_text = ''
    found_tax = dict()
    for protein_id, tax in org_dict.items():
        tax_id, tax_name = tax
        try:
            org_list = found_dict[tax_id]
        except:
            org_list = [map_symb2org[parent] for parent in [tax_id] + hierarchy_dict[tax_id]]
            found_tax.setdefault(tax_id, org_list)
        org_text += protein_id + '\t' + ','.join(org_list[::-1]) + '\n'
    with gzip.open(out_folder + 'map_target_taxonomy.tsv.gz', 'wt') as f:
        f.write(org_text)

        
def main(in_folder, tax_folder, out_folder):
    org_dict = fasta_organism(in_folder)
    hierarchy_dict = process_nodes_dmp(tax_folder)
    map_symb2org = process_names_dmp(tax_folder)
    map_names(hierarchy_dict, map_symb2org, org_dict, out_folder)
        
    
def argument_parser():
    parser= argparse.ArgumentParser(description='obtain taxonomy lineage from sequence fasta file')
    parser.add_argument('-i', '--in_folder', type=str, help='processed fasta file folder')
    parser.add_argument('-o', '--out_folder', type=str, help='feature file folder')
    parser.add_argument('-t', '--tax_folder', type=str, help='taxonomy folder')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args =argument_parser()
    main(args.in_folder, args.tax_folder, args.out_folder)
