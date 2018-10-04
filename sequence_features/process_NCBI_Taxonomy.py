#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from common_processing import *
import tarfile
import sys
import glob

 
def untar(ftp_link, out_folder):
    tar = tarfile.open(out_folder + ftp_link.split("/")[-1])
    tar.extractall(path=out_folder)
    tar.close()
    
    
def process_nodes_dmp(out_folder):
    """
    extract data from file:nodes.dmp
    
    create 2 map_tables:
    map_organism2organism
    map_organism2rank
    """
    map_organism2organism = ""
    map_organism2rank = ""
    parent_tax_dict = dict()
    tax_tree_dict = dict()
    with open(out_folder + 'nodes.dmp', 'rb') as f:
        for line in f:
            tax_id, parent_tax_id, rank, embl_code, division_id, inherited_div_flag, genetic_code_id, inherited_gc_flag, mitochondrial_genetic_code_id, inherited_mgc_flag, genbank_hidden_flag, hidden_subtree_root_flag, comments = line.split("\t|\t")
            map_organism2rank += str(tax_id) + "\t" + rank + "\n"
            parent_tax_dict.setdefault(tax_id, parent_tax_id)
    for tax_id, parent_tax_id in parent_tax_dict.iteritems():
        tax_tree_dict.setdefault(tax_id, []).append(parent_tax_id)
        while parent_tax_dict[tax_tree_dict[tax_id][-1]] != tax_tree_dict[tax_id][-1]:
            tax_tree_dict[tax_id].append(parent_tax_dict[tax_tree_dict[tax_id][-1]])
    for tax_id, parent_tax_ids in tax_tree_dict.iteritems():
        map_organism2organism += '{}\t{}\t{}\n'.format(tax_id, tax_id, 0)
        for level, parent_tax_id in enumerate(parent_tax_ids):
            map_organism2organism += '{}\t{}\t{}\n'.format(tax_id, parent_tax_id, level+1)
    with open(out_folder + "map_organism2organism.tsv", "wb") as f:
        f.write(map_organism2organism)
    with open(out_folder + "map_organism2rank.tsv", "wb") as f:
        f.write(map_organism2rank)

        
def process_names_dmp(out_folder):
    """ 
    extract data from file:names.dmp
    map_symbol2organism
    name type included: scientific name, synonym, acronym, anamorph, misspelling, misnomer, common name, 
    """
    map_symbol2organism = ''
    non_unique_name = set()
    with open(out_folder + "names.dmp", "rb") as f:
        for line in f:
            tax_id, name_txt, unique_name, name_class = line.split("\t|\t")
            map_symbol2organism += "{}\t{}\t{}".format(tax_id, name_txt, name_class.split("|")[0].replace("\t", "\n")) 
    with open(out_folder + "map_symbol2organism.tsv", "wb") as f:
        f.write(map_symbol2organism)


def argument_parser():

    parser = argparse.ArgumentParser(description="download the Taxonomy PubMed from ftp")
    parser.add_argument("-f", "--ftp_link", type=str, help="ftp url link to the file")
    parser.add_argument("-o", "--out_folder", type=str, help="target folder of downloaded file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    print "processing Taxonomy data"
    ftp_download(args.ftp_link, args.out_folder)
    untar(args.ftp_link, args.out_folder)
    process_nodes_dmp(args.out_folder)
    process_names_dmp(args.out_folder)
