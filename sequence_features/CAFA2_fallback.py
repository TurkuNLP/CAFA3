#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import Counter
import gzip
import glob
import cPickle as pickle
import operator
import argparse

# 1. take go_id from uniprot_go order by its annotation frequency
# 2. take blast score of each cafa protein with maximum score


def load_go(go_file):
    with gzip.open(go_file, "rb") as f:
        uniprot_go = [line.strip("\n").split("\t")[0] for line in f]
    return set(uniprot_go)


def load_predict(in_file):
    with gzip.open(in_file, "rb") as f:
        all_predict = [line for i, line in enumerate(f) if i != 1]
    return all_predict


def blast_hit(in_folder, go_file, out_folder, cond): 
    uniprot_go = load_go(go_file)
    print len(uniprot_go), list(uniprot_go)[0:10]
    for filename in glob.glob(in_folder + cond):        
        if filename.split("/")[-1] not in [""]:            
            with gzip.open(filename, "rb") as f:
                all_line = [line for line in f if line.split("\t")[4] in uniprot_go]
            with gzip.open(out_folder + filename.split("/")[-1], "wb") as f:
                f.write("".join(all_line))

                
def load_go_protein(go_propagate):
    prot_dict = dict()
    with gzip.open(go_propagate, "rb") as f:
        uniprot_go = [line.strip("\n").split("\t")[0:2] for line in f]
    for item in uniprot_go:
        prot_dict.setdefault(item[0], []).append(item[1])
    return prot_dict


def load_prot_function_mappings(swissprot_evidence_filename, topn):
    prot_function_mappings = {}
    topn_set = None
    if topn > 0:
        func_counts = {}
        with gzip.open(swissprot_evidence_filename, 'rt') as file:
            for line in file:
                line = line.strip()
                line_parts = line.split('\t')
                if line_parts[1] in func_counts:
                    func_counts[line_parts[1]] += 1
                else:
                    func_counts[line_parts[1]] = 1
        topn_set = set()
        for i, (func, count) in enumerate(sorted(func_counts.items(), key=lambda x: x[1], reverse=True)):
            if i < topn:
                topn_set.add(func)
            else:
                break
    return topn_set


def transfer_go(in_folder, topN, prot_dict, top_5000, out_file, protein_list, cond):
    go_match = dict()
    for filename in glob.glob(in_folder + cond):        
        blast_match = dict()
        with gzip.open(filename, "rb") as f:
            for line in f:
                lines = line.split("\t")
                if lines[0] != lines[4] and lines[0] in protein_list:
                    blast_match.setdefault(lines[0], []).append(lines[4])
        for protein, v in blast_match.iteritems():            
            k = 0
            for i, item in enumerate(v):
                try:
                    if k < topN:
                        for goid in prot_dict[item]:
                            if goid not in top_5000:
                                go_match.setdefault(protein, []).append(goid)
                    k += 1
                except:
                    pass
    new_dict = dict()
    for k, v in go_match.iteritems():
        new_dict.setdefault(k, Counter(v))
    text = []
    for k, v in new_dict.iteritems():
        for goid, count_p in v.iteritems():
            if count_p > 2:
                text.append([k, goid, str(count_p)])
    new_list = sorted(text, key=operator.itemgetter(0, 2), reverse=True)
    new_text = []
    for item in new_list:
        new_text.append("\t".join(item))
    with gzip.open(in_folder + "trial/" + out_file, "wb") as f:
        f.write("\n".join(new_text))


def transfer_go_cafa(in_folder, topN, prot_dict, top_5000, out_file, map_dict, sub_folder):
    print protein_list
    go_match = dict()    
    for filename in glob.glob(in_folder + sub_folder + '*'):        
        blast_match = dict()
        with gzip.open(filename, "rb") as f:
            for line in f:
                lines = line.split("\t")
                
                if lines[0] != lines[4]:
                    blast_match.setdefault(lines[0], []).append(lines[4])        
        for protein, v in blast_match.iteritems():
            print protein, v
            k = 0
            for i, item in enumerate(v):
                try:                    
                    if k < topN:                        
                        for goid in prot_dict[item]:
                            if goid in top_5000:
                                print goid
                                go_match.setdefault(protein, []).append(goid)
                    k += 1
                except:
                    pass
    new_dict = dict()
    for k, v in go_match.iteritems():        
        new_dict.setdefault(k, Counter(v))
    text = []
    for k, v in new_dict.iteritems():
        for goid, count_p in v.iteritems():
            if count_p > 2:
                if goid in map_dict.keys():
                    text.append([k, goid, str(count_p)])
    new_list = sorted(text, key=operator.itemgetter(0, 2), reverse=True)
    new_text = []
    for item in new_list:        
        # print item, new_dict[item[0]]
        new_text.append("\t".join(item))
    with gzip.open(in_folder + "trial/" + out_file, "wb") as f:
        f.write("\n".join(new_text))
        

def argument_parser():
    parser = argparse.ArgumentParser(description="process fallback system following CAFA2 submission")
    parser.add_argument("-g", "--go_evidence", type=str, default="/home/kahaka/CAFA3/data_save/Swissprot_evidence.tsv.gz", help="full path to go evidence file from Swissprot")
    parser.add_argument("-i", "--go_propagate", type=str, default="/home/kahaka/CAFA3/data_save/Swissprot_propagated.tsv.gz", help="full path to go propagated file got from Kai")
    parser.add_argument("-o", "--go_pair", type=str, default="/home/kahaka/CAFA3/data_save/parent_term2term.txt.gz", help="full path to go-pairs file got from Kai")    
    parser.add_argument("-d", "--data_folder", type=str, default="/home/sukaew/CAFA_QA/fallback/", help="folder for fallback file to be stored")
    parser.add_argument("-p", "--prot_file", type=str, default="/home/kahaka/projects/cafa_pi/data/develtest.txt.gz", help="input file from Jari")
    parser.add_argument("-f", "--out_file", type=str, default="Jari_develtest.tsv.gz", help="output filename")
    parser.add_argument("-b", "--blast_features", type=str, default="/home/sukaew/CAFA_QA/features/CAFA2/CAFA3_features/", help="Blast feature folder")
    parser.add_argument("-s", "--sub_folder", type=str, default="targetFiles/", help="sub folder where ")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    # with gzip.open("/home/sukaew/CAFA3/fallback/devel.txt.gz", "rb") as f:
    #     protein_list = [line.strip("\n") for i, line in enumerate(f)]
    # protein_set = set(protein_list)
    # transfer_go("/home/sukaew/CAFA3/fallback/", 5, prot_dict, set([]), "Jari_devel.tsv.gz", protein_set, "Swissprot*.features_tsv.gz")

    # with gzip.open("/home/sukaew/CAFA3/fallback/test.txt.gz", "rb") as f:
    #     protein_list = [line.strip("\n") for i, line in enumerate(f)]
    # protein_set = set(protein_list)
    # transfer_go("/home/sukaew/CAFA3/fallback/", 5, prot_dict, set([]), "Blast_test.tsv.gz", protein_set, "Swissprot*.features_tsv.gz")
    
    # blast_hit("/home/sukaew/CAFA_PI/features/CAFA2/training_features/", args.go_evidence, "/home/sukaew/CAFA_PI/fallback/training/", "*features*")
    # blast_hit("/home/sukaew/CAFA_PI/features/CAFA2/CAFA3_features/", args.go_evidence, "/home/sukaew/CAFA_QA/fallback/targetFiles/", "*features*")

    # blast_hit(args.blast_features, args.go_evidence, args.data_folder + args.sub_folder, "*features*")

    top_5000 = load_prot_function_mappings(args.go_propagate, 5000)
    
    prot_dict = load_go_protein(args.go_propagate)

    map_dict = {}
    with gzip.open(args.go_pair, 'rb') as f:
        for line in f:
            lines = line.strip('\n').split('\t') 
            if lines[0] in top_5000:
                map_dict.setdefault(lines[1], lines[0])
                map_dict.setdefault(lines[0], lines[0])
    
    with gzip.open(args.prot_file, 'rb') as f:        
        protein_list = set([line.strip('\n') for line in f])

    # "Jari_devel.tsv.gz", "Jari_test.tsv.gz"
    transfer_go_cafa(args.data_folder, 5, prot_dict, top_5000, args.out_file, map_dict, args.sub_folder) # "target*.features_tsv.gz"
