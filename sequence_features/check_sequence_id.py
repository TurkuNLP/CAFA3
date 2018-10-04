#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import argparse
from Bio import SeqIO
import gzip
import csv


def check_modify_fasta_file(seq_count, protein_pair, all_seq, in_file):
    for record in SeqIO.parse(in_file, 'fasta'):
        ori_record_id = record.id
        record.id = 'protein_' + str(seq_count)
        protein_pair.append([ori_record_id, record.id])
        seq_count += 1
        all_seq.append(record)
    return seq_count, protein_pair, all_seq


def write_protein_seq(all_seq, out_folder):
    max_seq = 1000
    new_seq = []
    j = 0
    for record in all_seq:        
        if len(new_seq) == max_seq:
            SeqIO.write(new_seq, out_folder + 'protein_seq_{}.fasta'.format(str(j)), 'fasta')
            j += 1
            new_seq = [record]
        new_seq.append(record)
    if new_seq != []:
        SeqIO.write(new_seq, out_folder + 'protein_seq_{}.fasta'.format(str(j)), 'fasta')
        
    
def write_protein_pair(protein_pair, out_folder):
    with gzip.open(out_folder + 'protein_id_pair.txt.gz', 'wt') as f:
        csv_writer = csv.writer(f, delimiter='\t')
        csv_writer.writerows(protein_pair)

    
def main(seq_folder, out_folder, main_folder):
    seq_count = 0
    all_seq = []
    protein_pair = [['ori_protein_id', 'new_protein_id']]
    for filename in glob.glob(seq_folder + '*'):
        seq_count, protein_pair, all_seq = check_modify_fasta_file(seq_count, protein_pair, all_seq, filename)        
    write_protein_seq(all_seq, out_folder)
    write_protein_pair(protein_pair, main_folder)

    
def argument_parser():
    parser = argparse.ArgumentParser(description="change protein id so there is consistencies in protein names throughout all analysis")
    parser.add_argument("-s", "--seq_folder", type=str, help="folder where sequences reside")
    parser.add_argument("-o", "--out_folder", type=str, help="folder where new combined sequences go")
    parser.add_argument("-m", "--main_folder", type=str, help="folder where file list goes")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    main(args.seq_folder, args.out_folder, args.main_folder)
