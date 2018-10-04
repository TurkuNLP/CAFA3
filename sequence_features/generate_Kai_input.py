#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import glob
import argparse


def prepare_sequence(in_folder, out_folder):
    all_ids = []
    all_seq = []
    for filename in glob.glob(in_folder + '*'):
        with open(filename, 'r') as f:
            all_lines = f.read()
        for seq in all_lines.split('>'):
            seq_split = seq.split("\n", 1)
            seq_len = "".join(seq_split[1:])
            if seq_len != '':
                all_ids.append('>' + filename.split('/')[-1].split('.')[0] + '\n')
                all_seq.append('>' + filename.split('/')[-1].split('.')[0] + '\n' + seq_len.replace('\n', '') + '\n')
    with gzip.open(out_folder + 'target.all.ids.gz', 'wb') as f:
        f.write(''.join(all_ids))
    with gzip.open(out_folder + 'target.all.fasta.gz', 'wb') as f:
        f.write(''.join(all_seq))
    

def argument_parser():
    parser = argparse.ArgumentParser(description="process fallback system")
    parser.add_argument("-i", "--in_folder", type=str, default="/home/sukaew/CAFA_QA/targetFiles/sequences/", help="folder for analyzed files")
    parser.add_argument("-o", "--out_folder", type=str, default="/home/sukaew/CAFA_QA/targetFiles/", help="folder for analyzed files")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    prepare_sequence(args.in_folder, args.out_folder)
