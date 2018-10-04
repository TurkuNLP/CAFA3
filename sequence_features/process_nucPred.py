#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import os
import argparse
import glob
import gzip
import time
import requests
import webbrowser
from requests_toolbelt.multipart.encoder import MultipartEncoder
from Bio import SeqIO


def prepare_sequence(in_folder):
    all_seq = []
    non_anal = []
    for filename in glob.glob(in_folder + '*'):
        for record in SeqIO.parse(filename, 'fasta'):
            all_seq.append(record)
    print('{} out of {} will be analyzed'.format(str(len(all_seq)), str(len(all_seq) + len(non_anal))))
    return all_seq


def run_nucPred(all_seq, out_folder):
    """
    maximum number of sequences submitted = 1000
    """
    max_seq = 1000
    i = 0
    text = []
    while i < len(all_seq):
        print('processing nucPred sequence #', i, i + max_seq)
        try:
            SeqIO.write(all_seq[i: i+n], out_folder + 'temp.fasta', 'fasta')
        except:
            SeqIO.write(all_seq[i: len(all_seq)], out_folder + 'temp.fasta', 'fasta')
        i += max_seq
        text = nucPred_submission(text, out_folder)        
    return text

    
def nucPred_submission(text, out_folder):
    with open(out_folder + 'temp.fasta', 'rb') as f:
        multipart_data = MultipartEncoder(
            fields={'upload': (out_folder + 'temp.fasta', f, 'fasta'),
                    '.submit': '.submit'})        
        response = requests.post(
            url=u'https://nucpred.bioinfo.se/cgi-bin/batch.cgi',
            data=multipart_data,
            headers={u'Content-Type': multipart_data.content_type})
        for item in response.content.decode().split('NucPred-score\n')[1].split('</pre>')[0].split("\n"):            
            if item:
                text.append('\t'.join(item.split()))
    return text


def write_nucPred_result(out_folder, text):                
    if os.path.exists(out_folder + 'CAFA3_nucPred.tsv.gz'):
        with gzip.open(out_folder + 'CAFA3_nucPred.tsv.gz', 'a') as f:
            f.write("\n".join(text) + "\n")
    else:
        with gzip.open(out_folder + 'CAFA3_nucPred.tsv.gz', 'wt') as f:
            f.write("Sequence-ID\tNucPred-score\n" + "\n".join(text) + "\n")


def run_nucPred(all_seq, out_folder):
    """
    maximum number of sequences submitted = 1000
    """
    max_seq = 1000
    i = 0
    text = []
    while i < len(all_seq):
        print('processing nucPred sequence #', i, i + max_seq)
        try:
            SeqIO.write(all_seq[i: i+n], out_folder + 'temp.fasta', 'fasta')
        except:
            SeqIO.write(all_seq[i: len(all_seq)], out_folder + 'temp.fasta', 'fasta')
        i += max_seq
        text = nucPred_submission(text, out_folder)
    return text


def check_header_protein_count(all_seq, out_folder):
    protein_res = set([])
    for filename in glob.glob(out_folder + '*.gz'):
        with gzip.open(filename, 'rt') as f:
            csv_reader = csv.reader(f, delimiter = '\t')
            headers = next(csv_reader)
            for row in csv_reader:
                if row != '\n':
                    protein_res.add(row[0])
                    assert len(row) == len(headers) #rows and header is not the same number of columns
    assert protein_res == set([seq.id for seq in all_seq]) # the number of protein results are not equal to the number of queried proteins

    
def main(in_folder, out_folder):
    all_seq = prepare_sequence(in_folder)
    text = run_nucPred(all_seq, out_folder)    
    write_nucPred_result(out_folder, text)
    check_header_protein_count(all_seq, out_folder)

    
def argument_parser():
    parser = argparse.ArgumentParser(description="targetP prediction")
    parser.add_argument("-i", "--in_folder", type=str, help="input sequence folder")
    parser.add_argument("-o", "--out_folder", type=str, help="output folder")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    main(args.in_folder, args.out_folder)
