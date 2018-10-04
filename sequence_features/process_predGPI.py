#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import csv
import os
import glob
import gzip
import time
import argparse
import requests
import webbrowser
# from requests_toolbelt.multipart.encoder import MultipartEncoder
from requests_toolbelt import MultipartEncoder
from Bio import SeqIO


def check_sequence(in_folder):
    """
    minimum length of sequence = 40 residues 
    """
    all_seq = []
    non_anal = []
    for filename in glob.glob(in_folder + '*'):
        for record in SeqIO.parse(filename, 'fasta'):
            if len(record.seq) > 40:
                all_seq.append(record)
            else:
                non_anal.append(record.id)
                print(record.id, 'sequence length is less than 40')
    print('{} out of {} will be analyzed'.format(str(len(all_seq)), str(len(all_seq) + len(non_anal))))    
    return all_seq


def run_predGPI(all_seq, out_folder):
    """
    maximum number of sequences submitted = 500
    """
    max_seq = 500
    i = 0
    text = []
    while i < len(all_seq):
        print('processing predGPI sequence #', i, i + max_seq)
        try:
            SeqIO.write(all_seq[i: i+n], out_folder + 'temp.fasta', 'fasta')
        except:
            SeqIO.write(all_seq[i: len(all_seq)], out_folder + 'temp.fasta', 'fasta')
        i += max_seq
        text = predGPI_submission(text, out_folder)
    return text


def predGPI_submission(text, out_folder):
    with open(out_folder + 'temp.fasta', 'rb') as f:
        print('running GPI prediction') 
        multipart_data = MultipartEncoder(
            fields={
                'upfile': f,
                'tipo_hmm': '0',
                'Display': 'Display'})        
        response = requests.post(
            url='http://gpcr2.biocomp.unibo.it/cgi-bin/predictors/gpi/download_fasta_1.4.cgi',
            data=multipart_data,
            headers={'Content-Type': multipart_data.content_type})        
        for item in response.content.decode().split(">")[1:]:
            protein_id, FPrate, OMEGA = item.split("\n")[0].split(' | ')
            try:
                text.append("\t".join([protein_id.split(' ', 1)[0], FPrate.split(':')[1], OMEGA.split('-')[1]]))
            except:
                text.append("\t".join([protein_id, FPrate.split(':')[0], OMEGA.split('-')[1]]))
        time.sleep(10)
    return text


def write_predGPI_result(out_folder, text):
    if os.path.exists(out_folder + 'new_CAFA3_predGPI.tsv.gz'):
        with gzip.open(out_folder + 'new_CAFA3_predGPI.tsv.gz', 'a') as f:
            f.write("\n".join(text) + "\n")
    else:
        with gzip.open(out_folder + 'new_CAFA3_predGPI.tsv.gz', 'wt') as f:
            f.write("protein_id\tFPrate\tOMEGA\n" + "\n".join(text) + "\n")


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
    all_seq = check_sequence(in_folder)
    text = run_predGPI(all_seq, out_folder)    
    write_predGPI_result(out_folder, text)
    check_header_protein_count(all_seq, out_folder)
                        
                        
def argument_parser():
    parser = argparse.ArgumentParser(description="predicting GPI anchor")
    parser.add_argument("-i", "--in_folder", type=str, help="full path to the folder")
    parser.add_argument("-o", "--out_folder", type=str, help="full path to the folder")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    main(args.in_folder, args.out_folder)
