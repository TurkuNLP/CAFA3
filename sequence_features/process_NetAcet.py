#!/usr/bin/env python
# -*- coding: utf-8 -*-

import csv
import glob
import urllib
import gzip
import time
import argparse
import requests
import webbrowser
from requests_toolbelt.multipart.encoder import MultipartEncoder
from Bio import SeqIO


def prepare_sequence(in_folder):
    """
    minimum length of sequence = 40 residues 
    maximum length of sequence = 4,000 residues
    """
    all_seq = []
    non_anal = []
    for filename in glob.glob(in_folder + '*'):
        for record in SeqIO.parse(filename, 'fasta'):
            if len(record.seq) > 40 and len(record.seq) <= 4000:
                all_seq.append(record)
            else:
                non_anal.append(record.id)
                print(record.id, 'sequence length is not in [40, 4000]')   
    print('NetAcet: {} out of {} will be analyzed'.format(str(len(all_seq)), str(len(all_seq) + len(non_anal))))
    return all_seq


def NetAcet_submission(result_link, out_folder):    
    with open(out_folder + 'temp.fasta', 'rb') as f:
        multipart_data = MultipartEncoder(
            fields={'SEQSUB': (out_folder + 'temp.fasta', f, 'fasta'), 'configfile': '/usr/opt/www/pub/CBS/services/NetAcet-1.0/NetAcet.cf', 'inns': 'button'})
        response = requests.post(url='http://www.cbs.dtu.dk/cgi-bin/webface2.fcgi', data=multipart_data, headers={'Content-Type': multipart_data.content_type})
        status = response.content.decode().split('<span name="status">')[1].split('</span></H1>')[0]
        link = response.content.decode().split('a href="')[1].split('">click here')[0]
        result_link.append(link)            
    return result_link


def run_NetAcet(all_seq, out_folder):
    """
    maximum number of sequences = 2000
    maximum number of residues of all sequences = 200000
    """
    max_seq = 2000
    max_residues = 200000
    result_link = []    
    all_residues = 0
    new_seq = []
    for record in all_seq:        
        all_residues += len(record.seq)        
        if all_residues > max_residues or len(new_seq) == max_seq:
            SeqIO.write(new_seq, out_folder + 'temp.fasta', 'fasta')
            result_link = NetAcet_submission(result_link, out_folder)
            all_residues = len(record_seq)
            new_seq = [record]
        new_seq.append(record)        
    if new_seq != []:
        SeqIO.write(new_seq, out_folder + 'temp.fasta', 'fasta')
        result_link = NetAcet_submission(result_link, out_folder)
    with gzip.open(out_folder + "temp.lnk.gz", 'wt') as f:
        f.write("\n".join(result_link) + "\n")


def read_links(out_folder):
    all_text = []
    with gzip.open(out_folder + "temp.lnk.gz", 'rt') as f:
        result_link = f.readlines()
    for link in result_link:        
        status = urllib.request.urlopen(link).read().decode().split('<span name="status">')[1].split('</span></H1>')[0]
        while status == 'queued':
            time.sleep(20)
            try:
                status = urllib.request.urlopen(link).read().decode().split('<span name="status">')[1].split('</span></H1>')[0]            
                print('{}: {}'.format(link, status))
            except:
                status = 'active'
                time.sleep(10)
        result = urllib.request.urlopen(link).read().decode()
        all_text.append(result.split("-"*68)[1].split("<font face=")[0].strip("\n"))                
    return '\n'.join(all_text)


def write_NetAcet_result(text, out_folder):
    """ 
    """
    loc_dict = {'1': 0, '2': 4, '3': 8}
    netacet = ['\t'.join(['protein_id', 'val_1', 'val_2', 'val_3', 'val_4', 'val_5', 'val_6', 'val_7', 'val_8', 'val_9', 'val_10', 'val_11', 'val_12'])]
    pred_dict = dict()     
    value_list = 4*['0.0']
    empty_list = 4*['-']
    for line in text.split('\n'):        
        text_list = value_list*3
        if " no Ala, Gly, Ser or Thr at positions 1-3" in line:
            protein_id = line.split("\t")[0]
        else:                
            protein_id = line.split()[0]
        pred_dict.setdefault(protein_id, []).append(line)
    for k, v in pred_dict.items():    
        text_list = value_list*3
        for item in v:
            if "no Ala, Gly, Ser or Thr at positions 1-3" not in item:
                item_l = item.split()
                i = loc_dict[item_l[1]]
                text_list[i:i+4] = item_l[4:8]
            netacet.append("\t".join([k] + text_list))
    with gzip.open(out_folder + 'target_sequence_NetAcet_extract.tsv.gz', "wt") as f:
        f.write("\n".join(netacet) + "\n")

        
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
    run_NetAcet(all_seq, out_folder)
    text = read_links(out_folder)
    write_NetAcet_result(text, out_folder)
    check_header_protein_count(all_seq, out_folder)
    

def argument_parser():
    parser = argparse.ArgumentParser(description="targetP prediction")
    parser.add_argument("-i", "--in_folder", type=str, help="inpit sequence folder")
    parser.add_argument("-o", "--out_folder", type=str, help="output folder")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = argument_parser()
    main(args.in_folder, args.out_folder) 
