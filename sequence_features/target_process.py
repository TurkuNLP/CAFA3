#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
from os.path import basename
import os
import shutil
import sys
import subprocess
import zipfile
import tarfile
from Bio import SeqIO
import argparse

import check_sequence_id as cs
import process_predGPI as gpi
import process_NetAcet as na
import process_nucPred as np
import process_Taxonomy as tx

def runCommand(cmd_list):
    try:
        subprocess.check_call(cmd_list, shell=True)
    except subprocess.CalledProcessError:
        print(cmd_list)
        pass # handle errors in the called executable
    except OSError:
        print(cmd_list)
        pass # executable not found

    
def create_dir(input_dir):
    if not os.path.exists(input_dir):
        os.mkdir(input_dir)
    

def check_input_sequence(fname, target_dir):
    if not os.path.isdir(fname):        
        if (fname.endswith("tar.gz")):
            tar = tarfile.open(fname, "r:gz")
            tar.extractall(path=target_dir)
            tar.close()
            seq_dir = fname.replace(".tar.gz", "/")
        elif (fname.endswith("tar")):
            tar = tarfile.open(fname, "r:")
            tar.extractall(path=target_dir_)
            tar.close()            
            seq_dir = fname.replace(".tar", "/")
        elif (fname.endswith("zip")):
            zip_ref = zipfile.ZipFile(fname, 'r')
            ret = zip_ref.testzip()
            if ret is not None:
                print("First bad file in zip: %s" % ret)
                sys.exit(1)
            else:
                zip_ref.extractall(target_dir)
                zip_ref.close()
                seq_dir = fname.replace(".zip", "/")
        else:
            shutil.copyfile(fname, target_dir)
    else:
        seq_dir = fname
    print(seq_dir)
    return seq_dir


def check_fasta(seq_file):
    with open(seq_file, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        if not any(fasta):
            print('remove non fasta file:', seq_file)
            os.remove(seq_file)  # False when `fasta` is empty, i.e. wasn't a FASTA file

            
def copy_feature_files(out_folder):
    feat_dict = {'CAFA2/CAFA3/': 'blastp62_features/',
                     'blast_result_features/': 'blastp45_features/',
                     'deltaBlast_result_features/': 'deltaBlast_features/',
                     'interproscan_result_features/': 'interproscan_features/',
                     'Taxonomy/': 'Taxonomy/',
                     'predGPI/': 'predGPI/',
                     'NetAcet/': 'NetAcet/',
                     'nucPred/': 'nucPred/'}    
    for feat in feat_dict.keys():
        if feat == 'CAFA2/CAFA3/':
            create_dir(out_folder + 'features/' + 'CAFA2/')            
        create_dir(out_folder + 'features/' + feat)
        print(out_folder + 'targetFiles/' + feat_dict[feat])
        for filename in glob.glob(out_folder + 'targetFiles/' + feat_dict[feat] + '*.gz'):
            # print(filename)
            shutil.copy2(filename, out_folder +'features/' + feat)
        

def main(out_folder, ori_seq):    
    CAFA_scripts = os.getcwd()
    for subdir in glob.glob(out_folder + '*/'):
        shutil.rmtree(subdir)

    input_dir = out_folder + 'targetFiles/'
    create_dir(input_dir)
    tax_folder = os.path.join((os.getcwd() + '/taxonomy/'))
    
    input_seq = input_dir + 'sequences/'
    create_dir(input_seq)
    seq_dir = check_input_sequence(out_folder + ori_seq, out_folder)
    for seq_file in glob.glob(seq_dir + '*'):
        check_fasta(seq_file)
    cs.main(seq_dir, input_seq, out_folder)
    with open(input_dir + 'fasta_filelist', 'wt') as f:
        f.write('\n'.join(glob.glob(input_seq + '*')) + '\n')    

    features = ['NetAcet/', 'nucPred/', 'Taxonomy/', 'blastp62/', 'blastp62_features/', 'blastp45/', 'blastp45_features/', 'deltaBlast/', 'deltaBlast_features/', 'interproscan/', 'interproscan_features/', 'Taxonomy', 'predGPI/']
    for feature in features:
        create_dir(input_dir + feature)
        if feature == 'predGPI/':
            gpi.main(input_seq, input_dir + feature)
            os.remove(input_dir + feature + 'temp.fasta')
        elif feature == 'NetAcet/':
            na.main(input_seq, input_dir + feature)
            os.remove(input_dir + feature + 'temp.fasta')
            os.remove(input_dir + feature + 'temp.lnk.gz')
        elif feature == 'nucPred/':
            np.main(input_seq, input_dir + feature)
            os.remove(input_dir + feature + 'temp.fasta')
        elif feature == 'Taxonomy/':
            tx.main(input_seq, tax_folder, input_dir + feature)

    subprocess.call('bash ' + CAFA_scripts + '/target_process.sh', shell=True)
    subprocess.call('python3 new_BLAST_results_FeatureGenerator.py -i {} -o {}'.format(input_dir + 'blastp62/', input_dir + 'blastp62_features/'), shell=True)
    subprocess.call('python3 new_BLAST_results_FeatureGenerator.py -i {} -o {}'.format(input_dir + 'blastp45/', input_dir + 'blastp45_features/'), shell=True)
    subprocess.call('python3 deltaBLAST_results_FeatureGenerator.py -i {} -o {}'.format(input_dir + 'deltaBlast/', input_dir + 'deltaBlast_features/'), shell=True)
    subprocess.call('python3 process_interproscan.py -i {} -o {} -m {}'.format(input_dir + 'interproscan/', input_dir + 'interproscan_features/', os.getcwd() + '/GO_mapping/'), shell=True)

    input_features = out_folder + 'features/'
    create_dir(input_features)    
    copy_feature_files(out_folder)

    
def argument_parser():
    parser = argparse.ArgumentParser(description="targetP prediction")
    parser.add_argument("-o", "--out_folder", type=str, help="output of the folder")
    parser.add_argument("-s", "--ori_seq", type=str, help="input sequence from the user")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    """
    missing the fallback system
    # mkdir $main_dir/fallback
    # mkdir $main_dir/fallback/targetFiles
    # mkdir $main_dir/fallback/trial

    # python CAFA2_fallback.py -d $main_dir/fallback/ -g /home/kahaka/CAFA3/data_save/Swissprot_evidence.tsv.gz -i /home/kahaka/CAFA3/data_save/Swissprot_propagated.tsv.gz -o /home/kahaka/CAFA3/data_save/parent_term2term.txt.gz -p $main_dir/cafa-predictions.tsv.gz -f Jari_develtest.tsv.gz -b /home/sukaew/CAFA_QA/features/CAFA2/CAFA3_features/ -s targetFiles/

    # python generate_Kai_input.py -i $input_seq/ -o $input_dir/
    """
    args = argument_parser()
    main(args.out_folder, args.ori_seq)
