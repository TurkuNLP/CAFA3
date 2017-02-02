"""
Sanity check between Kai's format and official format.
"""
import os
import sys
import gzip
import csv
from collections import defaultdict

def check(kai_file, farrokh_path):
    kai_pred = read_predictions(kai_file)
    
    files = [os.path.join(farrokh_path, 'TurkuBioNLP1_2_ProteinCentric', f) for f in os.listdir(os.path.join(farrokh_path, 'TurkuBioNLP1_2_ProteinCentric'))]
    files += [os.path.join(farrokh_path, 'TurkuBioNLP1_2_TermCentricAndMoonlighting', f) for f in os.listdir(os.path.join(farrokh_path, 'TurkuBioNLP1_2_TermCentricAndMoonlighting'))]
    
    submission_pred = build_submission_pred(files)
    
    print sorted(kai_pred) == sorted(submission_pred)
    
def build_submission_pred(files):
    res = []
    for f in files:
        lines = open(f).readlines()[3:]
        for i, li in enumerate(lines):
            if li == 'END':
                break
            if i % 100000 == 0:
                print i
            try:
                pid, gid, conf = li.strip().split('\t')
            except:
                import pdb; pdb.set_trace()
            res.append((pid, gid, conf))
    return res

def read_predictions(path):
    pred = []
    with gzip.open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for i, row in enumerate(reader):
            if i == 0: # Skip header line
                continue
            if i % 100000 == 0:
                print i
            prot_id = row[5]
            go_id = row[2]
            conf = row[4]
            pred.append((prot_id, go_id, conf))

    return pred

if __name__ == '__main__':
    check(sys.argv[1], sys.argv[2])