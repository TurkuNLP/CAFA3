"""
Evaluation for sanity check.
This evaluates the REAL score, not just top N.
"""
import csv
import os
import gzip
import sys
from collections import defaultdict
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.metrics import precision_recall_fscore_support

from train import get_annotation_dict, read_split_ids


def evaluate(pred_path):
    ann_path = './data/CAFA_PI_Swissprot_propagated.tsv.gz'
    ann_file = gzip.open(ann_path)
    ann_data = ann_file.readlines()
    ann_dict = get_annotation_dict(ann_data)
    
    devel_path = './data/devel.txt.gz'
    devel_ids = read_split_ids(devel_path, unique=True)
    
    test_path = './data/test.txt.gz'
    test_ids = read_split_ids(test_path, unique=True)

    devel_ann = {k: v for k, v in ann_dict.items() if k in devel_ids}
    
    test_ann = {k: v for k, v in ann_dict.items() if k in test_ids}
    
    devel_pred = read_predictions(os.path.join(pred_path, 'devel_pred.tsv.gz'))
    
    test_pred = read_predictions(os.path.join(pred_path, 'test_pred.tsv.gz'))

    # Not the real order of the annotations, but it doesn't matter as long as the gold/pred are consistent
    lb = MultiLabelBinarizer(sparse_output=True).fit(ann_dict.values())
    # import pdb; pdb.set_trace()
    print 'Devel results:'
    print _eval(lb, devel_ids, devel_pred, devel_ann)
    
    print 'Test results:'
    print _eval(lb, test_ids, test_pred, test_ann)
    
def _eval(lb, ids, pred, gold):
    gold_vectors = lb.transform([gold.get(i, []) for i in ids])
    print gold_vectors.shape, gold_vectors.shape[0]*gold_vectors.shape[1], gold_vectors.count_nonzero()
    
    pred_vectors = lb.transform([pred.get(i, []) for i in ids])
    # import pdb; pdb.set_trace()
    return precision_recall_fscore_support(gold_vectors, pred_vectors, average='micro')
    
def read_predictions(path):

    pred_dict = defaultdict(list)
    with gzip.open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for i, row in enumerate(reader):
            if i == 0: # Skip header line
                continue
            prot_id = row[0]
            go_id = row[2]
            pred_dict[prot_id].append(go_id)

    return pred_dict

if __name__ == '__main__':
    evaluate(sys.argv[1])