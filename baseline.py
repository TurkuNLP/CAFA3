"""
Naive blast baseline for comparison.
"""
import gzip
from collections import defaultdict

def predict(sequences, blast_dict, annotation_path):
    ann_file = gzip.open(annotation_path)
    ann_data = ann_file.readlines()
    annotations = defaultdict(list)
    for line in ann_data:
        prot_id, go_id, evidence = line.strip().split('\t')
        annotations[prot_id].append(go_id)
    
    res = []
    for s in sequences:
        blast_hits = blast_dict[s]
        ann_set = set()
        for b, score in blast_hits:
            b_anns = annotations[b]
            ann_set.update(b_anns)
        res.append((s, list(ann_set)))
    #import pdb; pdb.set_trace()
    return res