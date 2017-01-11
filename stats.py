"""
Calculates simple statistics from the Swissprot and GO data.
"""
import gzip
import numpy
from collections import defaultdict
from itertools import izip
import sys

def stats(ann_path='./data/Swissprot_evidence.tsv.gz'):
    seq_file = gzip.open('./data/Swissprot_sequence.tsv.gz')
    seq_data = seq_file.readlines()
    
    ann_file = gzip.open(ann_path)
    ann_data = ann_file.readlines()
    
    annotations = defaultdict(int)
    evidences = defaultdict(int)
    proteins = defaultdict(int)
    
    for line in ann_data:
        prot_id, annotation, evidence = line.strip().split('\t')
        annotations[annotation] += 1
        evidences[evidence] += 1
        proteins[prot_id] +=1

        
    print 'Unique annotations: %s, occurrences: %s, mean: %s, median: %s, max: %s' % (len(annotations.keys()), sum(annotations.values()), numpy.mean(annotations.values()), numpy.median(annotations.values()), max(annotations.values()))
    print 'Unique proteins: %s' % len(proteins)
    print 'Annotations per protein, mean: %s, median: %s, max: %s' % (numpy.mean(proteins.values()), numpy.median(proteins.values()), max(proteins.values()))
    
    print 'Top N annotations cover (Nth only):'
    for i in [1, 100, 1000, 2500, 4000, 5000, 10000, 20000]:
        top = sum(sorted(annotations.values(), reverse=True)[:i])
        nth = sorted(annotations.values(), reverse=True)[i-1]
        print '%s: %s, (%s)' % (i, float(top)/sum(annotations.values()), nth)
    
    print 'Evidence counts:'
    for ev, count in sorted(evidences.items()):
        print '%s\t%s' % (ev, count)
    
    seq_lengths = []
    characters = set()
    
    for prot_id, seq in pairwise(seq_data):
        prot_id = prot_id.strip().strip('>')
        seq = seq.strip()
        characters.update(seq)
        seq_lengths.append(len(seq))
        #if len(seq) > 30000:
        #    import pdb; pdb.set_trace()
    
    print 'Unique amino acids: %s' % len(characters)
    print 'Sequence length, mean: %s, median: %s, min: %s, max: %s' % (numpy.mean(seq_lengths), numpy.median(seq_lengths), min(seq_lengths), max(seq_lengths))
    print 'Acids: %s' % ''.join(sorted(characters))
    
    print 'Sequences of length N or smaller cover:'
    for i in [10, 100, 500, 1000, 2500, 5000, 10000]:
        seq_count = len([1 for s in seq_lengths if s <= i])
        print '%s\t%s' % (i, float(seq_count)/len(seq_lengths))

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)
    
if __name__ == '__main__':
    stats(sys.argv[1])