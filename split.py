"""
This script is for splitting the whole data to train/devel/test sets.
"""
import sys
import gzip

from stats import pairwise
from sklearn.model_selection import train_test_split

def split(seq_path):
    seq_file = gzip.open(seq_path)
    seq_data = seq_file.readlines()
    
    ids = [prot_id.strip().strip('>') for prot_id, seq in pairwise(seq_data)]
    
    train, devtest = train_test_split(ids, test_size=0.4) # 60% for training, 40% for devel and test
    
    devel, test = train_test_split(devtest, test_size=0.5) # 50% for devel (20% of the whole data), 50% for test
    
    splits = [train, devel, test]
    
    for split_name, split in zip(['train', 'devel', 'test'], splits):
        f = gzip.open('./data/%s.txt.gz' % split_name, 'w')
        for prot_id in split:
            f.write('%s\n' % prot_id)
        f.close()

if __name__ == '__main__':
    split(sys.argv[1])