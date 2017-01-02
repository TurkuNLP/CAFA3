from keras.layers import Input, LSTM, RepeatVector, Embedding, TimeDistributed, Dense, Convolution1D, MaxPooling1D, GlobalMaxPooling1D, Flatten, merge
from keras.models import Model
from keras.preprocessing import sequence
from keras.utils import np_utils

import gzip
import numpy as np
np.random.seed(1337)
import codecs
import cPickle as pickle
from collections import defaultdict

from stats import pairwise

timesteps = 5000 # maximum length of a sequence, the real max is 35K
latent_dim = 50 # Amino acid embedding size
char_set = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
vocab_size = len(char_set) + 1 # +1 for mask
char_dict = {c:i+1 for i,c in enumerate(char_set)} # Index 0 is left for padding

# TODO: Add fixed train/devel/test sets
# TODO: Add F-score evaluation and early stopping
# TODO: Save model + other needed files after training
# TODO: Try LSTM
# TODO: Try some fancy pooling approach
# TODO: Add normal features

def get_annotation_ids(annotation_path, top=None):
    """
    Maps GO ids to integers.
    """
    ann_file = gzip.open(annotation_path)
    ann_data = ann_file.readlines()
    annotations = defaultdict(int)
    for line in ann_data:
        annotations[line.strip().split('\t')[1]] += 1
    
    if top:
        annotations = [a[0] for a in sorted(annotations.items(), key=lambda x: x[1], reverse=True)[:top]]
    else:
        annotations = annotations.keys()
    
    
    ann_ids = {a: i for i, a in enumerate(annotations)}
    reverse_ann_ids = {i: a for a, i in ann_ids.items()}
    
    return ann_ids, reverse_ann_ids

def get_annotation_dict(annotation_data):
    """
    Maps sequence ids to a list of GO ids.
    """
    ann_dict = defaultdict(list)
    for line in annotation_data:
        prot_id, annotation, evidence = line.strip().split('\t')
        ann_dict[prot_id].append(annotation)
    return ann_dict

def generate_data(seq_path, ann_path, ann_ids):
    """
    Generates NN compatible data.
    """
    seq_file = gzip.open(seq_path)
    seq_data = seq_file.readlines()
    
    ann_file = gzip.open(ann_path)
    ann_data = ann_file.readlines()
    ann_dict = get_annotation_dict(ann_data)
    
    x = []
    y = []
    
    for i, (prot_id, seq) in enumerate(pairwise(seq_data)):
        if i % 10000 == 0:
            print i
        prot_id = prot_id.strip().strip('>')
        seq = seq.strip()
        seq_id_list = [char_dict[s] for s in seq]
        annotations = ann_dict[prot_id]
        ann_id_list = [ann_ids[a] for a in annotations if a in ann_ids]
        y_v = np.zeros((len(ann_ids)), dtype='int')
        y_v[ann_id_list] = 1
        x.append(seq_id_list)
        y.append(y_v)
        
    #import pdb; pdb.set_trace()
    x = sequence.pad_sequences(x, timesteps)
    y = np.array(y)
    
    return {'sequence': x, 'labels': y}

def train():
    print 'Generating training data'
    ann_ids, reverse_ann_ids = get_annotation_ids('./data/Swissprot_propagated.tsv.gz', top=5000)
    train_data = generate_data('./data/Swissprot_sequence.tsv.gz', './data/Swissprot_propagated.tsv.gz', ann_ids)
    
    print 'Building model'
    inputs = Input(shape=(timesteps, ), name='sequence')
    embedding = Embedding(vocab_size, latent_dim, mask_zero=False)(inputs)
    
    convs = []
    for i in [3, 9, 27, 81]:
        encoded = Convolution1D(50, i, border_mode='valid')(embedding)
        encoded = GlobalMaxPooling1D()(encoded)
        convs.append(encoded)
    encoded = merge(convs, mode='concat')
    
    predictions = Dense(len(ann_ids), activation='sigmoid', name='labels')(encoded)
    
    model = Model(inputs, predictions)
    model.compile(optimizer='adam', loss='categorical_crossentropy', metrics=['accuracy'])
    print model.summary()
    
    print 'Training model'
    model.fit(train_data, train_data, nb_epoch=100, batch_size=128, validation_split=0.2)
            
    import pdb; pdb.set_trace()
    
if __name__ == '__main__':
    train()