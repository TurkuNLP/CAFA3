from keras.layers import Input, LSTM, RepeatVector, Embedding, TimeDistributed, Dense, Convolution1D, MaxPooling1D, GlobalMaxPooling1D, GlobalAveragePooling1D, Flatten, merge, Masking
from keras.models import Model
from keras.preprocessing import sequence
from keras.utils import np_utils
from keras.callbacks import EarlyStopping, ModelCheckpoint

import os
import gzip
import numpy as np
np.random.seed(1337)
import codecs
import cPickle as pickle
from collections import defaultdict

from stats import pairwise

timesteps = 2500 # maximum length of a sequence, the real max is 35K. 2.5K covers 99% of the sequences, 5K 99.9%
latent_dim = 50 # Amino acid embedding size
char_set = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
vocab_size = len(char_set) + 1 # +1 for mask
char_dict = {c:i+1 for i,c in enumerate(char_set)} # Index 0 is left for padding
use_features = False # False = only sequence is used for prediction
model_dir = './model/' # path for saving model + other required stuff

# TODO: Save model + other needed files after training: model, go_id mapping, blast id map
# TODO: Add normal features
# TODO: Get prediction statistics
# TODO: Word dropout?
# TODO: Test embedding sizes
# FIXME: CAFA targets have * character which should be added to the embeddings (OOV character)
# TODO: Make predictions for test + cafa targets


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

def read_split_ids(split_path):
    split_file = gzip.open(split_path)
    split_data = set([s.strip() for s in split_file])
    return split_data

def generate_data(split_path, seq_path, ann_path, ann_ids, batch_size=256):
    """
    Generates NN compatible data.
    """
    seq_file = gzip.open(seq_path)
    seq_data = seq_file.readlines()
    
    ann_file = gzip.open(ann_path)
    ann_data = ann_file.readlines()
    ann_dict = get_annotation_dict(ann_data)
    
    if split_path:
        split_ids = read_split_ids(split_path)
    
    while True:
        prot_ids = []
        x = []
        blast_x = []
        y = []
        
        for i, (prot_id, seq) in enumerate(pairwise(seq_data)):
            #if i % 10000 == 0:
            #    print i
            prot_id = prot_id.strip().strip('>')
            if split_path and prot_id not in split_ids:
                continue
            prot_ids.append(prot_id)
            seq = seq.strip()
            seq_id_list = [char_dict[s] for s in seq]
            annotations = ann_dict[prot_id]
            ann_id_list = [ann_ids[a] for a in annotations if a in ann_ids]
            y_v = np.zeros((len(ann_ids)), dtype='int')
            y_v[ann_id_list] = 1
            x.append(seq_id_list)
            if use_features:
                blast_x.append(generate_blast_features(prot_id))
            y.append(y_v)
        
            if len(prot_ids) == batch_size:
                x = sequence.pad_sequences(x, timesteps)
                blast_x = np.array(blast_x)
                y = np.array(y, dtype='int32')

                nn_data = {'sequence': x, 'labels': y, 'features': blast_x, 'prot_ids': np.array(prot_ids)}
                yield nn_data, nn_data
                
                prot_ids = []
                x = []
                blast_x = []
                y = []
        
        if len(prot_ids) > 0:
            x = sequence.pad_sequences(x, timesteps)
            blast_x = np.array(blast_x)
            y = np.array(y)
        
            nn_data = {'sequence': x, 'labels': y, 'features': blast_x, 'prot_ids': np.array(prot_ids)}
            yield nn_data, nn_data

def generate_blast_data():
    """
    Creates blast features for the given sequences.
    """
    print 'Reading blast info'
    blast_f = gzip.open('./data/Swissprot_blast_filtered.tsv.gz')
    blast_data = blast_f.readlines()
    blast_dict = defaultdict(list)
    blast_hits = set()
    #blast_prot = set([line.strip().split('\t')[0] for line in blast_data])
    for e, line in enumerate(blast_data):
        if e % 1000000 == 0:
            print e
        data = line.strip().split('\t')
        prot_id = data[0]
        hit = data[4]
        #if hit not in blast_prot:
        #    continue
        score = float(data[9])
        blast_dict[prot_id].append((hit, score))
        blast_hits.add(hit)
        
    blast_hit_ids = {hit: i for i, hit in enumerate(list(sorted(blast_hits)))}
    
    return blast_dict, blast_hit_ids

#if use_features:
blast_dict, blast_hit_ids = generate_blast_data()

def generate_blast_features(prot_id):
    x = np.zeros((len(blast_hit_ids)))
    for hit, score in blast_dict[prot_id]:
        # TODO: Get GO's to transfer
        x[blast_hit_ids[hit]] = score
    return x

def go_to_ids(predictions, ann_ids):
    y = []
    for p in predictions:
        ann_id_list = [ann_ids[a] for a in p if a in ann_ids]
        y_v = np.zeros((len(ann_ids)), dtype='int')
        y_v[ann_id_list] = 1
        y.append(y_v)
    return np.array(y)

def _data_size(path):
    return len(gzip.open(path).readlines())

def train():
    print 'Generating training data'
    ann_path = './data/Swissprot_propagated.tsv.gz'
    ann_ids, reverse_ann_ids = get_annotation_ids(ann_path, top=4000)
    batch_size=16
    #pretrain_data = generate_data(None, '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz', '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ann-train-data.tsv.gz', ann_ids, 256)
    #pretrain_size = _data_size('/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz')/2
    train_data = generate_data('./data/train.txt.gz', './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    train_size = _data_size('./data/train.txt.gz')
    devel_data = generate_data('./data/devel.txt.gz', './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    devel_size = _data_size('./data/devel.txt.gz')
    #test_data = generate_data('./data/test.txt.gz', './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids)
    #test_size = _data_size('./data/test.txt.gz')
    
    #print "Making baseline predictions"
    #import baseline
    #devel_baseline = baseline.predict(devel_data['prot_ids'], blast_dict, ann_path)
    #devel_baseline_ids = go_to_ids([b[1] for b in devel_baseline], ann_ids)
    #from sklearn import metrics
    #baseline_score = metrics.precision_recall_fscore_support(devel_data['labels'], devel_baseline_ids, average='micro')
    #print 'Baseline score: ', baseline_score
    
    #import pdb; pdb.set_trace()
    #for ii in [3, 6, 9, 15, 27, 50]:
    #    print '### Testing window size %s' % ii
    print 'Building model'
    inputs = Input(shape=(timesteps, ), name='sequence')
    input_list = [inputs]
    embedding = Embedding(vocab_size, latent_dim, mask_zero=False)(inputs)
    
    #mask = Masking()(embedding)
    
    convs = []
    
    # Stacked CNN experiments
    #encoded = Convolution1D(50, 3, border_mode='valid', activation='linear')(embedding)
    ##maxed = GlobalMaxPooling1D()(encoded)
    ##convs.append(maxed)
    #encoded = Convolution1D(50, 3, border_mode='valid', activation='linear')(encoded)
    ##maxed = GlobalMaxPooling1D()(encoded)
    ##convs.append(maxed)
    #encoded = Convolution1D(50, 3, border_mode='valid', activation='linear')(encoded)
    #encoded = GlobalMaxPooling1D()(encoded)
    #convs.append(maxed)
    
    for i in [3, 9, 27, 81]:
        encoded = Convolution1D(50, i, border_mode='valid', activation='relu')(embedding)
        encoded = GlobalMaxPooling1D()(encoded)
        convs.append(encoded)

        ## LSTM attention
        #lstm = LSTM(50)(mask)
        ##convs.append(lstm)
        #
        #from attention import Attention
        #att = Attention()([encoded, lstm])
        #convs.append(att)
    
    if use_features:
        feature_input = Input(shape=(len(blast_hit_ids), ), name='features')
        feature_encoding = Dense(300, activation='tanh')(feature_input) # Squeeze the feature vectors to a tiny encoding
        convs.append(feature_encoding)
        input_list.append(feature_input)
    #
    #encoded = feature_encoding
    encoded = merge(convs, mode='concat')
    
    predictions = Dense(len(ann_ids), activation='sigmoid', name='labels')(encoded)
    
    model = Model(input_list, predictions)
    model.compile(optimizer='adam', loss=weighted_binary_crossentropy, metrics=['accuracy', 'precision', 'recall', 'fmeasure'])
    print model.summary()
    
    print 'Training model'
    es_cb = EarlyStopping(monitor='val_fmeasure', patience=10, verbose=1, mode='max')
    cp_cb = ModelCheckpoint(filepath=os.path.join(model_dir, 'model.hdf5'), monitor='val_fmeasure', save_best_only=True,verbose=1)
    #model.fit(pretrain_data, pretrain_data, nb_epoch=100, batch_size=16, validation_data=[devel_data, devel_data], callbacks=[es_cb, cp_cb])
    model.fit_generator(train_data, samples_per_epoch=train_size, nb_epoch=100, validation_data=devel_data, nb_val_samples=devel_size, callbacks=[es_cb, cp_cb])

        
    import pdb; pdb.set_trace()
    
def weighted_binary_crossentropy(target, output):
    from keras.backend.common import _EPSILON
    from keras import backend as K
    from theano import tensor as T
    from theano.tensor import basic as tensor
    pos_weight = 3.0 # 73 is roughly the inverse ratio of positives examples
    output = T.clip(output, _EPSILON, 1.0 - _EPSILON)
    ce = -(pos_weight * target * tensor.log(output) + (1.0 - target) * tensor.log(1.0 - output))
    return K.mean(ce, axis=-1)
    
if __name__ == '__main__':
    train()