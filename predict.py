from keras.layers import Dropout, Input, LSTM, RepeatVector, Embedding, TimeDistributed, Dense, Convolution1D, MaxPooling1D, GlobalMaxPooling1D, GlobalAveragePooling1D, Flatten, merge, Masking
from keras.models import Model
from keras.preprocessing import sequence
from keras.utils import np_utils
from keras.callbacks import EarlyStopping, ModelCheckpoint

import os
import gzip
import json
import numpy as np
np.random.seed(1337)
import codecs
import cPickle as pickle
from collections import defaultdict

from stats import pairwise
from train import generate_data, _get_ids, read_feature_json, get_annotation_ids, get_annotation_dict, read_split_ids, _data_size, save_predictions, weighted_binary_crossentropy


# FIXME: These should always match to the model and should be saved as a config file
ann_limit = 5000 # Taking top N GO annotations only
timesteps = 2500 # maximum length of a sequence, the real max is 35K. 2.5K covers 99% of the sequences, 5K 99.9%
latent_dim = 50 # Amino acid embedding size
batch_size = 16
char_set = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
vocab_size = len(char_set) + 1 # +1 for mask
char_dict = {c:i+1 for i,c in enumerate(char_set)} # Index 0 is left for padding
use_features = False # False = only sequence is used for prediction
model_dir = './cnn2_model/' # path for saving model + other required stuff
out_dir = model_dir
#out_dir = './test_predictions/'
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

if use_features:
    json_vectorizer = pickle.load(open(os.path.join(model_dir, 'json_vectorizer.pkl', 'rb')))
    feature_selector = pickle.load(open(os.path.join(model_dir, 'feature_selector.pkl', 'rb')))
    json_feature_matrix, json_id_map, _, _ = read_feature_json(vectorizer=json_vectorizer, feature_selector=feature_selector)

def predict():
    print 'Generating data'
    ann_path = './data/Swissprot_propagated.tsv.gz'
    ann_ids = pickle.load(open(os.path.join(model_dir, 'ann_ids.pkl', 'rb')))
    reverse_ann_ids = pickle.load(open(os.path.join(model_dir, 'reverse_ann_ids.pkl', 'rb')))
    
    devel_path = './data/devel.txt.gz'
    devel_data = generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    devel_size = _data_size(devel_path)
    devel_ids = _get_ids(generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size, endless=False))
    test_path = './data/test.txt.gz'
    test_data = generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids)
    test_size = _data_size(test_path)
    test_ids = _get_ids(generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size, endless=False))
    
    print 'Loading model'

    from keras.models import load_model
    model = load_model(filepath=os.path.join(model_dir, 'model.hdf5'), custom_objects={"weighted_binary_crossentropy":weighted_binary_crossentropy})
    
    print model.summary()
    
    print "Making predictions"
    
    devel_score = model.evaluate_generator(devel_data, devel_size)
    test_score = model.evaluate_generator(test_data, test_size)
    print 'Devel l/a/p/r/f: ', devel_score
    print 'Test l/a/p/r/f: ', test_score
    
    # Reinstantiate the data generators, otherwise they are not correctly aligned anymore
    devel_data = generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    test_data = generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids)


    devel_pred = model.predict_generator(devel_data, devel_size)
    test_pred = model.predict_generator(test_data, test_size)
    
    save_predictions(os.path.join(out_dir, 'devel_pred.tsv.gz'), devel_ids, devel_pred, reverse_ann_ids)
    save_predictions(os.path.join(out_dir, 'test_pred.tsv.gz'), test_ids, test_pred, reverse_ann_ids)
    
    print 'Making CAFA target predictions'
    
    cafa_id_path = './data/target.all.ids.gz'
    cafa_seq_path = '/home/sukaew/CAFA3/CAFA3_targets/target.all.fasta.gz'
    cafa_data = generate_data(None, cafa_seq_path, ann_path, ann_ids, batch_size=batch_size, cafa_targets=True, verbose=True)
    cafa_size = _data_size(cafa_id_path)
    cafa_ids = _get_ids(generate_data(None, cafa_seq_path, ann_path, ann_ids, batch_size=batch_size, cafa_targets=True, verbose=False, endless=True))
    #import pdb; pdb.set_trace()
    cafa_pred = model.predict_generator(cafa_data, cafa_size)
    
    save_predictions(os.path.join(out_dir, 'cafa_targets.tsv.gz'), cafa_ids, cafa_pred, reverse_ann_ids, cafa_targets=True)
    #import pdb; pdb.set_trace()
    
    print 'All done.'
    
if __name__ == '__main__':
    predict()