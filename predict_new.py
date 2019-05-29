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
import pickle as pickle
from collections import defaultdict
import argparse

from stats import pairwise
from train import generate_data, _get_ids, read_feature_json, get_annotation_ids, get_annotation_dict, read_split_ids, _data_size, save_predictions, weighted_binary_crossentropy, use_features, SEQUENCE_PATH
import train

# FIXME: Settings of the model and train.py should always match
batch_size = 100
# model_dir = './develtestnotax_model/' # path for saving model + other required stuff
    #import pdb; pdb.set_trace()
    
def predict(model_dir, sequence_path, feature_path, out_path='./asdf_predictions/predictions.tsv.gz'):


    if use_features:
        json_vectorizer = pickle.load(open(os.path.join(model_dir, 'json_vectorizer.pkl'), 'rb'))
        feature_selector = pickle.load(open(os.path.join(model_dir, 'feature_selector.pkl'), 'rb'))
        json_feature_matrix, json_id_map, _, _ = read_feature_json(path=feature_path, vectorizer=json_vectorizer, feature_selector=feature_selector)
        train.json_vectorizer = json_vectorizer
        train.feature_selector = feature_selector
        train.json_feature_matrix = json_feature_matrix
        train.json_id_map = json_id_map
    
    print('Generating data')
    reverse_ann_ids = pickle.load(open(os.path.join(model_dir, 'reverse_ann_ids.pkl'), 'rb'))
    
    devel_data = generate_data(None, sequence_path, None, None, batch_size)
    devel_ids = _get_ids(generate_data(None, sequence_path, None, None, batch_size, endless=False))
    
    print('Loading model')

    from keras.models import load_model
    model = load_model(filepath=os.path.join(model_dir, 'model.h5'), custom_objects={"weighted_binary_crossentropy":weighted_binary_crossentropy})
    
    print(model.summary())
    
    print("Making predictions")


    devel_pred = model.predict_generator(devel_data, batch_size, verbose=1)
    
    if not os.path.exists(os.path.dirname(out_path)):
        os.makedirs(out_path)
    save_predictions(out_path, devel_ids, devel_pred, reverse_ann_ids)
    
    print('All done.')
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Predict GO terms for protein sequences.')
    parser.add_argument('model', help='Path to model directory')
    parser.add_argument('sequences', help='Path to sequence file')
    parser.add_argument('features', help='Path to feature json file')
    parser.add_argument('output', help='Path for output file (should be .tsv.gz)')
    
    args = parser.parse_args()
    # import pdb; pdb.set_trace()
    
    predict(args.model, args.sequences, args.features, args.output)
