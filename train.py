from keras.layers import Dropout, Input, LSTM, RepeatVector, Embedding, TimeDistributed, Dense, Convolution1D, MaxPooling1D, GlobalMaxPooling1D, GlobalAveragePooling1D, Flatten, concatenate, Masking
from keras.models import Model
from keras.preprocessing import sequence
from keras.utils import np_utils
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.optimizers import Adam
import keras

import os
import gzip
import json
import numpy as np
np.random.seed(1337)
import pickle as pickle
from collections import defaultdict

from stats import pairwise

SEQUENCE_PATH = './data/Swissprot_sequence.tsv.gz'
ann_path = './data/Swissprot_propagated.tsv.gz'

ann_limit = 5000 # Taking top N GO annotations only
timesteps = 2500 # maximum length of a sequence, the real max is 35K. 2.5K covers 99% of the sequences, 5K 99.9%
latent_dim = 50 # Amino acid embedding size
batch_size = 1000 # Warning: this is actually the number of batches in the new Keras API
char_set = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
vocab_size = len(char_set) + 1 # +1 for mask
char_dict = {c:i+1 for i,c in enumerate(char_set)} # Index 0 is left for padding
use_features = True # False = only sequence is used for prediction
model_dir = './cnn_only/' # path for saving model + other required stuff
if not os.path.exists(model_dir):
    os.makedirs(model_dir)

# TODO: Make predictions for CNN2 and Full3 models
# TODO: Evaluate prediction files
# TODO: Do a sanity check for CAFA target predictions

# TODO: Get prediction statistics
# TODO: Word dropout?
# TODO: Remove obvious GOs?

def get_annotation_ids(annotation_path, top=None):
    """
    Maps GO ids to integers.
    """
    ann_file = gzip.open(annotation_path, 'rt')
    ann_data = ann_file.readlines()
    annotations = defaultdict(int)
    for line in ann_data:
        # import pdb; pdb.set_trace()
        annotations[line.strip().split('\t')[1]] += 1
    
    if top:
        annotations = [a[0] for a in sorted(list(annotations.items()), key=lambda x: x[1], reverse=True)[:top]]
    else:
        annotations = list(annotations.keys())
    
    
    ann_ids = {a: i for i, a in enumerate(annotations)}
    reverse_ann_ids = {i: a for a, i in list(ann_ids.items())}
    
    return ann_ids, reverse_ann_ids

ann_ids, reverse_ann_ids = get_annotation_ids(ann_path, top=ann_limit)

def get_annotation_dict(annotation_data):
    """
    Maps sequence ids to a list of GO ids.
    """
    ann_dict = defaultdict(list)
    for line in annotation_data:
        prot_id, annotation, evidence = line.strip().split('\t')
        ann_dict[prot_id].append(annotation)
    return ann_dict

def read_split_ids(split_path, unique=True):
    split_file = gzip.open(split_path, 'rt')
    split_data = [s.strip().replace('>', '') for s in split_file]
    if unique:
        split_data = set(split_data)
    return split_data

def generate_data(split_path, seq_path, ann_path, ann_ids, batches=125, cafa_targets=False, verbose=False, endless=True):
    """
    Generates NN compatible data.
    """
    seq_file = gzip.open(seq_path, 'rt')
    seq_data = seq_file.readlines()
    
    sequence_count = len(seq_data) / 2
    
    if ann_path:
        ann_file = gzip.open(ann_path, 'rt')
        ann_data = ann_file.readlines()
        ann_dict = get_annotation_dict(ann_data)
    else:
        ann_dict = defaultdict(list)
    
    if split_path:
        split_ids = read_split_ids(split_path)
        sequence_count = len(split_ids)
    else:
        split_ids = None
    
    sequences = {seq_id.strip().replace('>', ''): seq for seq_id, seq in pairwise(seq_data)}
    
    if split_ids:
        sequence_ids = list(split_ids)
    else:
        sequence_ids = list(sequences.keys())
    print('Data size: %s' % len(sequence_ids))
    # import pdb; pdb.set_trace()
    while True:
        for i, batch in enumerate(np.array_split(sequence_ids, batches)):
            if verbose:
                print(i, len(batch))
            prot_ids = []
            x = []
            blast_x = []
            y = []
            for prot_id in batch:
                try:
                    seq = sequences[prot_id]
                except:
                    import pdb; pdb.set_trace()
                # if verbose and i % 10000 == 0:
                #     print 'Generator: ', i
                if cafa_targets:
                    cafa_id, prot_id = prot_id.strip().replace('>', '').split(' ')
                else:
                    prot_id = prot_id.strip().strip('>')
                # if split_path and prot_id not in split_ids:
                #     continue
                if cafa_targets:
                    prot_ids.append('%s %s' % (cafa_id, prot_id))
                else:
                    prot_ids.append(prot_id)
                seq = seq.strip()
                seq_id_list = [char_dict[s] for s in seq if s != '*']
                #seq_id_list = [aa_index_ids.get(s, 0) for s in seq]
                if ann_ids:
                    annotations = ann_dict[prot_id]
                    ann_id_list = [ann_ids[a] for a in annotations if a in ann_ids]
                    y_v = np.zeros((len(ann_ids)), dtype='int')
                    y_v[ann_id_list] = 1
                    y.append(y_v)
                x.append(seq_id_list)
                if use_features:
                    #blast_x.append(generate_blast_features(prot_id)) # These are our original blast features
                    blast_x.append(get_feature_vector(prot_id)) # Jari's feature vectors
                
            
            x = sequence.pad_sequences(x, timesteps)
            blast_x = np.array(blast_x)
            y = np.array(y, dtype='int32')
            
            nn_data = {'sequence': x, 'labels': y, 'features': blast_x, 'prot_ids': np.array(prot_ids)}
            yield nn_data, nn_data
            
        if not endless:
            break

# def read_aaindex():
#     aa_f = open('/home/sukaew/CAFA3/aaindex/aaindex_table.tsv')
#     data = aa_f.readlines()
#     amino_acids = data[0].strip().split('\t')[1:]
#     aa_index_ids = {a: i+1 for i, a in enumerate(amino_acids)}
#     
#     embedding = np.zeros((len(amino_acids)+1, len(data)-1)) # +1 for OOV/MASK no separation for now
#     for i, row in enumerate(data[1:]):
#         values = row.strip().split('\t')[1:]
#         for ii, v in enumerate(values):
#             if v == 'NA': # FIXME: What to do with these?
#                 v = 0.0
#             embedding[ii+1, i] = float(v)
#         
#     #import pdb; pdb.set_trace()
#     return aa_index_ids, embedding
# 
# aa_index_ids, aa_embedding = read_aaindex()

def read_feature_json(path='./data/examples.json.gz', vectorizer=None, feature_selector=None):
    print('Reading feature data')
    js = json.load(gzip.open(path, 'rt'))
    # import pdb; pdb.set_trace()
    filters = ['DUMMY']#['BLAST', 'DELTA', 'GPI', 'TAX', 'IPS']
    print("Excluding: ", filters)
    for i, d in enumerate(js['features']):
        if i % 10000 == 0:
            print(i)
        for f in filters:
            for key in list(d.keys()):
                if key.startswith('%s:' % f):
                    d.pop(key)
    
    
    from sklearn.feature_extraction import DictVectorizer
    if vectorizer:
        v = vectorizer
    else:
        v = DictVectorizer()
        
        # Only get features that exist in training examples to densify the feature space, all the rest are useless anyway
        train_ids = set(read_split_ids('./data/train.txt.gz', unique=False))
        train_features = []
        for prot_id, features in zip(js['ids'], js['features']):
            if prot_id in train_ids:
                train_features.append(features)
        # import pdb; pdb.set_trace()
        v.fit(train_features)

    feature_matrix = v.transform(js['features'])
    
    id_map = {pid: i for i, pid in enumerate(js['ids'])}
    
    from sklearn import preprocessing
    scaler = preprocessing.MaxAbsScaler().fit(feature_matrix)
    std_matrix = scaler.transform(feature_matrix)
    
    from sklearn.feature_selection import VarianceThreshold
    #import pdb; pdb.set_trace()
    if feature_selector:
        vt = feature_selector
    else:
        vt = VarianceThreshold(0.0001).fit(std_matrix)
    
    better_matrix = vt.transform(std_matrix)
    
    # better_matrix = std_matrix # Bypass feature selection
    #good_features = np.array(v.feature_names_)[np.where(vt.transform(std_matrix)==True)]
    
    
    return better_matrix, id_map, v, vt

if use_features:
    json_feature_matrix, json_id_map, json_vectorizer, feature_selector = read_feature_json()

def get_feature_vector(prot_id):
    prot_index = json_id_map[prot_id]
    feature_vector = json_feature_matrix[prot_index]
    return feature_vector.toarray()[0]

def generate_blast_data():
    """
    Creates blast features for the given sequences.
    """
    print('Reading blast info')
    blast_f = gzip.open('./data/Swissprot_blast_filtered.tsv.gz', 'rt')
    blast_data = blast_f.readlines()
    blast_dict = defaultdict(list)
    blast_hits = set()
    #blast_prot = set([line.strip().split('\t')[0] for line in blast_data])
    for e, line in enumerate(blast_data):
        if e % 1000000 == 0:
            print(e)
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
#    blast_dict, blast_hit_ids = generate_blast_data()

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
    return len(gzip.open(path, 'rt').readlines())

def _get_ids(data_gen):
    """
    Needs a data gen with endless=False
    """
    ids = [i[0]['prot_ids'] for i in data_gen]
    ids = list(np.concatenate(ids))
    return ids

class Evaluate(keras.callbacks.Callback):
    def __init__(self, data_path, steps, reverse_ann_ids, patience=50):
        self.data_path = data_path
        self.steps = steps
        self.reverse_ann_ids = reverse_ann_ids
        self.best = 0.0
        self.patience = patience
        self.wait = 0
        super(Evaluate, self).__init__()

    def on_epoch_end(self, epoch, logs={}):
        self.data = generate_data(self.data_path, SEQUENCE_PATH, ann_path, ann_ids, self.steps)
        gold = np.concatenate([next(self.data)[0]['labels'] for i in range(self.steps)]) # If we consume a full cycle of the generator, we should have gold labels aligned with the predictions
        # gold_labels = self._to_labels(gold)
        
        pred = self.model.predict_generator(self.data, steps=self.steps)
        pred_class = np.round(pred)
        # pred_labels = self._to_labels(pred)
        
        from sklearn.metrics import f1_score, precision_recall_fscore_support, classification_report
        
        print('')
        score = precision_recall_fscore_support(gold, pred_class, average='micro')
        print(score)
        if score[2] > self.best:
            self.best = score[2]
            # print(classification_report(gold, pred_class))
            self.model.save(os.path.join(model_dir, 'model.h5'))
            # import pdb; pdb.set_trace()
            self.wait = 0
        else:
            self.wait += 1
            if self.wait >= self.patience:
                self.model.stop_training = True
        print('\n\n')
        
    def _to_labels(self, pred):
        """
        Inverse mapping from n-hot vectors to a set of labels
        """
        res = []
        for i in range(len(pred)):
            pred_indices = np.round(pred[i]).nonzero()[0]
            seq_labels = []
            for pred_i in pred_indices:
                go_id = self.reverse_ann_ids[pred_i]
                seq_labels.append(go_id)
            res.append(seq_labels)
        return res
        
        

def train():
    print('Generating training data')

    #pretrain_data = generate_data(None, '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz', '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ann-train-data.tsv.gz', ann_ids, 256)
    #pretrain_size = _data_size('/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz')/2
    train_path = './data/train.txt.gz'
    train_data = generate_data(train_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size)
    train_size = _data_size(train_path)
    train_ids = read_split_ids(train_path, unique=False)
    # import pdb; pdb.set_trace()
    devel_path = './data/devel.txt.gz'
    devel_data = generate_data(devel_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size//10)
    devel_size = _data_size(devel_path)
    devel_ids = read_split_ids(devel_path, unique=False)
    
    test_path = './data/test.txt.gz'
    test_data = generate_data(test_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size//10)
    test_size = _data_size(test_path)
    test_ids = read_split_ids(test_path, unique=False)

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
    print('Building model')
    inputs = Input(shape=(timesteps, ), name='sequence')
    input_list = [inputs]
    embedding = Embedding(vocab_size, latent_dim, mask_zero=False)(inputs)
    embedding = Dropout(0.5)(embedding)
    #embedding = Embedding(aa_embedding.shape[0], aa_embedding.shape[1], mask_zero=False, weights=[aa_embedding], trainable=True)(inputs)
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
    
    for i in [3, 9, 27]:
        encoded = Convolution1D(400, i, padding='valid', activation='relu')(embedding)
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
        #feature_input = Input(shape=(len(blast_hit_ids), ), name='features')
        feature_input = Input(shape=(json_feature_matrix.shape[1], ), name='features') # For Jari's feature vectors
        dropout = Dropout(0.5)(feature_input)
        feature_encoding = Dense(300, activation='tanh')(dropout) # Squeeze the feature vectors to a tiny encoding
        convs.append(feature_encoding)
        input_list.append(feature_input)
    #
    #encoded = feature_encoding
    encoded = concatenate(convs)
    
    predictions = Dense(len(ann_ids), activation='sigmoid', name='labels')(encoded)
    
    model = Model(input_list, predictions)
    model.compile(optimizer=Adam(lr=0.0005), loss=weighted_binary_crossentropy, metrics=['accuracy'])
    print(model.summary())
    
    print('Training model')
    pickle.dump(ann_ids, open(os.path.join(model_dir, 'ann_ids.pkl') ,'wb'))
    pickle.dump(reverse_ann_ids, open(os.path.join(model_dir, 'reverse_ann_ids.pkl') ,'wb'))

    if use_features:
        # For Jari's features
        pickle.dump(json_id_map, open(os.path.join(model_dir, 'json_id_map.pkl') ,'wb'))
        pickle.dump(json_vectorizer, open(os.path.join(model_dir, 'json_vectorizer.pkl') ,'wb'))
        pickle.dump(feature_selector, open(os.path.join(model_dir, 'feature_selector.pkl') ,'wb'))

    es_cb = EarlyStopping(monitor='val_acc', patience=10, verbose=0, mode='max')
    cp_cb = ModelCheckpoint(filepath=os.path.join(model_dir, 'model.hdf5'), monitor='val_acc', mode='max', save_best_only=True,verbose=0)
    ev_cb = Evaluate(devel_path, 500, reverse_ann_ids)
    # next(devel_data)
    # import pdb; pdb.set_trace()
    model.fit_generator(train_data, steps_per_epoch=batch_size, nb_epoch=60, validation_data=devel_data, validation_steps=batch_size, callbacks=[ev_cb])

        # If using our own blast features
        #pickle.dump(blast_hit_ids, open(os.path.join(model_dir, 'blast_hit_ids.pkl') ,'wb'))
        

    #import pdb; pdb.set_trace()
    
#     print "Making predictions"
#     from keras.models import load_model
#     model = load_model(filepath=os.path.join(model_dir, 'model.h5'), custom_objects={"weighted_binary_crossentropy":weighted_binary_crossentropy})
#     
#     # # Reinstantiate the data generators, otherwise they are not correctly aligned anymore
#     # devel_data = generate_data(devel_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size)
#     # test_data = generate_data(test_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size)
#     # 
#     # devel_score = model.evaluate_generator(devel_data, devel_size)
#     # test_score = model.evaluate_generator(test_data, test_size)
#     # print 'Devel l/a/p/r/f: ', devel_score
#     # print 'Test l/a/p/r/f: ', test_score
# 
# 
#     # Reinstantiate the data generators, otherwise they are not correctly aligned anymore
#     devel_data = generate_data(devel_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size)
#     test_data = generate_data(test_path, SEQUENCE_PATH, ann_path, ann_ids, batch_size)
#     
#     devel_pred = model.predict_generator(devel_data, steps=batch_size)
#     test_pred = model.predict_generator(test_data, steps=batch_size)
#     
#     save_predictions(os.path.join(model_dir, 'devel_pred.tsv.gz'), devel_ids, devel_pred, reverse_ann_ids)
#     save_predictions(os.path.join(model_dir, 'test_pred.tsv.gz'), test_ids, test_pred, reverse_ann_ids)
#     
#     print 'Making CAFA target predictions'
#     
#     cafa_id_path = '/home/sukaew/CAFA_PI/targetFiles/sequences/target.all.ids.gz'
#     cafa_seq_path = '/home/sukaew/CAFA_PI/targetFiles/sequences/target.all.fasta.gz'
#     cafa_data = generate_data(None, cafa_seq_path, ann_path, ann_ids, batch_size, cafa_targets=True, verbose=False)
#     cafa_size = _data_size(cafa_id_path)
#     cafa_ids = _get_ids(generate_data(None, cafa_seq_path, ann_path, ann_ids, batch_size, cafa_targets=True, verbose=False, endless=False))
# #    cafa_ids = read_split_ids(cafa_id_path, unique=False)
#     #import pdb; pdb.set_trace()
#     cafa_pred = model.predict_generator(cafa_data, batch_size)
#     
#     save_predictions(os.path.join(model_dir, 'cafa_targets.tsv.gz'), cafa_ids, cafa_pred, reverse_ann_ids, cafa_targets=True)
    # #import pdb; pdb.set_trace()
    
    print('All done.')
    

def save_predictions(out_path, prot_ids, predictions, reverse_ann_ids, cafa_targets=False):
    """
    Saves predictions in tsv format
    """
    import csv
    with gzip.open(out_path, 'wt') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['id', 'label_index', 'label', 'predicted', 'confidence', 'cafa_ids']) # Header line
        
        for i, prot_id in enumerate(prot_ids):
            if i % 10000 == 0:
                print(i)
            pred_indices = np.round(predictions[i]).nonzero()[0]
            if len(pred_indices) > 1500:
                print('WARNING: Maximum GO amount exceeded!')
                import pdb; pdb.set_trace()
            for pred_i in pred_indices:
                go_id = reverse_ann_ids[pred_i]
                confidence = predictions[i, pred_i]
                if cafa_targets:
                    cafa_id, p_id = prot_id.split(' ')
                    writer.writerow([p_id, pred_i, go_id, 1, '%.2f' % confidence, cafa_id])
                else:
                    writer.writerow([prot_id, pred_i, go_id, 1, '%.2f' % confidence, ''])

    
def weighted_binary_crossentropy(target, output):
    from keras.backend.common import _EPSILON
    from keras import backend as K
    # from theano import tensor as T
    # from theano.tensor import basic as tensor
    pos_weight = 1 # 73 is roughly the inverse ratio of positives examples
    output = K.clip(output, _EPSILON, 1.0 - _EPSILON)
    ce = -(pos_weight * target * K.log(output) + (1.0 - target) * K.log(1.0 - output))
    return K.mean(ce, axis=-1)
    
if __name__ == '__main__':
    train()
