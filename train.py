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

ann_limit = 5000 # Taking top N GO annotations only
timesteps = 2500 # maximum length of a sequence, the real max is 35K. 2.5K covers 99% of the sequences, 5K 99.9%
latent_dim = 50 # Amino acid embedding size
batch_size = 16
char_set = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'
vocab_size = len(char_set) + 1 # +1 for mask
char_dict = {c:i+1 for i,c in enumerate(char_set)} # Index 0 is left for padding
use_features = False # False = only sequence is used for prediction
model_dir = './test_model/' # path for saving model + other required stuff
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

def read_split_ids(split_path, unique=True):
    split_file = gzip.open(split_path)
    split_data = [s.strip().replace('>', '') for s in split_file]
    if unique:
        split_data = set(split_data)
    return split_data

def generate_data(split_path, seq_path, ann_path, ann_ids, batch_size=256, cafa_targets=False, verbose=False, endless=True):
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
            if verbose and i % 10000 == 0:
                print 'Generator: ', i
            if cafa_targets:
                cafa_id, prot_id = prot_id.strip().replace('>', '').split(' ')
            else:
                prot_id = prot_id.strip().strip('>')
            if split_path and prot_id not in split_ids:
                continue
            if cafa_targets:
                prot_ids.append('%s %s' % (cafa_id, prot_id))
            else:
                prot_ids.append(prot_id)
            seq = seq.strip()
            seq_id_list = [char_dict[s] for s in seq if s != '*']
            #seq_id_list = [aa_index_ids.get(s, 0) for s in seq]
            annotations = ann_dict[prot_id]
            ann_id_list = [ann_ids[a] for a in annotations if a in ann_ids]
            y_v = np.zeros((len(ann_ids)), dtype='int')
            y_v[ann_id_list] = 1
            x.append(seq_id_list)
            if use_features:
                #blast_x.append(generate_blast_features(prot_id)) # These are our original blast features
                blast_x.append(get_feature_vector(prot_id)) # Jari's feature vectors
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
            
        if not endless:
            break

def read_aaindex():
    aa_f = open('/home/sukaew/CAFA3/aaindex/aaindex_table.tsv')
    data = aa_f.readlines()
    amino_acids = data[0].strip().split('\t')[1:]
    aa_index_ids = {a: i+1 for i, a in enumerate(amino_acids)}
    
    embedding = np.zeros((len(amino_acids)+1, len(data)-1)) # +1 for OOV/MASK no separation for now
    for i, row in enumerate(data[1:]):
        values = row.strip().split('\t')[1:]
        for ii, v in enumerate(values):
            if v == 'NA': # FIXME: What to do with these?
                v = 0.0
            embedding[ii+1, i] = float(v)
        
    #import pdb; pdb.set_trace()
    return aa_index_ids, embedding

aa_index_ids, aa_embedding = read_aaindex()

def read_feature_json(path='./data/examples.json.gz', vectorizer=None, feature_selector=None):
    print 'Reading Jaris feature data'
    js = json.load(gzip.open(path))
    
    filters = ['DUMMY']#['BLAST', 'DELTA', 'GPI', 'TAX', 'IPS']
    print "Excluding: ", filters
    for i, d in enumerate(js['features']):
        if i % 10000 == 0:
            print i
        for f in filters:
            for key in d.keys():
                if key.startswith('%s:' % f):
                    d.pop(key)
    
    
    from sklearn.feature_extraction import DictVectorizer
    if vectorizer:
        v = vectorizer
    else:
        v = DictVectorizer()

    feature_matrix = v.fit_transform(js['features'])
    
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
    return len(gzip.open(path).readlines())

def _get_ids(data_gen):
    """
    Needs a data gen with endless=False
    """
    ids = [i[0]['prot_ids'] for i in data_gen]
    ids = list(np.concatenate(ids))
    return ids

def train():
    print 'Generating training data'
    ann_path = './data/Swissprot_propagated.tsv.gz'
    ann_ids, reverse_ann_ids = get_annotation_ids(ann_path, top=ann_limit)

    #pretrain_data = generate_data(None, '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz', '/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ann-train-data.tsv.gz', ann_ids, 256)
    #pretrain_size = _data_size('/home/hanmoe/CAFA3/ngrams/4kai/assocI-min_len5-min_freq3-top_fun5k/ngram-id2seq.tsv.gz')/2
    train_path = './data/train.txt.gz'
    train_data = generate_data(train_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    train_size = _data_size(train_path)
    train_ids = read_split_ids(train_path, unique=False)
    
    devel_path = './data/devel.txt.gz'
    devel_data = generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    devel_size = _data_size(devel_path)
    devel_ids = read_split_ids(devel_path, unique=False)
    
    test_path = './data/test.txt.gz'
    test_data = generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
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
    print 'Building model'
    inputs = Input(shape=(timesteps, ), name='sequence')
    input_list = [inputs]
    embedding = Embedding(vocab_size, latent_dim, mask_zero=False)(inputs)
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
        #feature_input = Input(shape=(len(blast_hit_ids), ), name='features')
        feature_input = Input(shape=(json_feature_matrix.shape[1], ), name='features') # For Jari's feature vectors
        dropout = Dropout(0.5)(feature_input)
        feature_encoding = Dense(300, activation='tanh')(dropout) # Squeeze the feature vectors to a tiny encoding
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
    es_cb = EarlyStopping(monitor='val_fmeasure', patience=10, verbose=0, mode='max')
    cp_cb = ModelCheckpoint(filepath=os.path.join(model_dir, 'model.hdf5'), monitor='val_fmeasure', mode='max', save_best_only=True,verbose=0)
    #model.fit(pretrain_data, pretrain_data, nb_epoch=100, batch_size=16, validation_data=[devel_data, devel_data], callbacks=[es_cb, cp_cb])
    model.fit_generator(train_data, samples_per_epoch=train_size, nb_epoch=100, validation_data=devel_data, nb_val_samples=devel_size, callbacks=[es_cb, cp_cb])
    
    pickle.dump(ann_ids, open(os.path.join(model_dir, 'ann_ids.pkl') ,'wb'))
    pickle.dump(reverse_ann_ids, open(os.path.join(model_dir, 'reverse_ann_ids.pkl') ,'wb'))
    if use_features:
        # For Jari's features
        pickle.dump(json_id_map, open(os.path.join(model_dir, 'json_id_map.pkl') ,'wb'))
        pickle.dump(json_vectorizer, open(os.path.join(model_dir, 'json_vectorizer.pkl') ,'wb'))
        pickle.dump(feature_selector, open(os.path.join(model_dir, 'feature_selector.pkl') ,'wb'))

        # If using our own blast features
        #pickle.dump(blast_hit_ids, open(os.path.join(model_dir, 'blast_hit_ids.pkl') ,'wb'))
        

    #import pdb; pdb.set_trace()
    
    print "Making predictions"
    from keras.models import load_model
    model = load_model(filepath=os.path.join(model_dir, 'model.hdf5'), custom_objects={"weighted_binary_crossentropy":weighted_binary_crossentropy})

    # Reinstantiate the data generators, otherwise they are not correctly aligned anymore
    devel_data = generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    test_data = generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)

    devel_score = model.evaluate_generator(devel_data, devel_size)
    test_score = model.evaluate_generator(test_data, test_size)
    print 'Devel l/a/p/r/f: ', devel_score
    print 'Test l/a/p/r/f: ', test_score


    # Reinstantiate the data generators, otherwise they are not correctly aligned anymore
    devel_data = generate_data(devel_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)
    test_data = generate_data(test_path, './data/Swissprot_sequence.tsv.gz', ann_path, ann_ids, batch_size)

    devel_pred = model.predict_generator(devel_data, devel_size)
    test_pred = model.predict_generator(test_data, test_size)
    
    save_predictions(os.path.join(model_dir, 'devel_pred.tsv.gz'), devel_ids, devel_pred, reverse_ann_ids)
    save_predictions(os.path.join(model_dir, 'test_pred.tsv.gz'), test_ids, test_pred, reverse_ann_ids)
    
    print 'Making CAFA target predictions'
    
    cafa_id_path = './data/target.all.ids.gz'
    cafa_seq_path = '/home/sukaew/CAFA3/CAFA3_targets/target.all.fasta.gz'
    cafa_data = generate_data(None, cafa_seq_path, ann_path, ann_ids, batch_size=batch_size, cafa_targets=True, verbose=True)
    cafa_size = _data_size(cafa_id_path)
    cafa_ids = read_split_ids(cafa_id_path, unique=False)
    #import pdb; pdb.set_trace()
    cafa_pred = model.predict_generator(cafa_data, cafa_size)
    
    save_predictions(os.path.join(model_dir, 'cafa_targets.tsv.gz'), cafa_ids, cafa_pred, reverse_ann_ids, cafa_targets=True)
    #import pdb; pdb.set_trace()
    
    print 'All done.'
    

def save_predictions(out_path, prot_ids, predictions, reverse_ann_ids, cafa_targets=False):
    """
    Saves predictions in tsv format
    """
    import csv
    with gzip.open(out_path, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(['id', 'label_index', 'label', 'predicted', 'confidence', 'cafa_ids']) # Header line
        
        for i, prot_id in enumerate(prot_ids):
            if i % 10000 == 0:
                print i
            pred_indices = np.round(predictions[i]).nonzero()[0]
            if len(pred_indices) > 1500:
                print 'WARNING: Maximum GO amount exceeded!'
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
    from theano import tensor as T
    from theano.tensor import basic as tensor
    pos_weight = 3.0 # 73 is roughly the inverse ratio of positives examples
    output = T.clip(output, _EPSILON, 1.0 - _EPSILON)
    ce = -(pos_weight * target * tensor.log(output) + (1.0 - target) * tensor.log(1.0 - output))
    return K.mean(ce, axis=-1)
    
if __name__ == '__main__':
    train()