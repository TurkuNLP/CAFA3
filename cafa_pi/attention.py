from __future__ import print_function
import numpy
import numpy as np
import random

from keras import backend as K
from keras.engine.topology import Layer

import theano.tensor as T
import theano

class Attention(Layer):

    def __init__(self, return_sequence=False, **kwargs):

        self.return_sequences = return_sequence
        self.supports_masking = True
        super(Attention, self).__init__(**kwargs)

    def build(self, input_shape):
        input_dim = input_shape[0][-1]
        self.x_input_shape = input_shape
        self.W_h = K.variable(np.random.random((input_dim, input_dim)), name='W_h')
        self.W_y = K.variable(np.random.random((input_dim, input_dim)), name='W_y')
        self.w_a = K.variable(np.random.random((input_dim,)), name='w_a')

        self.trainable_weights = [self.W_y, self.w_a, self.W_h]

    def call(self, x, mask=None):

        full_input = x[0]
        encoding = x[1]

        whhn = T.dot(encoding, self.W_h.T) #Maybe T here as well?
        whhn_r = T.extra_ops.repeat(whhn, full_input.shape[1], axis=0).reshape(full_input.shape)

        M = T.tanh(T.dot(full_input, self.W_y.T) + whhn_r)
        a_dot = T.nnet.nnet.sigmoid((T.dot(M, self.w_a)))

        if mask is None:
            a = T.nnet.softmax(a_dot)
        else:
            mask = mask[0]
            #This is softmax with a mask, maybe usable somewhere else as well?
            if mask == None:
                mask = 1
            softmax_sum_with_mask = T.sum(T.exp(a_dot) * mask, axis=1)
            repeated_sum = T.extra_ops.repeat(softmax_sum_with_mask, a_dot.shape[1]).reshape(a_dot.shape)
            a = ((T.exp(a_dot)/repeated_sum)) * mask

        ar = T.extra_ops.repeat(a, full_input.shape[-1], axis=1).reshape(full_input.shape)

        if self.return_sequences:
            return full_input * ar
        else:
            return T.sum(full_input * ar, axis=1)

    def call_softmax(self, x, mask=None):

        full_input = x[0]
        encoding = x[1]

        whhn = T.dot(encoding, self.W_h.T)
        whhn_r = T.extra_ops.repeat(whhn, full_input.shape[1], axis=0).reshape(full_input.shape)

        M = T.tanh(T.dot(full_input, self.W_y.T) + whhn_r)
        a_dot = T.dot(M, self.w_a)
        mask = mask[0]

        if mask is None:
            a = T.nnet.softmax(a_dot)
        else:
            #This is softmax with a mask, maybe usable somewhere else as well?
            softmax_sum_with_mask = T.sum(T.exp(a_dot) * mask, axis=1)
            repeated_sum = T.extra_ops.repeat(softmax_sum_with_mask, a_dot.shape[1]).reshape(a_dot.shape)
            a = ((T.exp(a_dot)/repeated_sum)) * mask

        return a

    def get_output_shape_for(self, input_shape):

        if self.return_sequences:
            return self.x_input_shape[0]
        else:
            return self.x_input_shape[1]

    def compute_mask(self, input, mask):
        if self.return_sequences:
            return mask[0]
        else:
            return None

    def get_config(self):
        base_config = super(Attention, self).get_config()
        #config = {'output_dim' : self.output_dim}
        return base_config

def main():

    window = 5
    vec_size = 20
    vec_len = 20

    x = []
    y = []
    for ax in range(90000):
        #
        len_r = random.randint(2, window - 1)
        cx = numpy.concatenate([numpy.ones(len_r), numpy.zeros(window-len_r)])
        #label
        label = random.randint(0,1)
        if label > 0:
            cx[random.randint(0,len_r - 1)] = 2
        else:
            cx[random.randint(0,len_r - 1)] = 3          
        x.append(cx)
        y.append(label)

    x = numpy.array(x, dtype=np.int32)
    y = numpy.array(y)

    dev_x = x[8000:]
    dev_y = y[9000:]

    x = x[:8000]
    y = y[:8000]

    from keras.layers import Dense, Dropout, Activation, Merge, Input, merge
    from keras.layers.embeddings import Embedding
    from keras.models import Model
    from keras.layers.core import Dense, Activation, Dropout, RepeatVector, Merge, TimeDistributedDense, Flatten
    from keras.layers.recurrent import GRU

    #Input
    x_input = Input(shape=(window, ), name='x_input', dtype='int32')

    #Embeddings
    char_emb = Embedding(4, vec_size, input_length=window, mask_zero=True)
    emb_out = char_emb(x_input)

    #Attention
    simple_attn = Attention()

    #Let's make an encoding
    rnn = GRU(vec_size)
    encoding = rnn(emb_out)

    #Attention takes as its input the input sequence and an encoding of it
    emb_attn = simple_attn([emb_out, encoding])

    #Output dense
    l_dense = Dense(1, activation='sigmoid')
    d_out = l_dense(emb_attn)
    model = Model(input=[x_input], output=d_out)
    model.compile(optimizer='adam', loss='mse', metrics=['accuracy'])

    model.fit(x,y, nb_epoch=2)


    t_input = T.imatrix()
    testf = theano.function([t_input], simple_attn.call_softmax( [char_emb.call(t_input), rnn.call(char_emb.call(t_input)) ],  [char_emb.compute_mask(t_input)] ) )

    np.set_printoptions(precision=3)
    print()
    print('Testing Attention')
    print()
    for example in dev_x[:5]:
        print('Example', example)
        print ('Attention')
        print (testf(np.array([example],dtype=np.int32)), np.sum(testf(np.array([example],dtype=np.int32)), axis=1) )
        print ()
        print ()

    from keras.models import model_from_json
    print ('Testing JSON load & save')
    json_mdl = model.to_json()
    print (json_mdl)
    model_clone = model_from_json(json_mdl, custom_objects={"Attention": Attention})
    model_clone.compile(optimizer='adam', loss='mse', metrics=['accuracy'])
    model_clone.fit(x, y, nb_epoch=2)

if __name__ == "__main__": main()