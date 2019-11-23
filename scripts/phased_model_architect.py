import tensorflow as tf
import numpy as np

def conv2d(x, W, b, strides=[1, 1, 1, 1],pad_type='VALID'):
    x = tf.nn.conv2d(x, W, strides=strides, padding=pad_type)
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)

def conv_net(x,ref, weights, biases,keep):
    conv1_1 = conv2d(x, weights['wc1_1'], biases['bc1_1'],pad_type='SAME')
    conv1_2 = conv2d(x, weights['wc1_2'], biases['bc1_2'],pad_type='SAME')
    conv1_3 = conv2d(x, weights['wc1_3'], biases['bc1_3'],pad_type='SAME')
    
    merge_conv1=tf.concat([conv1_1, conv1_2,conv1_3],3)
    
    pool_output=tf.nn.max_pool(merge_conv1,ksize=[1,5,1,1],strides=[1,1,1,1],padding='SAME')
    conv2 = conv2d(pool_output, weights['wc2'], biases['bc2'],strides=[1,2,1,1])
    
    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],strides=[1,2,1,1])

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    
    drop_out_flat =flat_nn# tf.nn.dropout(flat_nn, keep)
    
    fc1 = tf.add(tf.matmul(drop_out_flat, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    
    drop_out_fc1=tf.nn.dropout(fc1, keep)
   
    fa = tf.add(tf.matmul(tf.concat([drop_out_fc1],1), weights['w-tr']), biases['w-tr'])
    fa = tf.nn.selu(fa)
    out_allele = tf.add(tf.matmul(fa, weights['out-all']), biases['out-all'])

    return out_allele

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    
    for i in range(2):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(1)))
       
    a,b,c=8,8,8
    d,e,f=32,48,32
    g,h=16,16
    i=0
    '''
    a,b,c=8,8,8
    d,e,f=32,48,32
    g,h=16,16
    i=0
    '''
    weights = {
        'wc1_1': tf.get_variable('W0_1', shape=(1,5,depth,a), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc1_2': tf.get_variable('W0_2', shape=(8,1,depth,b), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.get_variable('W0_3', shape=(3,3,depth,c), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.get_variable('W1', shape=(3,3,a+b+c,d), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,d,e), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*e,f), initializer=tf.contrib.layers.xavier_initializer()), 
        'w-tr': tf.get_variable('W4', shape=(i+f,g), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('W5', shape=(g,4), initializer=tf.contrib.layers.xavier_initializer())
    }
    
    biases = {
        'bc1_1': tf.get_variable('B0_1', shape=(a), initializer=tf.contrib.layers.xavier_initializer()), 
        'bc1_2': tf.get_variable('B0_2', shape=(b), initializer=tf.contrib.layers.xavier_initializer()),
        'bc1_3': tf.get_variable('B0_3', shape=(c), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(d), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(e), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(f), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('B4', shape=(g), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('B5', shape=(4), initializer=tf.contrib.layers.xavier_initializer())
    }
    
    x = tf.placeholder("float", [None,None,n_input[-2],n_input[-1]])
    allele,ref=tf.placeholder("float", [None, 4]),tf.placeholder("float", [None, 4])
    keep = tf.placeholder(tf.float32)
    
    fc_layer = conv_net(x,ref,weights, biases, keep)
    
    gamma=0.001
    trade_off=1.0
    reg_loss=tf.add_n([tf.nn.l2_loss(t) for t in weights.values()])
                      
    cost=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=fc_layer, labels=allele)) + tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    
    correct_prediction= tf.cast(tf.equal(tf.argmax(fc_layer, 1), tf.argmax(allele, 1)), tf.float32)
    accuracy = tf.reduce_sum(correct_prediction)
        
    tensors=(x,allele,ref,fc_layer,cost,optimizer,keep,correct_prediction, accuracy)
    return weights,biases,tensors