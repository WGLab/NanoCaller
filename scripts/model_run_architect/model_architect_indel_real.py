import tensorflow as tf
import numpy as np
        
def conv2d(x, W, b, strides=1,pad_type='VALID'):
    x = tf.nn.conv2d(x, W, strides=strides, padding=pad_type)
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)




def conv_net(x,weights, biases,rate):
    conv1_1 = conv2d(x, weights['wc1_1'], biases['bc1_1'],strides=[1,1,1,1],pad_type='SAME')
    conv1_2 = conv2d(x, weights['wc1_2'], biases['bc1_2'],strides=[1,1,1,1],pad_type='SAME')
    conv1_3 = conv2d(x, weights['wc1_3'], biases['bc1_3'],strides=[1,1,1,1],pad_type='SAME')
    
    merge_conv1=tf.concat([conv1_1, conv1_2,conv1_3],3)
    conv2 = conv2d(merge_conv1, weights['wc2'], biases['bc2'],strides=[1,1,2,1])
    
    #pooled=tf.nn.max_pool(merge_conv1,ksize=[1,5,1,1],strides=[1,1,1,1],padding='SAME')
    #conv2 = conv2d(pooled, weights['wc2'], biases['bc2'],strides=[1,1,2,1])
    
    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],strides=[1,1,2,1])

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])

    fc1 = tf.add(tf.matmul(flat_nn, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    drop_out_fc1=tf.nn.dropout(fc1, rate=rate)
    
    '''fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd2']), biases['bd2'])
    fc2 = tf.nn.selu(fc2)
    drop_out_fc2=tf.nn.dropout(fc2, rate)'''
    
    fc3 = tf.add(tf.matmul(drop_out_fc1, weights['wd3']), biases['bd3'])
    fc3 = tf.nn.selu(fc3)
    
    out = tf.add(tf.matmul(fc3, weights['out']), biases['out'])
    
    return out

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    #h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    #w_in=int(np.ceil(float(w_in-4+1)/float(1)))
    #h_in=28
       
    h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    w_in=int(np.ceil(float(w_in-3+1)/float(2)))

    h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    w_in=int(np.ceil(float(w_in-3+1)/float(2)))
    
    a,b,c=8,8,8
    d,e,f=32,48,32
    q=32
    g,h=24,16
    j=8
    i=0

    weights = {
        'wc1_1': tf.compat.v1.get_variable('W0_1', shape=(1,5,depth,a), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc1_2': tf.compat.v1.get_variable('W0_2', shape=(5,1,depth,b), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.compat.v1.get_variable('W0_3', shape=(5,5,depth,c), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.compat.v1.get_variable('W1', shape=(2,3,a+b+c,d), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.compat.v1.get_variable('W2', shape=(2,3,d,e), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.compat.v1.get_variable('W3', shape=(h_in*w_in*e,f), initializer=tf.contrib.layers.xavier_initializer()),
        #'wd2': tf.get_variable('W3_2', shape=(f,q), initializer=tf.contrib.layers.xavier_initializer()),
        'wd3': tf.compat.v1.get_variable('W3_3', shape=(f,g), initializer=tf.contrib.layers.xavier_initializer()),
      
        'out': tf.compat.v1.get_variable('W13', shape=(g,3), initializer=tf.contrib.layers.xavier_initializer())
    }
    
    biases = {
        'bc1_1': tf.compat.v1.get_variable('B0_1', shape=(a), initializer=tf.contrib.layers.xavier_initializer()), 
        'bc1_2': tf.compat.v1.get_variable('B0_2', shape=(b), initializer=tf.contrib.layers.xavier_initializer()),
        'bc1_3': tf.compat.v1.get_variable('B0_3', shape=(c), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.compat.v1.get_variable('B1', shape=(d), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.compat.v1.get_variable('B2', shape=(e), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.compat.v1.get_variable('B3', shape=(f), initializer=tf.contrib.layers.xavier_initializer()),
        #'bd2': tf.get_variable('B3_2', shape=(q), initializer=tf.contrib.layers.xavier_initializer()),
        'bd3': tf.compat.v1.get_variable('B3_3', shape=(g), initializer=tf.contrib.layers.xavier_initializer()),
        'out': tf.compat.v1.get_variable('B13', shape=(3), initializer=tf.contrib.layers.xavier_initializer()),
    }
    
    
    x0,x1 = tf.compat.v1.placeholder("float", [None]+n_input), tf.compat.v1.placeholder("float", [None]+n_input),
    
    label_1=tf.compat.v1.placeholder("float", [None, 3])
    label_2=tf.compat.v1.placeholder("float", [None, 3])
    
    
    rate = tf.compat.v1.placeholder(tf.float32)

    score_0 = conv_net(x0, weights, biases, rate)
    score_1 = conv_net(x1, weights, biases, rate)

    gamma=0.0001
    reg_loss=tf.add_n([tf.nn.l2_loss(t) for t in weights.values()])
    
    prob_0=tf.nn.softmax(logits=score_0)
    prob_1=tf.nn.softmax(logits=score_1)

    prediction_0=tf.argmax(prob_0, 1)
    prediction_1=tf.argmax(prob_1, 1)

    
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=score_0, labels=label_1))+tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=score_1, labels=label_2))+ tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    
    optimizer = tf.compat.v1.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    

    
    prediction_0_accuracy= tf.equal(prediction_0, tf.argmax(label_1, 1))
    prediction_1_accuracy= tf.equal(prediction_1, tf.argmax(label_2, 1))
    
  
    
    accuracy_0 = tf.reduce_sum(tf.cast(prediction_0_accuracy, tf.float32))
    accuracy_1 = tf.reduce_sum(tf.cast(prediction_1_accuracy, tf.float32))
   
    accuracy=tf.reduce_sum(tf.cast(tf.logical_and(prediction_0_accuracy,prediction_1_accuracy), tf.float32))
    
    

    tensors=(x0, x1, label_1,label_2, accuracy, accuracy_0, accuracy_1, cost, optimizer, prob_0,prob_1, rate)

    
    
    return weights,biases,tensors
