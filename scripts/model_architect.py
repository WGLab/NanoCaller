import tensorflow as tf
import numpy as np

def conv2d(x, W, b, strides=1,pad_type='VALID'):
    x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding=pad_type)
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)




def conv_net(x,ref, weights, biases,keep):
    conv1_1 = conv2d(x, weights['wc1_1'], biases['bc1_1'],1,pad_type='SAME')
    conv1_2 = conv2d(x, weights['wc1_2'], biases['bc1_2'],1,pad_type='SAME')
    conv1_3 = conv2d(x, weights['wc1_3'], biases['bc1_3'],1,pad_type='SAME')
    
    merge_conv1=tf.concat([conv1_1, conv1_2,conv1_3],3)
    print(tf.shape(merge_conv1))
    pool_output = tf.cond(tf.equal(tf.shape(merge_conv1)[1], 100),
                      lambda:tf.nn.max_pool(merge_conv1,ksize=[1,18,1,1],strides=[1,3,1,1],padding='VALID'),
                      lambda: tf.nn.max_pool(merge_conv1,ksize=[1,5,1,1],strides=[1,1,1,1],padding='VALID'))
    
    
    conv2 = conv2d(pool_output, weights['wc2'], biases['bc2'],2)
    
    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],2)

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    
    drop_out_flat =flat_nn# tf.nn.dropout(flat_nn, keep)
    
    fc1 = tf.add(tf.matmul(drop_out_flat, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    
    drop_out_fc1=tf.nn.dropout(fc1, keep)
    
    fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd-gtp']), biases['bd-gtp'])
    fc2 = tf.nn.selu(fc2)
    out_gt= tf.add(tf.matmul(fc2, weights['out-gtp']), biases['out-gtp'])
    
    
    fa = tf.add(tf.matmul(tf.concat([drop_out_fc1],1), weights['w-tr']), biases['w-tr'])
    fa = tf.nn.selu(fa)
    out_allele = tf.add(tf.matmul(fa, weights['out-all']), biases['out-all'])

    return out_gt,out_allele

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    
    for i in range(2):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(2)))
       
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
        'wc1_2': tf.get_variable('W0_2', shape=(33,1,depth,b), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.get_variable('W0_3', shape=(3,3,depth,c), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.get_variable('W1', shape=(3,3,a+b+c,d), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,d,e), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*e,f), initializer=tf.contrib.layers.xavier_initializer()), 
        'w-tr': tf.get_variable('W4', shape=(i+f,g), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('W5', shape=(g,4), initializer=tf.contrib.layers.xavier_initializer()),
        'wd-gtp': tf.get_variable('W6', shape=(f,h), initializer=tf.contrib.layers.xavier_initializer()), 
        'out-gtp': tf.get_variable('W7', shape=(h,2), initializer=tf.contrib.layers.xavier_initializer()), 
    }
    
    biases = {
        'bc1_1': tf.get_variable('B0_1', shape=(a), initializer=tf.contrib.layers.xavier_initializer()), 
        'bc1_2': tf.get_variable('B0_2', shape=(b), initializer=tf.contrib.layers.xavier_initializer()),
        'bc1_3': tf.get_variable('B0_3', shape=(c), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(d), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(e), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(f), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('B4', shape=(g), initializer=tf.contrib.layers.xavier_initializer()),
        'out-all': tf.get_variable('B5', shape=(4), initializer=tf.contrib.layers.xavier_initializer()),
        'bd-gtp': tf.get_variable('B6', shape=(h), initializer=tf.contrib.layers.xavier_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
    }
    
    x,y = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, 2])
    allele,ref=tf.placeholder("float", [None, 4]),tf.placeholder("float", [None, 4])
    keep = tf.placeholder(tf.float32)
    
    pred,fc_layer = conv_net(x,ref,weights, biases, keep)
    
    gt_likelihood=tf.nn.softmax(logits=pred)
    allele_likelihood=tf.nn.softmax(logits=fc_layer)
    
    cost_gt=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=pred,\
    labels=y))
    cost_allele=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=fc_layer, labels=allele))
    gamma=0.001
    trade_off=1.0
    reg_loss=tf.add_n([tf.nn.l2_loss(t) for t in weights.values()])
    #+[tf.constant(0.005)*tf.nn.l2_loss(weights['w-tr'])+tf.constant(0.005)*tf.nn.l2_loss(weights['out-all'])])
                      
    cost = cost_allele+cost_gt+tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    
    correct_prediction_gt = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    accuracy_gt = tf.reduce_sum(tf.cast(correct_prediction_gt, tf.float32))
    
    correct_prediction_allele = tf.equal(tf.argmax(fc_layer, 1), tf.argmax(allele, 1))
    correct_prediction=tf.cast(correct_prediction_allele, tf.float32)
    accuracy_allele = tf.reduce_sum(correct_prediction)
    
    accuracy=tf.reduce_sum(tf.cast(tf.logical_and(correct_prediction_allele,correct_prediction_gt), tf.float32))
    
    t1=(x,y,allele,ref,fc_layer,pred,cost,optimizer,cost_gt,cost_allele,keep)
    t2=(correct_prediction, correct_prediction_gt, correct_prediction_allele, accuracy, accuracy_gt, accuracy_allele, gt_likelihood, allele_likelihood)
    return weights,biases,t1,t2