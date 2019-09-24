import tensorflow as tf
import numpy as np
        
def conv2d(x, W, b, strides=1,pad_type='VALID'):
    x = tf.nn.conv2d(x, W, strides=strides, padding=pad_type)
    x = tf.nn.bias_add(x, b)
    return tf.nn.selu(x)




def conv_net(x,A_ref,G_ref,T_ref,C_ref,weights, biases,keep):
    conv1_1 = conv2d(x, weights['wc1_1'], biases['bc1_1'],strides=[1,1,1,1],pad_type='SAME')
    conv1_2 = conv2d(x, weights['wc1_2'], biases['bc1_2'],strides=[1,1,1,1],pad_type='SAME')
    conv1_3 = conv2d(x, weights['wc1_3'], biases['bc1_3'],strides=[1,1,1,1],pad_type='SAME')
    
    merge_conv1=tf.concat([conv1_1, conv1_2,conv1_3],3)
    conv2 = conv2d(merge_conv1, weights['wc2'], biases['bc2'],strides=[1,1,2,1])
    
    #pooled=tf.nn.max_pool(merge_conv1,ksize=[1,5,1,1],strides=[1,1,1,1],padding='SAME')
    #conv2 = conv2d(pooled, weights['wc2'], biases['bc2'],strides=[1,1,2,1])
    
    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],strides=[1,1,2,1])

    flat_nn = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    
    drop_out_flat =flat_nn# tf.nn.dropout(flat_nn, keep)
    
    fc1 = tf.add(tf.matmul(drop_out_flat, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.selu(fc1)
    
    drop_out_fc1=tf.nn.dropout(fc1, keep)
    
    '''fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd2']), biases['bd2'])
    fc2 = tf.nn.selu(fc2)
    
    drop_out_fc2=tf.nn.dropout(fc2, keep)'''
    
    fa = tf.add(tf.matmul(drop_out_fc1, weights['w-tr']), biases['w-tr'])
    fa = tf.nn.selu(fa)
    
    out_A = tf.add(tf.matmul(tf.concat([fa,A_ref],1), weights['out-A']), biases['out-A'])
    out_G = tf.add(tf.matmul(tf.concat([fa,G_ref],1), weights['out-G']), biases['out-G'])
    out_T = tf.add(tf.matmul(tf.concat([fa,T_ref],1), weights['out-T']), biases['out-T'])
    out_C = tf.add(tf.matmul(tf.concat([fa,C_ref],1), weights['out-C']), biases['out-C'])
    
    #fc2 = tf.add(tf.matmul(tf.concat([drop_out_fc1, out_A, out_T, out_G, out_C],1), weights['wd-gtp_0']), biases['bd-gtp_0'])
    
    fc2 = tf.add(tf.matmul(drop_out_fc1, weights['wd-gtp_0']), biases['bd-gtp_0'])
    gt_fc1 = tf.nn.selu(fc2)
    
    
    gt_fc2=tf.nn.selu(tf.add(tf.matmul(gt_fc1, weights['wd-gtp']), biases['bd-gtp']))
    out_gt= tf.add(tf.matmul(gt_fc2, weights['out-gtp']), biases['out-gtp'])
    
    
    
    
    return out_gt,out_A,out_G,out_T,out_C

def get_tensors(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    
    #h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    #w_in=int(np.ceil(float(w_in-4+1)/float(1)))
    #h_in=28
       
    h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    w_in=int(np.ceil(float(w_in-3+1)/float(2)))

    h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    w_in=int(np.ceil(float(w_in-3+1)/float(2)))
    
    a,b,c=16,16,16
    d,e,f=32,64,48
    q=32
    g,h=16,16
    j=8
    i=0

    weights = {
        'wc1_1': tf.get_variable('W0_1', shape=(1,4,depth,a), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc1_2': tf.get_variable('W0_2', shape=(4,1,depth,b), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.get_variable('W0_3', shape=(3,3,depth,c), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.get_variable('W1', shape=(2,3,a+b+c,d), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(2,3,d,e), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*e,f), initializer=tf.contrib.layers.xavier_initializer()),
        #'wd2': tf.get_variable('W3_2', shape=(f,q), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('W4', shape=(f,g), initializer=tf.contrib.layers.xavier_initializer()),
        'wd-gtp_0': tf.get_variable('W5', shape=(f,h), initializer=tf.contrib.layers.xavier_initializer()),
        'wd-gtp': tf.get_variable('W6', shape=(h,j), initializer=tf.contrib.layers.xavier_initializer()), 
        'out-gtp': tf.get_variable('W7', shape=(j,2), initializer=tf.contrib.layers.xavier_initializer()),
        
        'out-A': tf.get_variable('W8', shape=(g+1,2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-G': tf.get_variable('W9', shape=(g+1,2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-T': tf.get_variable('W10', shape=(g+1,2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-C': tf.get_variable('W11', shape=(g+1,2), initializer=tf.contrib.layers.xavier_initializer())
    }
    
    biases = {
        'bc1_1': tf.get_variable('B0_1', shape=(a), initializer=tf.contrib.layers.xavier_initializer()), 
        'bc1_2': tf.get_variable('B0_2', shape=(b), initializer=tf.contrib.layers.xavier_initializer()),
        'bc1_3': tf.get_variable('B0_3', shape=(c), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(d), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(e), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(f), initializer=tf.contrib.layers.xavier_initializer()),
        #'bd2': tf.get_variable('B3_2', shape=(q), initializer=tf.contrib.layers.xavier_initializer()),
        'w-tr': tf.get_variable('B4', shape=(g), initializer=tf.contrib.layers.xavier_initializer()),
        'bd-gtp_0': tf.get_variable('B5', shape=(h), initializer=tf.contrib.layers.xavier_initializer()),
        'bd-gtp': tf.get_variable('B6', shape=(j), initializer=tf.contrib.layers.xavier_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
        
        'out-A': tf.get_variable('B8', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-G': tf.get_variable('B9', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-T': tf.get_variable('B10', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
        'out-C': tf.get_variable('B11', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
    }
    
    
    x,GT_label = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, 2])
    
    A_label,G_label,T_label,C_label=tf.placeholder("float", [None, 2]), tf.placeholder("float", [None, 2]), tf.placeholder("float", [None, 2]), tf.placeholder("float", [None, 2])
    
    A_ref,G_ref,T_ref,C_ref=tf.placeholder("float", [None, 1]), tf.placeholder("float", [None, 1]), tf.placeholder("float", [None, 1]), tf.placeholder("float", [None, 1])
    
    keep = tf.placeholder(tf.float32)
    
    GT_score,A_score,G_score,T_score,C_score = conv_net(x,A_ref,G_ref,T_ref,C_ref, weights, biases, keep)

    gamma=0.001
    reg_loss=tf.add_n([tf.nn.l2_loss(t) for t in weights.values()])
    
    
    
    
    
    
    prob_GT=tf.nn.softmax(logits=GT_score)
    prob_A=tf.nn.softmax(logits=A_score)
    prob_G=tf.nn.softmax(logits=G_score)
    prob_T=tf.nn.softmax(logits=T_score)
    prob_C=tf.nn.softmax(logits=C_score)
    
    cost_GT=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=GT_score, labels=GT_label))
    cost_A=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=A_score, labels=A_label))
    cost_G=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=G_score, labels=G_label))
    cost_T=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=T_score, labels=T_label))
    cost_C=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=C_score, labels=C_label))
    
    cost = cost_A+cost_G+cost_T+cost_C+tf.constant(1.0)*cost_GT+tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    prediction_GT=tf.argmax(prob_GT, 1)
    prediction_A=tf.argmax(prob_A, 1)
    prediction_G=tf.argmax(prob_G, 1)
    prediction_T=tf.argmax(prob_T, 1)
    prediction_C=tf.argmax(prob_C, 1)
    
    '''
    prob_GT=tf.nn.softmax(logits=GT_score)
    
    prob_A=tf.square(tf.subtract(A_label,A_score))
    prob_G=tf.square(tf.subtract(G_label,G_score))
    prob_T=tf.square(tf.subtract(T_label,T_score))
    prob_C=tf.square(tf.subtract(C_label,C_score))
    
    cost_GT=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=GT_score, labels=GT_label))
    cost_A=tf.reduce_mean(prob_A)
    cost_G=tf.reduce_mean(prob_G)
    cost_T=tf.reduce_mean(prob_T)
    cost_C=tf.reduce_mean(prob_C)
    
    cost = cost_A+cost_G+cost_T+cost_C+tf.constant(1.1)*cost_GT+tf.constant(gamma)*tf.reduce_mean(reg_loss)
    
    prediction_GT=tf.argmax(prob_GT, 1)
    prediction_A=tf.argmin(prob_A, 1)
    prediction_G=tf.argmin(prob_G, 1)
    prediction_T=tf.argmin(prob_T, 1)
    prediction_C=tf.argmin(prob_C, 1)
    '''
    
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    

    
    prediction_accuracy_A= tf.equal(prediction_A, tf.argmax(A_label, 1))
    prediction_accuracy_G= tf.equal(prediction_G, tf.argmax(G_label, 1))
    prediction_accuracy_T= tf.equal(prediction_T, tf.argmax(T_label, 1))
    prediction_accuracy_C= tf.equal(prediction_C, tf.argmax(C_label, 1))
    prediction_accuracy_GT=tf.equal(prediction_GT, tf.argmax(GT_label, 1))
    
    accuracy_GT = tf.reduce_sum(tf.cast(prediction_accuracy_GT, tf.float32))
    accuracy_A = tf.reduce_sum(tf.cast(prediction_accuracy_A, tf.float32))
    accuracy_G = tf.reduce_sum(tf.cast(prediction_accuracy_G, tf.float32))
    accuracy_T = tf.reduce_sum(tf.cast(prediction_accuracy_T, tf.float32))
    accuracy_C = tf.reduce_sum(tf.cast(prediction_accuracy_C, tf.float32))
    
    
    tmp1=tf.logical_and(prediction_accuracy_GT, prediction_accuracy_A)
    tmp2=tf.logical_and(prediction_accuracy_G, prediction_accuracy_T)
    tmp3=tf.logical_and(tmp1,tmp2)
    accuracy=tf.reduce_sum(tf.cast(tf.logical_and(prediction_accuracy_C,tmp3), tf.float32))
    
    tensors=(x,GT_label,A_label, G_label, T_label, C_label,GT_score, A_score, G_score, T_score, C_score, accuracy_GT, accuracy_A,  accuracy_G,  accuracy_T,  accuracy_C, prediction_accuracy_GT, prediction_accuracy_A,  prediction_accuracy_G,  prediction_accuracy_T,  prediction_accuracy_C, prediction_GT, prediction_A,  prediction_G,  prediction_T,  prediction_C, accuracy, cost, optimizer,  cost_GT, cost_A, cost_G, cost_T, cost_C,A_ref,G_ref,T_ref,C_ref,prob_GT,prob_A,prob_G,prob_T,prob_C,keep)

    
    
    return weights,biases,tensors

'''
def get_tensors_2(n_input,learning_rate=0):
    h_in,w_in,depth=n_input
    
    #h_in=int(np.ceil(float(h_in-2+1)/float(1)))
    #w_in=int(np.ceil(float(w_in-4+1)/float(1)))
    #h_in=28
    for i in range(2):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(2)))
       
    
    a,b,c=8,8,8
    d,e,f=32,48,32
    g,h=32,16
    i=0
    a,b,c=8,8,8
    d,e,f=32,48,32
    g,h=16,16
    i=0
    weights = {
        'wc1_1': tf.get_variable('W0_1', shape=(1,5,depth,a), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc1_2': tf.get_variable('W0_2', shape=(33,1,depth,b), initializer=tf.contrib.layers.xavier_initializer()),
        'wc1_3': tf.get_variable('W0_3', shape=(3,3,depth,c), initializer=tf.contrib.layers.xavier_initializer()),
        'wc2': tf.get_variable('W1', shape=(3,3,a+b+c,d), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,d,e), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*e,f), initializer=tf.contrib.layers.xavier_initializer()), 
        'w-tr': tf.get_variable('W4', shape=(i+f,g), initializer=tf.contrib.layers.xavier_initializer()),
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
        'bd-gtp': tf.get_variable('B6', shape=(h), initializer=tf.contrib.layers.xavier_initializer()),
        'out-gtp': tf.get_variable('B7', shape=(2), initializer=tf.contrib.layers.xavier_initializer()),
        
    }

    
    x,GT_label = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, 2])

    
    keep = tf.placeholder(tf.float32)
    
    GT_score= conv_net(x,weights, biases, keep)
    
    cost=tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=GT_score, labels=GT_label))
    
    prediction_GT=tf.argmax(prob_GT, 1)
    prediction_accuracy_GT=tf.equal(prediction_GT, tf.argmax(GT_label, 1))
    
    accuracy = tf.reduce_mean(tf.cast(prediction_accuracy_GT, tf.float32))
    
    
    tensors=(x,GT_label, accuracy, cost, optimizer, keep)

    
    
    return weights,biases,tensors'''