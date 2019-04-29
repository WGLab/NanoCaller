import sys,pysam, time,os,re,copy,argparse
from collections import Counter
import pandas as pd
import numpy as np
import multiprocessing as mp
from pysam import VariantFile
import sklearn
from sklearn.model_selection import train_test_split
import tensorflow as tf
import generate_candidate_pileups as gcp
from matplotlib import pyplot as plt


def conv2d(x, W, b, strides=1):
    x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding='VALID')
    x = tf.nn.bias_add(x, b)
    return tf.nn.relu(x)

def conv_net(x, weights, biases):  
    conv1 = conv2d(x, weights['wc1'], biases['bc1'],2)

    conv2 = conv2d(conv1, weights['wc2'], biases['bc2'],2)

    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],2)

    fc1 = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    fc1 = tf.add(tf.matmul(fc1, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.relu(fc1)

    out = tf.add(tf.matmul(fc1, weights['out']), biases['out'])
    return out
    
def genotype_caller(train_path,test_path):

    labelmap={'hom-ref':0,'hom-alt':1,'het':2}

    train_pileups=gcp.read_pileups_from_file(train_path)
    test_pileups=gcp.read_pileups_from_file(test_path)
    #pileups=[[pos,gtype,rnames,p_mat]]

    x_train=np.array([x[-1] for x in train_pileups.values()])
    y_train=np.array([labelmap[x[1]] for x in train_pileups.values()])
    y_train=np.eye(3)[y_train]

    x_test=np.array([x[-1] for x in test_pileups.values()])
    y_test=np.array([labelmap[x[1]] for x in test_pileups.values()])
    y_test=np.eye(3)[y_test]

    training_iters = 30
    learning_rate = 0.001 
    batch_size = 20
    n_input = [32, 101, 8]

    n_classes = 3

    tf.reset_default_graph()

    weights = {
        'wc1': tf.get_variable('W0', shape=(3,3,8,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc2': tf.get_variable('W1', shape=(3,3,32,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,64,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(3*11*128,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'out': tf.get_variable('W6', shape=(128,n_classes), initializer=tf.contrib.layers.xavier_initializer()), 
    }
    biases = {
        'bc1': tf.get_variable('B0', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'out': tf.get_variable('B4', shape=(n_classes), initializer=tf.contrib.layers.xavier_initializer()),
    }

    x = tf.placeholder("float", [None]+n_input)
    y = tf.placeholder("float", [None, n_classes])

    pred = conv_net(x, weights, biases)

    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))

    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    init = tf.global_variables_initializer()



    with tf.Session() as sess:
        sess.run(init) 
        train_loss = []
        test_loss = []
        train_accuracy = []
        test_accuracy = []

        summary_writer = tf.summary.FileWriter('./Output', sess.graph)
        for i in range(training_iters):
            for batch in range(len(x_train)//batch_size):
                batch_x = x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))]
                batch_y = y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))]    
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y})
                loss, acc = sess.run([cost, accuracy], feed_dict={x: batch_x,y: batch_y})
            print("Iter " + str(i) + ", Loss= " + \
                          "{:.6f}".format(loss) + ", Training Accuracy= " + \
                          "{:.5f}".format(acc))
            print("Optimization Finished!")

            # Calculate accuracy for all 10000 mnist test images
            test_acc,valid_loss = sess.run([accuracy,cost], feed_dict={x: x_test,y : y_test})
            train_loss.append(loss)
            test_loss.append(valid_loss)
            train_accuracy.append(acc)
            test_accuracy.append(test_acc)
            print("Testing Accuracy:","{:.5f}".format(test_acc))
        #saver.save(sess, save_path='./test-case.ckpt')
        summary_writer.close()

        plt.plot(range(len(train_loss)), train_loss, 'b', label='Training loss')
        plt.plot(range(len(train_loss)), test_loss, 'r', label='Test loss')
        plt.title('Training and Test loss')
        plt.xlabel('Epochs ',fontsize=16)
        plt.ylabel('Loss',fontsize=16)
        plt.legend()
        plt.figure()
        plt.show()
        fig1 = plt.gcf()
        fig1.savefig('loss.png',dpi=1000)

        plt.plot(range(len(train_loss)), train_accuracy, 'b', label='Training Accuracy')
        plt.plot(range(len(train_loss)), test_accuracy, 'r', label='Test Accuracy')
        plt.title('Training and Test Accuracy')
        plt.xlabel('Epochs ',fontsize=16)
        plt.ylabel('Loss',fontsize=16)
        plt.legend()
        plt.figure()
        plt.show()
        fig2 = plt.gcf()
        fig2.savefig('accuracy.png',dpi=1000)