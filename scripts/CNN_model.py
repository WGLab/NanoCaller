import time,os,copy,argparse,utils,gzip
import pandas as pd
import numpy as np
import multiprocessing as mp

import tensorflow as tf
import generate_candidate_pileups as gcp


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

def get_data(file_path,window=None):
    pileups=utils.read_pileups_from_file(file_path)
    #pileups=[[pos,gtype,rnames,p_mat]]

    x=np.array([x[-1] for x in pileups.values()])
    y=np.array([x[1] for x in pileups.values()])
    y=np.eye(3)[y]
    if window:
        n=int((x.shape[2]-1)/2)
        x=x[:,:,n-window:n+1+window,:]
    return x,y
    
    
def get_tensors(n_input,n_classes,learning_rate=0):
    h_in,w_in,depth=n_input
    for i in range(3):
        h_in=int(np.ceil(float(h_in-3+1)/float(2)))
        w_in=int(np.ceil(float(w_in-3+1)/float(2)))
        
    weights = {
        'wc1': tf.get_variable('W0', shape=(3,3,depth,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc2': tf.get_variable('W1', shape=(3,3,32,64), initializer=tf.contrib.layers.xavier_initializer()), 
        'wc3': tf.get_variable('W2', shape=(3,3,64,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'wd1': tf.get_variable('W3', shape=(h_in*w_in*128,128), initializer=tf.contrib.layers.xavier_initializer()), 
        'out': tf.get_variable('W6', shape=(128,n_classes), initializer=tf.contrib.layers.xavier_initializer()), 
    }
    biases = {
        'bc1': tf.get_variable('B0', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'out': tf.get_variable('B4', shape=(n_classes), initializer=tf.contrib.layers.xavier_initializer()),
    }

    x,y = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, n_classes])

    pred = conv_net(x, weights, biases)
    cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred, labels=y))
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    correct_prediction = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
    
    t1=(x,y,pred,cost,optimizer,correct_prediction,accuracy)
    return weights,biases,t1
    
def genotype_caller(params,input_type='path',training_data=None,testing_data=None):
    tf.reset_default_graph()
    if input_type=='path':
        x_train,y_train=get_data(params['train_path'],params['window'])
        print('train data received')
        x_test,y_test=get_data(params['test_path'],params['window'])
        print('test data received')
    
    else:
        x_train,y_train=training_data
        x_test,y_test=testing_data
        
    training_iters, learning_rate, batch_size, n_classes,plot= params['iters'],\
    params['rate'], params['size'], params['classes'],params['plot']
    
    n_input=[i for i in x_train.shape[1:]]

    weights,biases,t1=get_tensors(n_input,n_classes,learning_rate)
    x,y,pred,cost,optimizer,correct_prediction,accuracy=t1
    
    
    testing_summary=[tf.summary.scalar('Testing_Loss', cost),tf.summary.scalar('Testing_Accuracy', accuracy)]
    
    
    summary_list=[tf.summary.scalar('Training_Loss', cost),\
                  tf.summary.scalar('Training_Accuracy', accuracy),\
                  tf.summary.histogram('Cnv1_weight',weights['wc1']),\
                  tf.summary.histogram('Cnv2_weight',weights['wc2']),\
                  tf.summary.histogram('Cnv3_weight',weights['wc3']),\
                  tf.summary.histogram('Full_weight',weights['wd1']),\
                  tf.summary.histogram('Out_weight',weights['out'])]

    merged = tf.summary.merge(summary_list)
    testing_merged=tf.summary.merge(testing_summary)
    
    saver = tf.train.Saver()
    init = tf.global_variables_initializer()
    count=0
    with tf.Session() as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        writer = tf.summary.FileWriter('./Output', sess.graph)
        print('starting training')
        for i in range(training_iters):
            for batch in range(len(x_train)//batch_size):
                batch_x = x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))]
                batch_y = y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))]    
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y})
                loss, acc,merge_summary = sess.run([cost, accuracy,merged], feed_dict={x: batch_x,y: batch_y})
                writer.add_summary(merge_summary,count)
                count+=1
            print("Iter " + str(i) + ", Loss= " + \
                          "{:.6f}".format(loss) + ", Training Accuracy= " + \
                          "{:.5f}".format(acc))
            print("Optimization Finished!")

            # Calculate accuracy for all 10000 mnist test images
            test_acc,valid_loss,test_summary = sess.run([accuracy,cost,testing_merged],\
                                           feed_dict={x: x_test,y : y_test})
            writer.add_summary(test_summary,i)
                # repeat steps 4 for the histogram summary            
            stats.append([loss,valid_loss,acc,test_acc])
            print("Testing Accuracy:","{:.5f}".format(test_acc))
        
        stats=np.array(stats)
        
        score,test_acc,valid_loss = sess.run([pred,accuracy,cost],\
                                           feed_dict={x: x_test,y : y_test})
        
        res=np.array([np.argmax(score,axis=1),np.argmax(y_test,axis=1)]).transpose()

        rec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,1]==1)
        prec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,0]==1)

        rec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,1]==2)
        prec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,0]==2)

        print('Hom-Alt: Precision=%.2f,  Recall=%.2f' %(prec1,rec1))
        print('Het: Precision=%.2f,  Recall=%.2f' %(prec2,rec2))
        
        saver.save(sess, save_path=model_path)
        writer.close()
    if plot:
        utils.plot_training_stats(stats)
    return (score,stats,[prec1,rec1,prec2,rec2])

def test_model(model_path,test_path,n_classes,window=None):
    
    parser.add_argument("-n", "--classes", help="Number of classes",type=int)
    parser.add_argument("-w", "--window", help="Window size",type=int)
    parser.add_argument("-model", "--model", help="Model path")
    parser.add_argument("-test", "--test", help="Test path")

    args = parser.parse_args()
    
    model_path,test_path,n_classes,window=args.model,args.test,args.classes,args.window
    
    tf.reset_default_graph()
    x_test,y_test=get_data(test_path,window)
    n_input=[i for i in x_test.shape[1:]]

    weights,biases,t1=get_tensors(n_input,n_classes)
    x,y,pred,cost,optimizer,correct_prediction,accuracy=t1
    
    sess = tf.Session()
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    score,test_acc,valid_loss = sess.run([pred,accuracy,cost],\
                                           feed_dict={x: x_test,y : y_test})
        
    res=np.array([np.argmax(score,axis=1),np.argmax(y_test,axis=1)]).transpose()

    rec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,1]==1)
    prec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,0]==1)

    rec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,1]==2)
    prec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,0]==2)

    print('Hom-Alt: Precision=%.2f,  Recall=%.2f' %(prec1,rec1))
    print('Het: Precision=%.2f,  Recall=%.2f' %(prec2,rec2))
    
def train_model(parser):
    parser.add_argument("-r", "--rate", help="Learning rate",type=float)
    parser.add_argument("-i", "--iterations", help="Training iterations",type=int)
    parser.add_argument("-s", "--size", help="Batch size",type=int)
    parser.add_argument("-n", "--classes", help="Number of classes",type=int)
    parser.add_argument("-w", "--window", help="Window size",type=int)
    parser.add_argument("-train", "--train", help="Train path")
    parser.add_argument("-test", "--test", help="Test path")
    parser.add_argument("-p", "--plot", help="Show plots",type=int)
    
    args = parser.parse_args()
        
    in_dict={'rate':args.rate,'iters':args.iterations,'size':args.size,'classes':args.classes,'window':args.window,'plot':args.plot,'train_path':args.train,'test_path':args.test,'output':args.output}
    
    genotype_caller(in_dict)

    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mode", help="Mode")
    t=time.time()
    
    if parser.parse_args().mode=='train':
        train_model(parser)
    else:
        test_model(parser)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)