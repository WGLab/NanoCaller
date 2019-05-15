import time,os,copy,argparse,utils,gzip
import pandas as pd
import numpy as np
import multiprocessing as mp
import multi_gcp as mgcp
import tensorflow as tf
import generate_candidate_pileups as gcp


def conv2d(x, W, b, strides=1):
    x = tf.nn.conv2d(x, W, strides=[1, strides, strides, 1], padding='VALID')
    x = tf.nn.bias_add(x, b)
    return tf.nn.relu(x)

def conv_net(x,ref, weights, biases):  
    conv1 = conv2d(x, weights['wc1'], biases['bc1'],2)

    conv2 = conv2d(conv1, weights['wc2'], biases['bc2'],2)

    conv3 = conv2d(conv2, weights['wc3'], biases['bc3'],2)

    fc1 = tf.reshape(conv3, [-1, weights['wd1'].get_shape().as_list()[0]])
    fc1 = tf.add(tf.matmul(fc1, weights['wd1']), biases['bd1'])
    fc1 = tf.nn.relu(fc1)

    out = tf.add(tf.matmul(fc1, weights['out']), biases['out'])
    
    fc2=tf.concat([out,fc1,ref],1)
    fc2 = tf.add(tf.matmul(fc2, weights['wd2']), biases['bd2'])
    fc2 = tf.nn.relu(fc2)
    
    out_prob= tf.add(tf.matmul(fc2, weights['out_prob']), biases['out_prob'])
    
    return out,out_prob

def get_data(file_path,dims,window=None):
    pileups=utils.read_pileups_from_file(file_path,dims)
    #pileups=(int(pos),p_mat,int(gtype[0]),int(allele),int(ref),rnames)
    
    pos=np.array([x[0] for x in pileups.values()])
    
    mat=np.array([x[1] for x in pileups.values()])
    
    gt=np.array([x[2] for x in pileups.values()])
    gt=np.eye(3)[gt]
    
    allele=np.array([x[3] for x in pileups.values()])
    allele=np.eye(4)[allele]
    
    ref=np.array([x[4] for x in pileups.values()])
    ref=np.eye(4)[ref]
    

    
    
    if window:
        n=int((x.shape[2]-1)/2)
        x=x[:,:,n-window:n+1+window,:]
    return (pos.astype(int),mat.astype(np.int8),gt.astype(bool),allele.astype(bool),ref.astype(bool))
    
    
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
        'wd2': tf.get_variable('W7', shape=(135,32), initializer=tf.contrib.layers.xavier_initializer()), 
        'out_prob': tf.get_variable('W8', shape=(32,4), initializer=tf.contrib.layers.xavier_initializer()), 

    
    }
    biases = {
        'bc1': tf.get_variable('B0', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'bc2': tf.get_variable('B1', shape=(64), initializer=tf.contrib.layers.xavier_initializer()),
        'bc3': tf.get_variable('B2', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'bd1': tf.get_variable('B3', shape=(128), initializer=tf.contrib.layers.xavier_initializer()),
        'out': tf.get_variable('B4', shape=(n_classes), initializer=tf.contrib.layers.xavier_initializer()),
        'bd2': tf.get_variable('B5', shape=(32), initializer=tf.contrib.layers.xavier_initializer()),
        'out_prob': tf.get_variable('B6', shape=(4), initializer=tf.contrib.layers.xavier_initializer()),

    
    }

    x,y = tf.placeholder("float", [None]+n_input), tf.placeholder("float", [None, n_classes])
    allele,ref=tf.placeholder("float", [None, 4]),tf.placeholder("float", [None, 4])
    
    pred,fc_layer = conv_net(x,ref,weights, biases)
    
    gt_likelihood=tf.nn.softmax(logits=pred)
    allele_likelihood=tf.nn.softmax(logits=fc_layer)
    
    cost = tf.add(tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=pred,\
    labels=y)),tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=fc_layer, labels=allele)))
    
    optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)
    
    
    correct_prediction_gt = tf.equal(tf.argmax(pred, 1), tf.argmax(y, 1))
    accuracy_gt = tf.reduce_mean(tf.cast(correct_prediction_gt, tf.float32))
    
    
    correct_prediction_allele = tf.equal(tf.argmax(fc_layer, 1), tf.argmax(allele, 1))
    accuracy_allele = tf.reduce_mean(tf.cast(correct_prediction_allele, tf.float32))
    
    t1=(x,y,allele,ref,fc_layer,pred,cost,optimizer)
    t2=(correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)
    return weights,biases,t1,t2


'''def test_model(params):
    model_path,test_path,n_classes,window= params['model'], params['test_path'], params['classes'], params['window']
    
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
    print('Het: Precision=%.2f,  Recall=%.2f' %(prec2,rec2))'''

    
    
def genotype_caller_skinny(params,input_type='path',data=None):
    tf.reset_default_graph()
    n_input=params['dims']
    if input_type=='path':
        x_train,y_train,train_allele,train_ref= get_data(params['train_path'],n_input,params['window'])
        print('train data received')

    else:
        x_train,y_train,train_allele,train_ref=data['train_data']
        
    training_iters, learning_rate, batch_size, n_classes,plot= params['iters'],\
    params['rate'], params['size'], params['classes'],params['plot']

    n_input=[i for i in x_train.shape[1:]]

    weights,biases,t1,t2=get_tensors(n_input,n_classes,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2

    init = tf.global_variables_initializer()
    saver = tf.train.Saver(write_version=tf.train.SaverDef.V2)


    with tf.Session() as sess:
        sess.run(init)
        sess.run(tf.local_variables_initializer())
        stats=[]
        print('starting training')
        for i in range(training_iters):
            for batch in range(len(x_train)//batch_size):
                batch_x = x_train[batch*batch_size:min((batch+1)*batch_size,len(x_train))]
                batch_y = y_train[batch*batch_size:min((batch+1)*batch_size,len(y_train))]    
                batch_ref = train_ref[batch*batch_size:min((batch+1)*batch_size,len(train_ref))]
                batch_allele = train_allele[batch*batch_size:min((batch+1)*batch_size,len(train_allele))]
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele})
        saver.save(sess, save_path=params['model'])
    #(int(pos),p_mat,int(gtype[0]),int(allele),int(ref)
def genotype_caller(params,input_type='path',data=None):
    tf.reset_default_graph()
    n_input=params['dims']
    if input_type=='path':
        _,x_train,y_train,train_allele,train_ref= get_data(params['train_path'],n_input,params['window'])
        print('train data received')
        _,x_test,y_test,test_allele,test_ref= get_data(params['test_path'],n_input,params['window'])
        print('test data received')

    else:
        x_train,y_train,train_allele,train_ref=data['train_data']
        x_test,y_test,test_allele,test_ref=data['test_data']

    training_iters, learning_rate, batch_size, n_classes,plot= params['iters'],\
    params['rate'], params['size'], params['classes'],params['plot']

    n_input=[i for i in x_train.shape[1:]]

    weights,biases,t1,t2=get_tensors(n_input,n_classes,learning_rate)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2

    init = tf.global_variables_initializer()
    saver = tf.train.Saver()


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
                batch_ref = train_ref[batch*batch_size:min((batch+1)*batch_size,len(train_ref))]
                batch_allele = train_allele[batch*batch_size:min((batch+1)*batch_size,len(train_allele))]
                # Run optimization op (backprop).
                    # Calculate batch loss and accuracy
                opt = sess.run(optimizer, feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele})
                loss, acc_gt,acc_allele = sess.run([cost, accuracy_gt,accuracy_allele], feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele})
            print("Iter " + str(i) + ", Loss= " + \
                          "{:.6f}".format(loss) + ", Training GT Accuracy= " + \
                          "{:.5f}".format(acc_gt)+",\n Training Allele Accuracy= " + \
                          "{:.5f}".format(acc_allele))
            print("Optimization Finished!\n")

            # Calculate accuracy for all 10000 mnist test images
            test_acc_gt,test_acc_allele,valid_loss = sess.run([accuracy_gt,accuracy_allele,cost],\
                                           feed_dict={x: x_test,y : y_test,ref:test_ref,allele:test_allele})
                # repeat steps 4 for the histogram summary            
            stats.append([loss,valid_loss,acc_gt,test_acc_gt,acc_allele,test_acc_allele])
            print("Testing GT Accuracy= "+"{:.5f}".format(test_acc_gt)+",\n Testing Allele Accuracy= "+"{:.5f}".format(test_acc_allele))
            print('--------------------------------------------------\n')
        stats=np.array(stats)
        fc_layer_train,score_train,loss, acc = sess.run([fc_layer,pred,cost, accuracy_gt],\
                                                                      feed_dict={x: x_train,y: y_train,ref:train_ref,allele:train_allele})
        fc_layer_test,score_test,test_acc,valid_loss = sess.run([fc_layer,pred,accuracy_gt,cost],\
                                           feed_dict={x: x_test,y : y_test,ref:test_ref,allele:test_allele})


        res=np.array([np.argmax(score_test,axis=1),np.argmax(y_test,axis=1)]).transpose()

        rec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,1]==1)
        prec1=100*sum(np.all(res[:,:]==1,axis=1))/sum(res[:,0]==1)

        rec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,1]==2)
        prec2=100*sum(np.all(res[:,:]==2,axis=1))/sum(res[:,0]==2)

        print('Hom-Alt: Precision=%.2f,  Recall=%.2f' %(prec1,rec1))
        print('Het: Precision=%.2f,  Recall=%.2f' %(prec2,rec2))

        saver.save(sess, save_path=params['model'])
        writer.close()
    totes=np.hstack([tmp,tmp2])
    m=np.sum((totes[:,0]==totes[:,1])&(totes[:,2]==totes[:,3])&(totes[:,3]!=0))/np.sum(totes[:,3]!=0)
    print('Recall=%.4f' %m)

    n=np.sum((totes[:,0]==totes[:,1])&(totes[:,2]==totes[:,3])&(totes[:,2]!=0))/np.sum(totes[:,2]!=0)
    print('Precision=%.4f' %n)

    f1=2*n*m/(n+m)
    print('F1=%.4f'%f1)
    with open('cnn_result','w') as print_file:
        for item in stats:
            print_file.write("%s\n" % item)
        print_file.write('Recall= %.4f , Precision= %.4f , F1= %.4f' %(m,n,f1))
    return (fc_layer_test,score_test)

def test_model(params):
    model_path,test_path,n_input,n_classes,chrom,vcf_path= params['model'], params['test_path'],params['dims'], 3,params['chrom'],params['vcf_path']
    tf.reset_default_graph()
    pos,x_test,y_test,test_allele,test_ref= get_data(params['test_path'],n_input)
    print('test data received')

    n_input=[i for i in x_test.shape[1:]]
    
    weights,biases,t1,t2=get_tensors(n_input,n_classes,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    
    rev_mapping={0:'A',1:'G',2:'T',3:'C'}
    gt_map={1:'1|1',2:'0|1'}
    
    batch_size=1000

    with open(vcf_path,'w') as f:

        f.write('##fileformat=VCFv4.2\n')
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        c='##contig=<ID=%s>\n' %chrom
        f.write('##contig=<ID=%s>\n' %chrom)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Consensus Genotype across all datasets with called genotype">\n')
        f.write('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE\n')
        ttt=[]
        for batch in range(len(x_test)//batch_size+1):
                    batch_pos = pos[batch*batch_size:min((batch+1)*batch_size,len(pos))]
                    batch_x = x_test[batch*batch_size:min((batch+1)*batch_size,len(x_test))]
                    batch_y = y_test[batch*batch_size:min((batch+1)*batch_size,len(y_test))]    
                    batch_ref = test_ref[batch*batch_size:min((batch+1)*batch_size, len(test_ref))]
                    batch_allele = test_allele[batch*batch_size:min((batch+1)*batch_size, len(test_allele))]
                    # Run optimization op (backprop).
                        # Calculate batch loss and accuracy


                    fc_layer_batch,score_batch,gl,qual = sess.run([fc_layer,pred,gt_likelihood,allele_likelihood],\
                                               feed_dict={x: batch_x,y: batch_y,ref:batch_ref,allele:batch_allele})
                    
                    gl=10*np.log(np.min((1-gl),axis=1))
                    qual=10*np.log(np.min((1-qual),axis=1))
                    all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)

                    mat=np.hstack([batch_pos[:,np.newaxis],np.argmax(batch_ref,axis=1)\
                                             [:,np.newaxis],all_pred[:,np.newaxis],gt_pred[:,np.newaxis],qual[:,np.newaxis],\
                                   gl[:,np.newaxis]])
                    mat=mat[(mat[:,3]!=0) & (mat[:,1]!=mat[:,2])]
                    for j in range(len(mat)):
                        v=mat[j]
                        s='%s\t%d\t.\t%s\t%s\t%.3f\tPASS\t.\tGT:GL\t%s:%.3f\n' %\
                        (chrom,v[0],rev_mapping[v[1]],rev_mapping[v[2]],v[4],gt_map[v[3]],v[5])
                        f.write(s)
                    ttt.append(qual2)
    return np.vstack(ttt)  

        
def test_model_bam(params):
    model_path,n_classes,n_input= params['model'], params['classes'], params['dims']
    start,end=params['start'],params['end']
    chrom=params['chrom']
    start=params['start']
    end=params['end']
    sam_path=params['sam_path']
    fasta_path=params['fasta_path']
    
    tf.reset_default_graph()
    
    
    weights,biases,t1,t2=get_tensors(n_input,n_classes,1)
    (x,y,allele,ref,fc_layer,pred,cost,optimizer)=t1
    (correct_prediction_gt,correct_prediction_allele,accuracy_gt,accuracy_allele,gt_likelihood,allele_likelihood)=t2
    
    init = tf.global_variables_initializer()
    sess = tf.Session()
    sess.run(init)
    sess.run(tf.local_variables_initializer())
    saver = tf.train.Saver()
    saver.restore(sess, model_path)
    total=[]
    batch_size=50
    rec_,prec_=[0,0],[0,0]
    #file=open('bam_vcf' , "w")
    
    for batch in range(start,end,1000):
        batch_pos,batch_x,batch_ref =mgcp.create_testing_pileup({'chrom':chrom,'start':batch,'end':batch+1000,\
                                                                 'depth':params['depth'],'window':params['window'],'sam_path':sam_path,'fasta_path':fasta_path})
               # Run optimization op (backprop).
            # Calculate batch loss and accuracy


        fc_layer_batch,score_batch = sess.run([fc_layer,pred],\
                                   feed_dict={x: batch_x,ref:batch_ref})

        all_pred,gt_pred=np.argmax(fc_layer_batch,axis=1),np.argmax(score_batch,axis=1)
        total.append([np.hstack([batch_pos[:,np.newaxis],np.argmax(batch_ref,axis=1)[:,np.newaxis],all_pred[:,np.newaxis],gt_pred[:,np.newaxis]])])
        
    return total    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-r", "--rate", help="Learning rate",type=float)
    parser.add_argument("-i", "--iterations", help="Training iterations",type=int)
    parser.add_argument("-s", "--size", help="Batch size",type=int)
    parser.add_argument("-n", "--classes", help="Number of classes",type=int)
    parser.add_argument("-w", "--window", help="Window size",type=int)
    parser.add_argument("-train", "--train", help="Train path")
    parser.add_argument("-test", "--test", help="Test path")
    parser.add_argument("-p", "--plot", help="Show plots",type=int)
    parser.add_argument("-model", "--model", help="Model output path")
    parser.add_argument("-m", "--mode", help="Mode")
    parser.add_argument("-dim", "--dimensions", help="Input dimensions")
    parser.add_argument("-vcf", "--vcf", help="VCF output path")
    parser.add_argument("-chrom", "--chrom", help="Chromosome")
    
    args = parser.parse_args()
    input_dims=[int(x) for x in args.dimensions.split(':')]
    t=time.time()
    
    
    if args.mode=='train':
        in_dict={'rate':args.rate,'iters':args.iterations,'size':args.size,\
                 'classes':args.classes,'window':args.window,'dims':input_dims,\
                 'plot':args.plot,'train_path':args.train,'test_path':args.test,'model':args.model}
        genotype_caller_skinny(in_dict)
    
    else:
        in_dict={'dims':input_dims,'test_path':args.test,'model':args.model,'chrom':args.chrom,'vcf_path':args.vcf}
        test_model(in_dict)
        
    elapsed=time.time()-t
    print ('Total Time Elapsed: %.2f seconds' %elapsed)
