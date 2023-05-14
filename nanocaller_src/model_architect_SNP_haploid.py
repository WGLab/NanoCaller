import tensorflow as tf

from tensorflow.keras.layers import Dense, Flatten, Conv2D, Dropout, Softmax
from tensorflow.keras import Model
from tensorflow.keras import regularizers

class haploid_SNP_model(Model):
    def __init__(self):
        super(haploid_SNP_model, self).__init__()
        kernel_regularizer=regularizers.l2(l2=1e-3)
        bias_regularizer=regularizers.l2(1e-3)
        activity_regularizer=regularizers.l2(1e-3)
        
        self.conv1_1 = Conv2D(16, kernel_size=[1,5], strides=[1,1], activation='selu', name='C1_1', use_bias=True, padding='same', kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        self.conv1_2 = Conv2D(16, kernel_size=[5,1], strides=[1,1], activation='selu', name='C1_2', use_bias=True, padding='same', kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        self.conv1_3 = Conv2D(16, kernel_size=[5,5], strides=[1,1], activation='selu', name='C1_3', use_bias=True, padding='same', kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        
        self.conv2 = Conv2D(32, kernel_size=[2,3], strides=[1,2], activation='selu', name='C2', use_bias=True, padding='valid', kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        self.conv3 = Conv2D(64, kernel_size=[2,3], strides=[1,2], activation='selu', name='C3', use_bias=True, padding='valid', kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        
        
        self.flatten = Flatten()
        self.dropout = Dropout(0.5)
        
        
        self.fc1 = Dense(48, activation='selu',name='C4',use_bias=True, kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
                
        self.fc2 = Dense(16, activation='selu',name='C6',use_bias=True, kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        self.fc3 = Dense(4, activation='selu',name='C7',use_bias=True, kernel_regularizer=kernel_regularizer,bias_regularizer=bias_regularizer,activity_regularizer=activity_regularizer)
        
        self.softmax=Softmax()
        
    def call(self, inputs):
        x,  ref = inputs
        C1_1 = self.conv1_1(x)
        C1_2 = self.conv1_2(x)
        C1_3 = self.conv1_3(x)
        merge_conv_0=tf.concat([C1_1,C1_2,C1_3],3)
        
        
        C2   = self.conv2(merge_conv_0)
        C3   = self.conv3(C2)

        
        flat_nn = self.flatten(C3)
        
        fc1 =  self.fc1(flat_nn)
        drop_out_fc1 = self.dropout(fc1)
        fc2   = self.fc2(drop_out_fc1)
        fc3   = self.fc3(tf.concat([fc2, ref],1))
        out= self.softmax(fc3)
        
        return out

