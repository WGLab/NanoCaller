import tensorflow as tf

from tensorflow.keras.layers import Dense, Flatten, Conv2D, Dropout, Softmax
from tensorflow.keras import Model

class SNP_model(Model):
    def __init__(self):
        super(SNP_model, self).__init__()
        self.conv1_1 = Conv2D(16, kernel_size=[1,5], strides=[1,1], activation='selu', name='C1_1', use_bias=True, padding='same')
        self.conv1_2 = Conv2D(16, kernel_size=[5,1], strides=[1,1], activation='selu', name='C1_2', use_bias=True, padding='same')
        self.conv1_3 = Conv2D(16, kernel_size=[5,5], strides=[1,1], activation='selu', name='C1_3', use_bias=True, padding='same')
        
        self.conv2 = Conv2D(32, kernel_size=[2,3], strides=[1,2], activation='selu', name='C2', use_bias=True, padding='valid')
        self.conv3 = Conv2D(64, kernel_size=[2,3], strides=[1,2], activation='selu', name='C3', use_bias=True, padding='valid')
        
        
        self.flatten = Flatten()
        self.dropout = Dropout(0.5)
        
        
        self.fc1 = Dense(48, activation='selu',name='C4',use_bias=True)
        
        self.fa  = Dense(16, activation='selu',name='C5',use_bias=True)
            
        self.A = Dense(2, activation=None ,name='C9',use_bias=True)
        self.G = Dense(2, activation=None ,name='C10',use_bias=True)
        self.T = Dense(2, activation=None ,name='C11',use_bias=True)
        self.C = Dense(2, activation=None ,name='C12',use_bias=True)
        
        self.fc2 = Dense(16, activation='selu',name='C6',use_bias=True)
        self.fc3 = Dense(8, activation='selu',name='C7',use_bias=True)
        self.GT = Dense(2, activation=None,name='C8',use_bias=True)
        
        self.softmax=Softmax()
        
    def call(self, inputs):
        x, A_ref, G_ref, T_ref, C_ref = inputs
        C1_1 = self.conv1_1(x)
        C1_2 = self.conv1_2(x)
        C1_3 = self.conv1_3(x)
        merge_conv_0=tf.concat([C1_1,C1_2,C1_3],3)
        
        
        C2   = self.conv2(merge_conv_0)
        C3   = self.conv3(C2)

        
        flat_nn = self.flatten(C3)
        
        fc1 =  self.fc1(flat_nn)
        drop_out_fc1 = self.dropout(fc1)
        
        fa  =  self.fa(drop_out_fc1)
        
        out_A = self.softmax(self.A(tf.concat([fa,A_ref],1)))
        out_G = self.softmax(self.G(tf.concat([fa,G_ref],1)))
        out_T = self.softmax(self.T(tf.concat([fa,T_ref],1)))
        out_C = self.softmax(self.C(tf.concat([fa,C_ref],1)))
        
        fc2   = self.fc2(drop_out_fc1)
        fc3   = self.fc3(tf.concat([fc2,out_A,out_G,out_T,out_C],1))
        out_GT= self.softmax(self.GT(fc3))
        
        return out_A,out_G,out_T,out_C,out_GT

