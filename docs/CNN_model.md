# CNN_model.py

## genotype_caller(params,input_type='path',training_data=None,testing_data=None)
Inputs:

params - parameter dictionary containing following keys:
1. 'rate' : learning rate for neural net
2. 'iters' : number of epochs
3. 'size' : batch size
4. 'classes' : number of label classes
5. 'window' : how many base pairs on either side of candidate sites to include in pileup image
6. 'plot' : plot loss and accuracy graphs
7. 'train_path' : path to training pileups from generate_candidate_pileups
8. 'test_path'  : path to testing pileups from generate_candidate_pileups
9. 'model' : path to save tensorflow model

Outputs:

tuple - (score,stats,[prec1,rec1,prec2,rec2])

1. score : output from last fully connected layer before softmax function evaluated on test set
2. stats : list of [training accuracy,testing accuracy, training loss, testing loss]
3. prec1,rec1 : precision and recall for class 1 (hom-alt)
4. prec2,rec2 : precision and recall for class 2 (het)

This function accepts pileup files for training and testing, as well as hyperparameters for our CNN model. This neural network consists of
three convolutional layers and two fully connected layers, feeding into a softmax layer for predicting genotype. It also calculates 
precision and recall statistics for hom-alt and het genotypes, in addition to loss and accuiracy statistics. The trained model is saved is storage, as well as
tensorboard summaries in the subfolder 'Output' of the working directory.


## test_model(params)
Inputs:

params - parameter dictionary containing following keys:
1. 'classes' : number of label classes
2. 'window' : how many base pairs on either side of candidate sites to include in pileup image. Should be same as trhe one used to generate the saved model
3. 'test_path'  : path to testing pileups from generate_candidate_pileups
4. 'model' : path to restore tensorflow model

This function restores a a saved model and runs it on a given testing file to display precision and recall statistics.

### Usage
This module can be used in training and testing mode, which can be chosen by setting -m (--mode) flag to 'train' or 'test'.

From command line run the following for training mode:
`python CNN_model.py -r learning_rate -i epochs -s batch_size -n number_of_classes -w window_size -train path_to_training_file
-test path_to_testing_file -p binary_flag -model path_to_save_model -m train`

And the following for testing mode:
`python CNN_model.py -n number_of_classes -w window_size -test path_to_testing_file -model path_to_restore_model -m test`
