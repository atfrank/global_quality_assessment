import pandas as pd
import numpy as np
import sys
import timeit
from itertools import compress
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,recall_score
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras import metrics
from keras.optimizers import RMSprop,Adam,SGD
from keras.constraints import maxnorm
from keras import backend as K
import tensorflow as tf
from tensorflow.python.client import device_lib
from keras.models import load_model
from keras.wrappers.scikit_learn import KerasClassifier
from keras import regularizers



def create_model(neurons, activation, init, l2, lr, loss, opt):
  model = Sequential()
  model.add(Dense(neurons, input_dim = 58, activation = activation, 
  kernel_initializer = init, kernel_regularizer = regularizers.l2(l2)))
  model.add(Dense(1, activation = 'sigmoid'))
  if opt == 'Adam':
    optimizer = Adam(lr = lr)
  if opt == 'RMSprop':
    optimizer = RMSprop(lr = lr)
  model.compile(loss = loss, optimizer = optimizer, metrics = ['accuracy'])
  return model

def main():
  K.tensorflow_backend._get_available_gpus()
  DIR_PATH = '/home/kexin/projects/global_quality_assessment/neural_net_classifier/'  
  # LOAD NORMALIZED ERROR MATRIX
  errors = pd.read_table(DIR_PATH+'errors_normalized.txt',sep=' ')
  
  # LEAVE-ONE-OUT ANLAYSIS
  rnas = ["1XHP", "1YSV", "1Z2J", "1ZC5", "28SR", "2FDT", "2JWV", "2K66", "2KOC", "2L1V", 
          "2L3E", "2LBJ", "2LBL", "2LDL", "2LDT", "2LHP", "2LI4", "2LK3", "2LP9", "2LPA", 
          "2LPS", "2LQZ", "2LU0", "2LUB", "2LUN", "2LV0", "2M12", "2M21", "2M22", "2M24", 
          "2M4W", "2M5U", "2M8K", "2MEQ", "2MHI", "2MNC", "2MXL", "2N2O", "2N2P", "2N4L", 
          "2N6S", "2N6T", "2N6W", "2N6X", "2NCI", "2QH2", "2QH4", "2RVO", "2Y95", "4A4S", 
          "4A4T", "4A4U", "5A17", "5A18", "5IEM", "5KH8", "5KMZ", "5KQE", "5LSN", "5LWJ", 
          "5N5C", "5UF3", "5UZT", "5V16", "5V17", "5WQ1", "6EZ0"]
  
  # remove similar sequences
  sim = pd.read_table(DIR_PATH+'similarity.txt',sep=' ')
  threshold = 0.80
  rna = sys.argv[1]
  print('target RNA: ', rna)
  fil = sim[rna]>threshold
  remove_rnas = list(compress(rnas, fil.values))
  train_rnas = [e for e in rnas if e not in remove_rnas]
  print('sequences removed: ', remove_rnas)
  
  train = errors.loc[errors['id'].isin(train_rnas)]
  test = errors[errors['id']==rna]
  
  # scramble data just to be sure
  train = train.sample(len(train),random_state=1)
  test = test.sample(len(test),random_state=1)
  
  # make copy of train and test actual flag
  train_result = train[['id','model','flag']]
  test_result = test[['id','model','flag']]
  X_train = train.drop(columns=['id','model','flag','rmsd'])
  X_test = test.drop(columns=['id','model','flag','rmsd'])
  y_train = train['flag']
  y_test = test['flag']
  
  del train, test
  
  # BUILD CLASSIFIER
  # grid search parameters
  neurons = [24, 36, 48]
  activation = ['relu']
  init = ['normal']
  l2 = [0.1, 0.01]
  lr = [0.01, 0.001]
  loss = ['binary_crossentropy']
  opt = ['Adam', 'RMSprop']
  epochs = [25, 50]
  batch_size = [128, 256]
  # grid search cv 
  param_grid = dict(opt = opt, epochs = epochs, init = init, loss = loss, neurons = neurons, batch_size = batch_size, lr = lr, l2 = l2, activation = activation)
  with tf.device('/cpu:0'):
    model = KerasClassifier(build_fn=create_model,batch_size=batch_size,epochs=epochs,verbose=1)
  clf = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, n_jobs=1, 
                     return_train_score=True, scoring='accuracy', verbose = 2)
  grid_result = clf.fit(X_train,y_train)
  
  # write log file for neural network information
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  logfile = DIR_PATH+'log/nn_log_'+rna+'.txt'
  with open(logfile, 'a') as f:
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  
  # predict and format results
  y_pred_train = clf.predict(X_train).flatten().tolist()
  train_result.loc[:,'pred'] = y_pred_train
  y_pred_test = clf.predict(X_test).flatten().tolist()
  test_result.loc[:,'pred'] = y_pred_test
  # print confusion matrix
  with open(logfile, 'a') as f:
    print(classification_report(y_test, y_pred_test), file=f)
    print(confusion_matrix(y_test, y_pred_test), file=f)
  train_result.to_csv(DIR_PATH + 'pred/train_pred_'+rna+'.txt', sep = ' ', index = False)
  test_result.to_csv(DIR_PATH + 'pred/test_pred_'+rna+'.txt', sep = ' ', index = False)
  # save model
  best_model = clf.best_estimator_
  modelname = DIR_PATH + 'models/nn_classifier_'+rna+'.h5'
  best_model.model.save(modelname)
  
if __name__ == "__main__":
  main()
                 
  
  
  
  
  
  
  
