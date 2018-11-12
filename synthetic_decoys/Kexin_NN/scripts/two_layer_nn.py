import pandas as pd
import numpy as np
import sys
import timeit

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

def norm_errors(errors):
  cols_1 = ['id', 'model', 'flag', 'rmsd']
  cols_2 = ['C1..GUA', 'C2..GUA', 'C3..GUA', 'C4..GUA', 'C5..GUA', 'C8.GUA', 'H1..GUA', 
            'H2..GUA', 'H3..GUA', 'H4..GUA', 'H5..GUA', 'H5...GUA', 'H2.GUA', 'H5.GUA', 
            'H6.GUA', 'H8.GUA', 'C1..ADE', 'C3..ADE', 'C5..ADE', 'C2.ADE', 'C5.ADE', 
            'C6.ADE', 'C8.ADE', 'H2..ADE', 'H4..ADE', 'H5..ADE', 'H2.ADE', 'H6.ADE', 
            'H8.ADE', 'C1..CYT', 'C2..CYT', 'C3..CYT', 'C5.CYT', 'C6.CYT', 'C8.CYT', 
            'H1..CYT', 'H2..CYT', 'H3..CYT', 'H4..CYT', 'H5..CYT', 'H5...CYT', 'H2.CYT', 
            'H5.CYT', 'H6.CYT', 'C1..URA', 'C3..URA', 'C4..URA', 'C5..URA', 'C2.URA', 
            'C5.URA', 'C6.URA', 'C8.URA', 'H1..URA', 'H2..URA', 'H4..URA', 'H5...URA', 
            'H2.URA', 'H6.URA']
  # standard scaling
  errors_1 = errors[cols_1]
  errors_2 = errors[cols_2]
  scaler = StandardScaler()
  errors_2 = scaler.fit_transform(errors_2)
  return pd.DataFrame(np.hstack([errors_1, errors_2]))
  

def create_model(neurons, neurons_2, activation, init, l2, lr, loss, opt):
  model = Sequential()
  model.add(Dense(neurons, input_dim = 58, activation = activation, 
  kernel_initializer = init, kernel_regularizer = regularizers.l2(l2)))
  model.add(Dense(neurons_2, activation = activation,
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
  # load errors matrices
  errors = pd.read_table(DIR_PATH+'errors.txt',sep=' ')
  
  # fill missing values with mean column values
  # errors.fillna(errors.mean(), inplace=True)
  # fill missing values with 0
  errors.fillna(0, inplace=True)
  
  # normalize the error for each RNA (specificed by id) separately
  cols = errors.columns
  errors = errors.groupby('id').apply(norm_errors)
  errors.columns = cols
  
  # scramble data just to be sure
  errors = errors.sample(len(errors),random_state = 1)
  
  # SPLIT INTO TRAINING AND TESTING
  #X = errors.iloc[:,4:61].copy()
  #y = errors.iloc[:,2].copy()
  X_train_tmp, X_test_tmp, y_train, y_test = train_test_split(errors, errors.iloc[:,2], test_size=0.25, random_state=42)
  train_result = X_train_tmp[['id','model','flag']]
  test_result = X_test_tmp[['id','model','flag']]
  X_train = X_train_tmp.drop(columns=['id','model','flag','rmsd'])
  X_train = X_train.apply(pd.to_numeric)
  X_test = X_test_tmp.drop(columns=['id','model','flag','rmsd'])
  # convert data types
  X_train = X_train.apply(pd.to_numeric)
  y_train = y_train.astype(int)
  y_test = y_test.astype(int)
 
  
  # BUILD CLASSIFIER
  # grid search parameters
  neurons = [48]
  neurons_2 = [12,24,36]
  activation = ['relu']
  init = ['normal']
  l2 = [0.1, 0.01]
  lr = [0.01, 0.001]
  #loss = ['binary_crossentropy','logcosh']
  loss = ['binary_crossentropy']
  opt = ['Adam', 'RMSprop']
  epochs = [25, 50]
  batch_size = [128, 256]
  
  # grid search cv 
  param_grid = dict(opt = opt, epochs = epochs, init = init, loss = loss, neurons = neurons, neurons_2 = neurons_2, batch_size = batch_size, lr = lr, l2 = l2, activation = activation)
  with tf.device('/cpu:0'):
    model = KerasClassifier(build_fn=create_model,batch_size=batch_size,epochs=epochs,verbose=1)
  clf = GridSearchCV(estimator=model, param_grid=param_grid, cv=5, n_jobs=1, 
                     return_train_score=True, scoring='accuracy', verbose = 2)
  grid_result = clf.fit(X_train,y_train)
  
  # write log file for neural network information
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  logfile = DIR_PATH+'nn_log_two_layer.txt'
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
  train_result.to_csv(DIR_PATH + 'train_predictions_two_layer.txt', sep = ' ', index = False)
  test_result.to_csv(DIR_PATH + 'test_predictions_two_layer.txt', sep = ' ', index = False)
  # save model
  best_model = clf.best_estimator_
  modelname = DIR_PATH + 'nn_classifier_two_layer.h5'
  best_model.model.save(modelname)
  
if __name__ == "__main__":
  main()
                 
  
  
  
  
  
  
  
