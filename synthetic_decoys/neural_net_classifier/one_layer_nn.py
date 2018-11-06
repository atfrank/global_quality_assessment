import pandas as pd
import numpy as np
import sys
import timeit

from sklearn.preprocessing import StandardScaler
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
norm_errors <- function(errors){
  # normalizes the errors based on the mean
  errors$total <- rowSums(errors[,!(colnames(errors) %in% c("id", "model", "reference_flag" ,"rmsd"))], na.rm = TRUE)
  cols <- colnames(errors)[!(colnames(errors) %in% c("id", "model", "reference_flag","rmsd"))]
  for (col in cols){
    #errors[, col] <- errors[, col]/median(errors[, col])
    errors[, col] <- errors[, col]/median(errors[, col])
  }
  return(errors)
}
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
  

def create_model(neurons, activation, init, constraint, rr, lr, loss, opt):
  model = Sequential()
  model.add(Dense(neurons, input_dim = 76, activation = activation, 
  kernel_initializer = init, kernel_constraint = maxnorm(constraint), 
  kernel_regularizer = regularizers.l2(rr)))
  model.add(Dense(1, activation = 'sigmoid'))
  if opt == 'Adam':
    optimizer = Adam(lr = lr)
  if opt == 'RMSprop':
    optimizer = RMSprop(lr = lr)
  model.compile(loss = loss, optimizer = optimizer, metrics = ['accuracy'])
  return model

def main():
  DIR_PATH = '~/documents/github/global_quality_assessment/synthetic_decoys/neural_net_classifier/'
  # load errors matrices
  errors = pd.read_table(DIR_PATH+'errors.txt',sep=' ')
  
  # fill missing values with mean column values
  # errors.fillna(errors.mean(), inplace=True)
  # fill missing values with 0
  errors.fillna(0, inplace=True)
  
  # normalize the error for each RNA (specificed by id) separately
  cols = errors.columns
  errors = errors.groupby('id').apply(norm_errors)
  #errors.columns = cols
  
  # scramble data just to be sure
  errors = errors.sample(len(errors),random_state = 1)
  
  # SPLIT INTO TRAINING AND TESTING
  
  # BUILD CLASSIFIER
  