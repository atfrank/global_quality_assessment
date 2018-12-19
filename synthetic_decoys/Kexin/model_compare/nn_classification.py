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

def create_flag(row, m):
  if row['DI'] > m:
    return 1
  else:
    return 0

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
  DIR_PATH = '/home/kexin/projects/global_quality_assessment/RESTORE/'
  # LOAD NORMALIZED ERROR MATRIX
  errors = pd.read_table(DIR_PATH+'errors_with_DI.txt',sep=' ')
  errors = errors.sample(len(errors))
  errors = errors.drop(columns=['flag'])
  m = float(sys.argv[2])
  errors['flag'] = errors.apply(lambda row: create_flag(row, m), axis = 1) # self defined threshold

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
  rna_id = sys.argv[1]
  rna = rnas[int(rna_id)]
  print('target RNA: ', rna)
  fil = sim[rna]>threshold
  remove_rnas = list(compress(rnas, fil.values))
  train_rnas = [e for e in rnas if e not in remove_rnas]
  print('sequences removed: ', remove_rnas)

  # SPLIT TRAIN TEST
  train = errors.loc[errors['id'].isin(train_rnas)]
  test = errors[errors['id']==rna]
  # scramble data just to be sure
  train = train.sample(len(train))
  test = test.sample(len(test))
  test_result = test[['id','model','flag','DI','rmsd']]
  X_train = train.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  X_test = test.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  y_train = train['flag']
  y_test = test['flag']
  del train, test, errors

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
  grid = dict(opt = opt, epochs = epochs, init = init, loss = loss, neurons = neurons, batch_size = batch_size, lr = lr, l2 = l2, activation = activation)
  with tf.device('/cpu:0'):
    model = KerasClassifier(build_fn=create_model,verbose=2)
  clf_cv = GridSearchCV(estimator=model, param_grid=grid, cv=5, n_jobs=1,
                     return_train_score=True, scoring='accuracy', verbose = 2)
  grid_result = clf_cv.fit(X_train,y_train)

  # write log file for neural network information
  logfile = DIR_PATH+'model_compare/log/nnclf_'+rna+'_'+str(m)+'.txt'
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  with open(logfile, 'a') as f:
    print("classification threshold is ", m)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  best_model = clf_cv.best_estimator_

  # PREDICT
  pred_binary = best_model.predict(X_test).flatten().tolist()
  pred_proba = best_model.predict_proba(X_test).flatten().tolist()
  test_result['pred_binary'] = pred_binary
  test_result['pred_proba'] = pred_proba
  # print confusion matrix
  with open(logfile, 'a') as f:
    print(classification_report(y_test, y_pred_test), file=f)
    print(confusion_matrix(y_test, y_pred_test), file=f)
    print("score is: %f" % best_model.score(X_test,y_test), file=f)
  test_result.to_csv(DIR_PATH + 'model_compare/pred/nnclf_'+rna+'_'+str(m)+'.txt', sep = ' ', index = False)

if __name__ == "__main__":
  main()
