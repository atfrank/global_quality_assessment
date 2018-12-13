import pandas as pd
import numpy as np
import sys
import timeit
import pickle
from itertools import compress
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,recall_score
from sklearn.ensemble import RandomForestRegressor


def main():
  DIR_PATH = '/home/kexin/projects/global_quality_assessment/RESTORE/'  
  # LOAD NORMALIZED ERROR MATRIX
  errors = pd.read_table(DIR_PATH+'errors_with_DI.txt',sep=' ')
  
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
  
  train = errors.loc[errors['id'].isin(train_rnas)]
  test = errors[errors['id']==rna]
  
  # scramble data just to be sure
  train = train.sample(len(train),random_state=1)
  test = test.sample(len(test),random_state=1)
  
  # make copy of train and test actual flag
  train_result = train[['id','model','flag']]
  test_result = test[['id','model','flag']]
  X_train = train.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  X_test = test.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  y_train = train['DI']
  y_test = test['DI']
  
  del train, test
  
  # BUILD REGRESSOR
  # grid search parameters
  n_estimators = [200, 400, 600, 800, 1000]
  min_samples_split = [2, 5, 10]
  bootstrap = [True]
  max_features = ['sqrt', 'log2', None, ]
 
  # grid search cv 
  rf = RandomForestRegressor(oob_score = True, verbose = 1)
  param_grid = dict(n_estimators = n_estimators, min_samples_split = min_samples_split, bootstrap = bootstrap, max_features = max_features)
  clf = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, n_jobs=4, 
                     return_train_score=True, scoring='neg_mean_absolute_error', verbose = 2)
  grid_result = clf.fit(X_train,y_train)
  
  # write log file for neural network information
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  logfile = DIR_PATH+'tree_based/log/DI_rf_rg_log_'+rna+'.txt'
  with open(logfile, 'a') as f:
    print("OOB score: %f" % (clf.best_estimator_.oob_score_), file=f)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  
  # predict and format results
  y_pred_train = clf.predict(X_train).flatten().tolist()
  train_result.loc[:,'pred'] = y_pred_train
  y_pred_test = clf.predict(X_test).flatten().tolist()
  test_result.loc[:,'pred'] = y_pred_test
  # print confusion matrix
  #with open(logfile, 'a') as f:
    #print(classification_report(y_test, y_pred_test), file=f)
    #print(confusion_matrix(y_test, y_pred_test), file=f)
  train_result.to_csv(DIR_PATH + 'tree_based/pred/DI_rf_rg_train_pred_'+rna+'.txt', sep = ' ', index = False)
  test_result.to_csv(DIR_PATH + 'tree_based/pred/DI_rf_rg_test_pred_'+rna+'.txt', sep = ' ', index = False)
  
  # save best model
  best_model = clf.best_estimator_
  modelname = DIR_PATH + 'models/DI_rf_rg_model_'+rna+'.sav'
  pickle.dump(best_model, open(modelname, 'wb'))
  
if __name__ == "__main__":
  main()
                 
  
  
  
  
  
  
  
