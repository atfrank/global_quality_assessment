import pandas as pd
import numpy as np
import sys
import timeit
from itertools import compress
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,recall_score
from sklearn.ensemble import ExtraTreesClassifier
#import statsmodels.api as sm

def create_flag(row, m):
  if row['DI'] > m:
    return 1
  else:
    return 0

def main():
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
  test_result = test[['id', 'model', 'flag', 'DI','rmsd']]
  X_train = train.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  X_test = test.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  y_train = train['flag']
  y_test = test['flag']
  del train, test, errors

  # IMPLEMENT MODEL
  ert = ExtraTreesClassifier(oob_score = True, verbose = 2)
  # grid search parameters
  n_estimators = [800, 1000, 1200]
  min_samples_split = [2, 5]
  max_features = ['sqrt', 'log2', None]
  bootstrap = ['True']
  class_weight = ['balanced_subsample', 'balanced', None]
  grid = dict(n_estimators = n_estimators, min_samples_split = min_samples_split,
  bootstrap = bootstrap, max_features = max_features, class_weight = class_weight)
  ert_cv = GridSearchCV(ert, grid, cv = 5, verbose = 2, scoring = 'accuracy', n_jobs = 6)
  grid_result = ert_cv.fit(X_train, y_train)

  logfile = DIR_PATH+'model_compare/log/ertclf_'+rna+'_'+str(m)+'.txt'
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  with open(logfile, 'a') as f:
    print("classification threshold is ", m)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  best_model = ert_cv.best_estimator_
  # PREDICT
  pred_binary = best_model.predict(X_test)
  pred_proba = best_model.predict_proba(X_test)[:,1]
  test_result.loc[:,'pred_binary'] = pred_binary
  test_result.loc[:,'pred_proba'] = pred_proba
  with open(logfile, 'a') as f:
    print(classification_report(y_test, pred_binary), file=f)
    print(confusion_matrix(y_test, pred_binary), file=f)
    print("R^2 value is: %f" % best_model.score(X_test,y_test), file=f)
  test_result.to_csv(DIR_PATH + 'model_compare/pred/ertclf_'+rna+'_'+str(m)+'.txt', sep = ' ', index = False)

if __name__ == "__main__":
  main()
