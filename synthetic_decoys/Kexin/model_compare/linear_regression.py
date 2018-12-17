import pandas as pd
import numpy as np
import sys
import timeit
from itertools import compress
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,recall_score
from sklearn.linear_model import LinearRegression
#import statsmodels.api as sm

def main():
  DIR_PATH = '/home/kexin/projects/global_quality_assessment/RESTORE/'
  # LOAD NORMALIZED ERROR MATRIX
  errors = pd.read_table(DIR_PATH+'errors_with_DI.txt',sep=' ')
  errors = errors.sample(len(errors))
  errors = errors.drop(columns=['flag'])

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
  train = train.sample(len(train))
  test = test.sample(len(test))
  test_result = test[['id', 'model', 'DI','rmsd']]
  X_train = train.drop(columns=['id','model','rmsd','inf_all','DI'])
  X_test = test.drop(columns=['id','model','rmsd','inf_all','DI'])
  y_train = train['DI']
  y_test = test['DI']
  del train, test, errors

  # IMPLEMENT MODEL
  lm = LinearRegression(fit_intercept = True)
  grid = {"normalize":["True","False"]}
  lm_cv = GridSearchCV(lm, grid, cv = 5, verbose = 2, scoring = 'neg_mean_absolute_error')
  grid_result = lm_cv.fit(X_train, y_train)
  logfile = DIR_PATH+'model_compare/log/lm_'+rna+'.txt'
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  with open(logfile, 'a') as f:
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  best_model = lm_cv.best_estimator_
  # PREDICT
  pred = best_model.predict(X_test)
  test_result.loc[:,'pred'] = pred
  with open(logfile, 'a') as f:
    print("R^2 is: %f" % best_model.score(X_test, y_test), file = f)
  test_result.to_csv(DIR_PATH + 'model_compare/pred/lm_'+rna+'.txt', sep = ' ', index = False)
if __name__ == "__main__":
  main()
