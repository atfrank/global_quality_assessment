import pandas as pd
import numpy as np
import sys
import timeit
from itertools import compress
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report,confusion_matrix,accuracy_score,recall_score
from sklearn.linear_model import LogisticRegression
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
  # FEATURE SELECTION
  '''
  X = errors.drop(columns=['id','model','flag','rmsd','inf_all','DI'])
  y = errors['flag']
  logit_model=sm.Logit(y,X)
  result=logit_model.fit()
  print(result.summary2())

                           Results: Logit
================================================================
Model:              Logit            No. Iterations:   7.0000
Dependent Variable: flag             Pseudo R-squared: 0.295
Date:               2018-12-17 13:14 AIC:              6565.8423
No. Observations:   6700             BIC:              6960.8144
Df Model:           57               Log-Likelihood:   -3224.9
Df Residuals:       6642             LL-Null:          -4577.4
Converged:          1.0000           Scale:            1.0000
-----------------------------------------------------------------
             Coef.   Std.Err.     z      P>|z|    [0.025   0.975]
-----------------------------------------------------------------
C1pGUA      -0.3068    0.0483   -6.3493  0.0000  -0.4015  -0.2121
C2pGUA      -0.1833    0.0530   -3.4604  0.0005  -0.2871  -0.0795
C3pGUA      -0.2638    0.0485   -5.4356  0.0000  -0.3589  -0.1687
C4pGUA      -0.2810    0.0446   -6.2986  0.0000  -0.3684  -0.1935
C5pGUA      -0.1505    0.0425   -3.5435  0.0004  -0.2337  -0.0673
C8GUA       -0.2685    0.0509   -5.2793  0.0000  -0.3682  -0.1688
H1pGUA      -0.3344    0.0550   -6.0830  0.0000  -0.4422  -0.2267
H2pGUA      -0.3640    0.0567   -6.4236  0.0000  -0.4751  -0.2529
H3pGUA      -0.0644    0.0473   -1.3630  0.1729  -0.1571   0.0282
H4pGUA      -0.0415    0.0530   -0.7841  0.4330  -0.1454   0.0623
H5pGUA      -0.0720    0.0541   -1.3314  0.1830  -0.1781   0.0340
H5ppGUA     -0.0669    0.0544   -1.2302  0.2186  -0.1735   0.0397
H2GUA       -0.0620    0.0469   -1.3218  0.1862  -0.1540   0.0300
H5GUA       -0.2012    0.0536   -3.7519  0.0002  -0.3063  -0.0961
H6GUA       -0.0423    0.0584   -0.7238  0.4692  -0.1567   0.0722
H8GUA       -0.1948    0.0634   -3.0738  0.0021  -0.3190  -0.0706
C1pADE      -0.0957    0.0481   -1.9918  0.0464  -0.1899  -0.0015
C3pADE       0.0460    0.0492    0.9357  0.3494  -0.0504   0.1424
C5pADE      -0.2049    0.0454   -4.5099  0.0000  -0.2939  -0.1158
C2ADE        0.0004    0.0568    0.0062  0.9951  -0.1110   0.1118
C5ADE       -0.2128    0.0642   -3.3139  0.0009  -0.3387  -0.0870
C6ADE       -0.0617    0.0628   -0.9835  0.3254  -0.1848   0.0613
C8ADE       -0.1129    0.0522   -2.1633  0.0305  -0.2151  -0.0106
H2pADE      -0.3545    0.0509   -6.9645  0.0000  -0.4542  -0.2547
H4pADE      -0.2451    0.0440   -5.5734  0.0000  -0.3313  -0.1589
H5pADE      -0.2200    0.0476   -4.6201  0.0000  -0.3133  -0.1267
H2ADE       -0.1591    0.0485   -3.2816  0.0010  -0.2542  -0.0641
H6ADE        0.1594    0.0370    4.3016  0.0000   0.0868   0.2320
H8ADE       -0.1024    0.0352   -2.9071  0.0036  -0.1715  -0.0334
C1pCYT      -0.2046    0.0367   -5.5731  0.0000  -0.2765  -0.1326
C2pCYT      -0.0562    0.0339   -1.6567  0.0976  -0.1227   0.0103
C3pCYT      -0.1712    0.0360   -4.7577  0.0000  -0.2418  -0.1007
C5CYT       -0.0508    0.0369   -1.3756  0.1689  -0.1231   0.0216
C6CYT       -0.1303    0.0362   -3.5945  0.0003  -0.2013  -0.0592
C8CYT        0.1029    0.0360    2.8598  0.0042   0.0324   0.1734
H1pCYT      -0.0543    0.0358   -1.5153  0.1297  -0.1245   0.0159
H2pCYT       0.3007    0.0375    8.0194  0.0000   0.2272   0.3742
H3pCYT      -0.1950    0.0378   -5.1610  0.0000  -0.2690  -0.1209
H4pCYT       0.0558    0.0356    1.5681  0.1168  -0.0139   0.1255
H5pCYT       0.1028    0.0365    2.8195  0.0048   0.0313   0.1742
H5ppCYT     -0.1457    0.0414   -3.5162  0.0004  -0.2269  -0.0645
H2CYT       -0.1461    0.0429   -3.4063  0.0007  -0.2301  -0.0620
H5CYT       -0.1550    0.0407   -3.8056  0.0001  -0.2348  -0.0752
H6CYT       -0.0265    0.0384   -0.6912  0.4894  -0.1018   0.0487
C1pURA      -0.2634    0.0399   -6.5936  0.0000  -0.3416  -0.1851
C3pURA      -0.2591    0.0383   -6.7707  0.0000  -0.3341  -0.1841
C4pURA      -0.1193    0.0464   -2.5712  0.0101  -0.2103  -0.0284
C5pURA       0.0016    0.0416    0.0392  0.9687  -0.0800   0.0832
C2URA        0.0084    0.0462    0.1810  0.8563  -0.0821   0.0988
C5URA       -0.0971    0.0416   -2.3321  0.0197  -0.1787  -0.0155
C6URA       -0.0246    0.0485   -0.5064  0.6126  -0.1196   0.0705
C8URA       -0.0873    0.0462   -1.8889  0.0589  -0.1779   0.0033
H1pURA       0.1726    0.0498    3.4639  0.0005   0.0749   0.2702
H2pURA       0.1891    0.0444    4.2592  0.0000   0.1021   0.2761
H4pURA      -0.3039    0.0384   -7.9233  0.0000  -0.3791  -0.2287
H5ppURA     -0.1445    0.0370   -3.9096  0.0001  -0.2169  -0.0721
H2URA       -0.2864    0.0367   -7.8097  0.0000  -0.3582  -0.2145
H6URA       -0.5433    0.0414  -13.1127  0.0000  -0.6245  -0.4621
================================================================
  '''
  cols_remove = ['H3pGUA', 'H4pGUA', 'H5pGUA', 'H5ppGUA', 'H2GUA', 'H6GUA', 'C3pADE',
                 'C2ADE', 'C6ADE', 'C2pCYT', 'C5CYT', 'H1pCYT', 'H4pCYT', 'H6CYT',
                 'C5pURA', 'C2URA', 'C6URA', 'C8URA']
  errors = errors.drop(columns = cols_remove)


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
  logreg = LogisticRegression()
  grid = {"C":np.logspace(-3, 3, 7), "penalty": ["l1","l2"]}
  logreg_cv = GridSearchCV(logreg, grid, cv = 5, verbose = 2, scoring = 'accuracy')
  grid_result = logreg_cv.fit(X_train, y_train)
  logfile = DIR_PATH+'model_compare/log/logit_'+rna+'_'+str(m)+'.txt'
  means = grid_result.cv_results_['mean_test_score']
  stds = grid_result.cv_results_['std_test_score']
  params = grid_result.cv_results_['params']
  with open(logfile, 'a') as f:
    print("classification threshold is ", m)
    print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_), file=f)
    for mean, stdev, param in zip(means, stds, params):
      print("%f (%f) with: %r" % (mean, stdev, param), file=f)
  best_model = logreg_cv.best_estimator_
  # PREDICT
  pred_binary = best_model.predict(X_test)
  pred_proba = best_model.predict_proba(X_test)[:,1]
  test_result.loc[:,'pred_binary'] = pred_binary
  test_result.loc[:,'pred_proba'] = pred_proba
  with open(logfile, 'a') as f:
    print(classification_report(y_test, pred_binary), file=f)
    print(confusion_matrix(y_test, pred_binary), file=f)
  test_result.to_csv(DIR_PATH + 'model_compare/pred/logit_'+rna+'_'+str(m)+'.txt', sep = ' ', index = False)
if __name__ == "__main__":
  main()
