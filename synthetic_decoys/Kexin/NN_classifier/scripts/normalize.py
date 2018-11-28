import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
def norm_errors(errors):
  cols = errors.columns
  cols_1 = ['id', 'model', 'flag', 'rmsd']
  cols_2 = [e for e in cols if e not in cols_1]
  
  # standard scaling
  errors_1 = errors[cols_1]
  errors_2 = errors[cols_2]
  scaler = StandardScaler()
  errors_2 = scaler.fit_transform(errors_2)
  return pd.DataFrame(np.hstack([errors_1, errors_2]))
  
  
  
def main():
  DIR_PATH = '/home/kexin/projects/global_quality_assessment/'  
  # load error matrix
  errors = pd.read_table(DIR_PATH+'errors.txt',sep=' ')
  
  # fill missing values with mean column values
  # errors.fillna(errors.mean(), inplace=True)
  # fill missing values with 0
  errors.fillna(0, inplace=True)
  
  # normalize the error for each RNA (specificed by id) separately
  cols = errors.columns
  errors = errors.groupby('id').apply(norm_errors)
  errors.columns = cols
  errors.to_csv(DIR_PATH + 'errors_normalized.txt', sep = ' ', index = False)

if __name__ == "__main__":
  main()  