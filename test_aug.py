from attrib_utils import *
from data_utils import *
import pandas as pd 

df = pd.read_csv("protein_data.csv")
x, y = pull_mats_w_label('./dat')
#x_aug_1, y_aug_1 = augment_mat(x[0], y[0])
x_aug, y_aug = aug_all(x, y)
print(np.shape(x_aug))