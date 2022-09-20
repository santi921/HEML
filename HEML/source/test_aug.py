from HEML.utils.attrib_utils import *
from HEML.utils.data_utils import *
import pandas as pd 

df = pd.read_csv("../../data/protein_data.csv")
x, y = pull_mats_w_label('../../data/dat')
#x_aug_1, y_aug_1 = augment_mat(x[0], y[0])
x_aug, y_aug = aug_all(x, y)
print(np.shape(x_aug))