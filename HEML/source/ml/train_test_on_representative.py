
import numpy as np
from HEML.utils.data import *
from HEML.utils.attrib import *
from HEML.utils.model import *

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score


if __name__ == "__main__":

        aug = True
        pca_tf = True
        model = "xgb"
        
        #df = pd.read_csv("../../data/protein_data.csv")
        x, y = pull_mats_w_label(dir_data = "../../../data/protein_data.csv", dir_fields = "../../../data/cpet/")

        arr_min, arr_max,  = np.min(x), np.max(x)
        #x = (x - arr_min) / (arr_max - arr_min + 1e-18)
        x_sign = np.sign(x)
        # getting absolute value of every element
        x_abs = np.abs(x)
        # applying log1p
        x_log1p = np.log1p(x_abs)
        # getting sign back
        x = np.multiply(x_log1p, x_sign)
        y = [np.argmax(i) for i in y]
        
        (
            X_train, 
            X_test, 
            y_train, 
            y_test, 
        ) = train_test_split(x, y, test_size=0.2, random_state=11)

        if(pca_tf):
            X_train_untransformed = X_train
            X_test_untransformed = X_test
            
            _, pca_obj = pca(np.concatenate((X_train, X_test)), verbose = True, pca_comps=15) 
            X_train, pca_obj_train = pca(X_train, pca_obj)        
            X_test, pca_obj_test = pca(X_test, pca_obj)  
    
        print(X_test.shape)
        model_obj = XGBClassifier(
            eta=0.8885,
            gamma = 0.0756,
            max_depth = 3,
            subsample= 0.5897054280821705, 
            reg_lambda=0.0000011289136863384718, 
            alpha=0.00002,
            eval_metric = "mlogloss"
        )

        kf = KFold(n_splits=5, random_state = 11, shuffle = True)
        acc_train, acc_val, f1_val, auroc_val = [], [], [], []
        
        for ind_train, ind_val in kf.split(X=X_train, y=y_train): 
            x_train_temp = X_train[ind_train]
            x_val_temp = X_train[ind_val]
            y_train_temp = np.array(y_train)[ind_train].tolist()
            y_val_temp = np.array(y_train)[ind_val].tolist()
            
            if(aug):
                x_train_temp, y_train_temp = aug_all(
                    X_train_untransformed[ind_train], 
                    y_train_temp, 
                    xy = True, 
                    z = False
                    )
                x_train_temp, _ = pca(x_train_temp, pca_obj) 
                model_obj.fit(x_train_temp, y_train_temp) 
                y_train_pred = model_obj.predict(x_train_temp)
                acc_train.append(accuracy_score(y_train_pred, y_train_temp))

            else: 
                model_obj.fit(x_train_temp, y_train_temp) 
                y_train_pred = model_obj.predict(x_train_temp)
                acc_train.append(accuracy_score(y_train_pred, y_train_temp))

            y_val_pred = model_obj.predict(x_val_temp)
            acc_val.append(accuracy_score(y_val_temp, y_val_pred))
            f1_val.append(f1_score(y_val_temp, y_val_pred, average = "weighted"))
            if(model!="ridge"):
                auroc_val.append(roc_auc_score(y_val_temp, model_obj.predict_proba(x_val_temp), 
                        multi_class="ovr"))

        y_test_pred = model_obj.predict(X_test)    
        acc_test = accuracy_score(y_test, y_test_pred)
        f1_test = f1_score(y_test, y_test_pred, average = "weighted")
        print(X_test.shape)

        print("train acc:        {:.2f}".format(np.mean(np.array(acc_train))))
        print("validation acc:   {:.2f}".format(np.mean(np.array(acc_val))))
        print("test acc:         {:.2f}".format(acc_test))
        print("validation f1:    {:.2f}".format(np.mean(np.array(f1_val))))
        print("test f1:          {:.2f}".format(f1_test))
        
        if(model!="ridge"):
            roc_auc_test = roc_auc_score(y_test, model_obj.predict_proba(X_test), 
                        multi_class="ovr")            
            print("validation auroc: {:.2f}".format(np.mean(np.array(auroc_val))))
            print("test auroc:       {:.2f}".format(roc_auc_test))

        names, mat = [], []
        dir_fields = "../../../data/charges_md/cpet/"
        test_files = os.listdir(dir_fields)
        
        for file in test_files: 
            if("dat" in file):
                mat.append(mat_pull(dir_fields + file))
                names.append(file)
        x = np.array(mat)
        arr_min, arr_max,  = np.min(x), np.max(x)
        #x = (x - arr_min) / (arr_max - arr_min + 1e-18)
        x_sign = np.sign(x)
        # getting absolute value of every element
        x_abs = np.abs(x)
        # applying log1p
        x_log1p = np.log1p(x_abs)
        # getting sign back
        x = np.multiply(x_log1p, x_sign)
        print(x.shape)
            
        x, _ = pca(x, pca_obj)
        print(x.shape)
        print(x.shape[0])
        result = {}
        for i in range(x.shape[0]):
            name_pro = names[i].split("_")[3]
            if(name_pro not in result):
                result[name_pro] = [0, 0, 0]
            pred = model_obj.predict(x[i].reshape(1, -1))

            if(pred[0] == 0):
                result[name_pro][0] += 1
            elif(pred[0] == 1):
                result[name_pro][1] += 1
            else:
                result[name_pro][2] += 1
        print(result)
        