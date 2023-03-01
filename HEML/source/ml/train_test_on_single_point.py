import numpy as np
from HEML.utils.data import pull_mats_w_label, mat_pull
from HEML.utils.attrib import *
from HEML.utils.model import *
from HEML.utils.fields import pca, aug_all

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score


if __name__ == "__main__":

    aug = True
    pca_tf = True
    model = "xgb"

    # df = pd.read_csv("../../data/protein_data.csv")
    x, y = pull_mats_w_label(
        data_file="../../../data/protein_data.csv", dir_fields="../../../data/cpet/"
    )
    

    arr_min, arr_max, = np.min(x), np.max(x)
    # x = (x - arr_min) / (arr_max - arr_min + 1e-18)
    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)
    print(x.shape)
    y = [np.argmax(i) for i in y]
    
    # reduce x to just the center point of the matrix (midpoint_ind, midpoint_ind, midpoint_ind, 3)
    midpoint_ind = int(x.shape[1] / 2)
    x = x[:, midpoint_ind, midpoint_ind, midpoint_ind, :]
    # get magnitude of vector
    x = np.linalg.norm(x, axis=1).reshape(-1, 1)

    (
        X_train,
        X_test,
        y_train,
        y_test,
    ) = train_test_split(x, y, test_size=0.2, random_state=11)


    print(X_test.shape)
    
    model_obj = RandomForestClassifier()
    kf = KFold(n_splits=5, random_state=11, shuffle=True)
    acc_train, acc_val, f1_val, auroc_val = [], [], [], []

    for ind_train, ind_val in kf.split(X=X_train, y=y_train):
        x_train_temp = X_train[ind_train]
        x_val_temp = X_train[ind_val]
        y_train_temp = np.array(y_train)[ind_train].tolist()
        y_val_temp = np.array(y_train)[ind_val].tolist()

        model_obj.fit(x_train_temp, y_train_temp)
        y_train_pred = model_obj.predict(x_train_temp)
        acc_train.append(accuracy_score(y_train_pred, y_train_temp))

        y_val_pred = model_obj.predict(x_val_temp)
        acc_val.append(accuracy_score(y_val_temp, y_val_pred))
        f1_val.append(f1_score(y_val_temp, y_val_pred, average="weighted"))
        if model != "ridge":
            auroc_val.append(
                roc_auc_score(
                    y_val_temp, model_obj.predict_proba(x_val_temp), multi_class="ovr"
                )
            )

    y_test_pred = model_obj.predict(X_test)
    acc_test = accuracy_score(y_test, y_test_pred)
    f1_test = f1_score(y_test, y_test_pred, average="weighted")
    print(X_test.shape)

    print("train acc:        {:.2f}".format(np.mean(np.array(acc_train))))
    print("validation acc:   {:.2f}".format(np.mean(np.array(acc_val))))
    print("test acc:         {:.2f}".format(acc_test))
    print("validation f1:    {:.2f}".format(np.mean(np.array(f1_val))))
    print("test f1:          {:.2f}".format(f1_test))

    if model != "ridge":
        roc_auc_test = roc_auc_score(
            y_test, model_obj.predict_proba(X_test), multi_class="ovr"
        )
        print("validation auroc: {:.2f}".format(np.mean(np.array(auroc_val))))
        print("test auroc:       {:.2f}".format(roc_auc_test))
