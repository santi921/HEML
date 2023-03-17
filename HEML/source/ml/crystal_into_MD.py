import os 
import numpy as np
from imblearn.over_sampling import BorderlineSMOTE
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score


from HEML.utils.attrib import *
from HEML.utils.model import *
#from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier

from HEML.utils.fields import pca, aug_all
from HEML.utils.data import pull_mats_w_label, mat_pull
from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier


if __name__ == "__main__":

    aug = True
    pca_tf = True
    model = "xgb"

    # df = pd.read_csv("../../data/protein_data.csv")
    x, y = pull_mats_w_label(
        data_file="../../../data/protein_data.csv", dir_fields="../../../data/cpet/"
    )

    arr_min, arr_max, = np.min(
        x
    ), np.max(x)
    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)
    #x = (x - arr_min) / (arr_max - arr_min + 1e-18)

    y = [np.argmax(i) for i in y]

    (
        X_train,
        X_test,
        y_train,
        y_test,
    ) = train_test_split(x, y, test_size=0.2, random_state=11)

    if pca_tf:
        X_train_untransformed = X_train
        X_test_untransformed = X_test

        sm = BorderlineSMOTE(random_state=42)
        shape_before = X_train_untransformed.shape
        X_train_untransformed = X_train_untransformed.reshape(X_train_untransformed.shape[0], -1)
        X_train_untransformed, y_train = sm.fit_resample(X_train_untransformed, y_train)
        X_train_untransformed = X_train_untransformed.reshape(len(X_train_untransformed), shape_before[1], shape_before[2], shape_before[3], shape_before[4])
        _, pca_obj = pca(np.concatenate((X_train, X_train_untransformed)), verbose=True, pca_comps=20)
        X_train, pca_obj_train = pca(X_train_untransformed, pca_obj)
        X_test, pca_obj_test = pca(X_test, pca_obj)
        
        
        
    print(X_train.shape)
    print(X_test.shape)
    
    model_obj = XGBClassifier(
        eta=0.78,
        gamma=0.6,
        max_depth=7,
        subsample=0.80,
        reg_lambda=0.0000011289136863384718,
        alpha=0.0001,
        eval_metric="mlogloss",
    )
    #model_obj = BalancedRandomForestClassifier(
    #    n_estimators=500,
    #    max_depth=7,
    #)
    #model = RandomForestClassifier(
    #    n_estimators=700, 
    #    max_depth=7
    #)

    #model_obj = EasyEnsembleClassifier(
    #)
    kf = KFold(n_splits=5, random_state=11, shuffle=True)
    acc_train, acc_val, f1_val, auroc_val = [], [], [], []

    for ind_train, ind_val in kf.split(X=X_train, y=y_train):
        x_train_temp = X_train[ind_train]
        x_val_temp = X_train[ind_val]
        y_train_temp = np.array(y_train)[ind_train].tolist()
        y_val_temp = np.array(y_train)[ind_val].tolist()

        if aug:
            x_train_temp, y_train_temp = aug_all(
                X_train_untransformed[ind_train], y_train_temp, xy=True, z=False
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

    names, mat = [], []
    dir_fields = "../../../data/fields/"
    #dir_fields = "../../../data/charges_md/cpet/"
    #dir_fields = "../../../data/md_run_boxes/1ebe/"
    test_files = os.listdir(dir_fields)

    for file in test_files:
        if "dat" in file:
            mat.append(mat_pull(dir_fields + file))
            names.append(file)
    x = np.array(mat)

    # x = (x - arr_min) / (arr_max - arr_min + 1e-18)
    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)
    #x = (x - arr_min) / (arr_max - arr_min + 1e-18)

    print(x.shape)

    x, _ = pca(x, pca_obj)
    print(x.shape)
    print(x.shape[0])
    result = {}
    
    for i in range(x.shape[0]):
        print(names[i])
        name_pro = names[i].split(".")[0].split("_")[3]
        
        if name_pro not in result:
            result[name_pro] = [0, 0, 0]
        pred = model_obj.predict(x[i].reshape(1, -1))

        if pred[0] == 0:
            result[name_pro][0] += 1
        elif pred[0] == 1:
            result[name_pro][1] += 1
        else:
            result[name_pro][2] += 1
    print(result)
