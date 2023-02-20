
from operator import concat
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import wandb, argparse
from HEML.utils.data import *
from HEML.utils.attrib import *
from HEML.utils.model import *

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

class training: 
    def __init__(self, model, pca_tf=True, aug=True):
        self.aug = aug
        self.pca_tf = pca_tf
        self.model = model
        
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
            self.X_train, 
            self.X_test, 
            self.y_train, 
            self.y_test, 
        ) = train_test_split(x, y, test_size=0.2, random_state=11)

        if(pca_tf):
            self.X_train_untransformed = self.X_train
            self.X_test_untransformed = self.X_test
            
            _, self.pca_obj = pca(np.concatenate((self.X_train, self.X_test)), verbose = True, pca_comps=15) 
            self.X_train, self.pca_obj_train = pca(self.X_train, self.pca_obj)        
            self.X_test, self.pca_obj_test = pca(self.X_test, self.pca_obj)  


    def make_model(self, config):
        model_obj = construct_models(config = config, model = self.model)
        return model_obj


    def train(self):
        with wandb.init(project="HemeML_redux") as run:
            config = wandb.config
            model_obj = self.make_model(config)

            kf = KFold(n_splits=5, random_state = 11, shuffle = True)
            acc_train, acc_val, f1_val, auroc_val = [], [], [], []
            
            for ind_train, ind_val in kf.split(X=self.X_train, y=self.y_train): 
                x_train_temp = self.X_train[ind_train]
                x_val_temp = self.X_train[ind_val]
                y_train_temp = np.array(self.y_train)[ind_train].tolist()
                y_val_temp = np.array(self.y_train)[ind_val].tolist()
                
                if(self.aug):
                    x_train_temp, y_train_temp = aug_all(
                        self.X_train_untransformed[ind_train], 
                        y_train_temp, 
                        xy = True, 
                        z = False
                        )
                    x_train_temp, _ = pca(x_train_temp, self.pca_obj) 
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

            y_test_pred = model_obj.predict(self.X_test)    
            acc_test = accuracy_score(self.y_test, y_test_pred)
            f1_test = f1_score(self.y_test, y_test_pred, average = "weighted")


            print("train acc:        {:.2f}".format(np.mean(np.array(acc_train))))
            print("validation acc:   {:.2f}".format(np.mean(np.array(acc_val))))
            print("test acc:         {:.2f}".format(acc_test))
            print("validation f1:    {:.2f}".format(np.mean(np.array(f1_val))))
            print("test f1:          {:.2f}".format(f1_test))
            wandb.log({"train_acc": np.mean(np.array(acc_train))})
            wandb.log({"val_acc": np.mean(np.array(acc_val))})
            wandb.log({"test_acc": np.mean(np.array(acc_test))})
            wandb.log({"val_f1": np.mean(np.array(f1_val))})
            wandb.log({"test_f1": f1_test})
            if(model!="ridge"):
                roc_auc_test = roc_auc_score(self.y_test, model_obj.predict_proba(self.X_test), 
                            multi_class="ovr")            
                print("validation auroc: {:.2f}".format(np.mean(np.array(auroc_val))))
                print("test auroc:       {:.2f}".format(roc_auc_test))
                wandb.log({"val_auroc": np.mean(np.array(auroc_val))})
                wandb.log({"test_auroc": roc_auc_test})

        run.finish()


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='options for hyperparam tune')
    parser.add_argument(
        "-dataset",
        action="store",
        dest="dataset",
        default=1,
        help="dataset to use",
    )

    parser.add_argument(
        "-count",
        action="store",
        dest="count",
        default=100,
        help="number of hyperparams",
    )

    parser.add_argument(
        "-model",
        action="store",
        dest="model",
        default="rfc",
        help="model",
    )

    parser.add_argument(
        "--aug", 
        action = "store_true",
        dest="aug",
        default = False,
        help="augment data"
    )
    
    pca_tf = True
    method = "bayes" 
    dataset_name = "base"
    results = parser.parse_args()
    model = str(results.model)
    dataset_int = int(results.dataset)
    count = int(results.count)
    aug = bool(results.aug)
    
    sweep_config = {}
    dict_hyper = hyperparameter_dicts()
    sweep_config["parameters"] = dict_hyper[model]
    sweep_config["name"] = method + "_" + model + "_" + dataset_name
    sweep_config["method"] = method 
    if(method == "bayes"):
        sweep_config["metric"] = {"name": "val_f1", "goal": "maximize"}
    sweep_id = wandb.sweep(sweep_config, project = "HemeML")
    training_obj = training(model = model, pca_tf=pca_tf)
    wandb.agent(sweep_id, function=training_obj.train, count=count)
    

    



