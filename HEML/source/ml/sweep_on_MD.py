from operator import concat
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
import wandb, argparse
from HEML.utils.data import pull_mats_from_MD_folder, pull_mats_w_label
from HEML.utils.attrib import *
from HEML.utils.model import *
from HEML.utils.fields import pca, aug_all

import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score, roc_auc_score

"""class config:
    def __init__(
        self, 
        nestimators,
        max_depth,
    ):
        self.nestimators = nestimators
        self.max_depth = max_depth
        self.bootstrap = True"""


class training:
    def __init__(self, model, pca_tf=True, aug=True, test_crystal=False, test_md=False):
        self.aug = aug
        self.pca_tf = pca_tf
        self.model = model
        self.test_crystal = test_crystal
        self.test_md = test_md

        # df = pd.read_csv("../../data/protein_data.csv")
        x, y, _ = pull_mats_from_MD_folder(label_ind=3)

        (
            arr_min,
            arr_max,
        ) = np.min(
            x
        ), np.max(x)

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
        ) = train_test_split(x, y, test_size=0.2, random_state=0)
        # print(self.y_test.shape)
        if self.pca_tf:
            self.X_train_untransformed = self.X_train
            self.X_test_untransformed = self.X_test

            _, self.pca_obj = pca(
                np.concatenate((self.X_train, self.X_test)), verbose=True, pca_comps=25
            )
            self.X_train, self.pca_obj_train = pca(self.X_train, self.pca_obj)
            self.X_test, self.pca_obj_test = pca(self.X_test, self.pca_obj)

        if self.test_md:
            x_md_test, y_md_test, _ = pull_mats_from_MD_folder(
                root_dir="../../../data/fields_test/",
                data_file="../../../data/protein_data.csv",
                label_ind=3,
            )
            x_sign = np.sign(x_md_test)
            # getting absolute value of every element
            x_abs = np.abs(x_md_test)
            # applying log1p
            x_log1p = np.log1p(x_abs)
            # getting sign back
            x_md_test = np.multiply(x_log1p, x_sign)
            # print('pull mats on md test {} {}'.format(len(x_md_test), len(y_md_test)))
            self.x_md_test, _ = pca(x_md_test, self.pca_obj)
            self.y_md_test = [np.argmax(i) for i in y_md_test]

        if test_crystal:
            x_crystal, y_crystal = pull_mats_w_label(
                data_file="../../../data/protein_data.csv",
                dir_fields="../../../data/cpet/",
            )

            x_sign = np.sign(x_crystal)
            x_abs = np.abs(x_crystal)
            x_log1p = np.log1p(x_abs)
            x_crystal = np.multiply(x_log1p, x_sign)

            self.x_crystal, _ = pca(x_crystal, self.pca_obj)
            self.y_crystal = [np.argmax(i) for i in y_crystal]

    def make_model(self, config):
        model_obj = construct_models(config=config, model=self.model)
        return model_obj

    def train(self):
        with wandb.init(project="HemeML_MD") as run:
            config = wandb.config

            model_obj = self.make_model(config)

            kf = KFold(n_splits=5, random_state=0, shuffle=True)
            acc_train, acc_val, f1_val, auroc_val = [], [], [], []

            print("kfold testing...")
            track_kfold = 0
            for ind_train, ind_val in kf.split(X=self.X_train, y=self.y_train):
                print("track kfold: " + str(track_kfold))
                track_kfold += 1
                x_train_temp = self.X_train[ind_train]
                x_val_temp = self.X_train[ind_val]
                y_train_temp = np.array(self.y_train)[ind_train].tolist()
                y_val_temp = np.array(self.y_train)[ind_val].tolist()

                if self.aug:
                    # print("aug...")
                    x_train_temp, y_train_temp = aug_all(
                        self.X_train_untransformed[ind_train],
                        y_train_temp,
                        xy=True,
                        z=False,
                    )
                    # print(x_train_temp.shape)
                    # print(len(y_train_temp))
                    # print(y_train_temp[0])
                    x_train_temp, _ = pca(x_train_temp, self.pca_obj)
                    model_obj.fit(x_train_temp, y_train_temp)
                    # print("done  fitting model")

                    y_train_pred = model_obj.predict(x_train_temp)
                    # print("done predicting")
                    acc_train.append(accuracy_score(y_train_pred, y_train_temp))
                    # print("adding accuracy")
                else:
                    model_obj.fit(x_train_temp, y_train_temp)
                    y_train_pred = model_obj.predict(x_train_temp)
                    acc_train.append(accuracy_score(y_train_pred, y_train_temp))

                y_val_pred = model_obj.predict(x_val_temp)
                # print("val predicted")
                acc_val.append(accuracy_score(y_val_temp, y_val_pred))
                # print("val accuraccy")
                f1_val.append(f1_score(y_val_temp, y_val_pred, average="weighted"))
                # print("val f1")
                # if model != "ridge":
                #    auroc_val.append(
                #        roc_auc_score(
                #            y_val_temp,
                #            model_obj.predict_proba(x_val_temp),
                #            multi_class="ovr",
                #        )
                #    )
            print("final testing...")
            y_test_pred = model_obj.predict(self.X_test)
            acc_test = accuracy_score(self.y_test, y_test_pred)
            f1_test = f1_score(self.y_test, y_test_pred, average="weighted")

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

            """ if model != "ridge":
                roc_auc_test = roc_auc_score(
                    self.y_test, model_obj.predict_proba(self.X_test), multi_class="ovr"
                )
                print("validation auroc: {:.2f}".format(np.mean(np.array(auroc_val))))
                print("test auroc:       {:.2f}".format(roc_auc_test))
                wandb.log({"val_auroc": np.mean(np.array(auroc_val))})
                wandb.log({"test_auroc": roc_auc_test})"""

            if self.test_crystal:
                print("test crystal structures...")
                print(
                    "crystal acc:      {:.2f}".format(
                        accuracy_score(
                            self.y_crystal, model_obj.predict(self.x_crystal)
                        )
                    )
                )
                wandb.log(
                    {
                        "crystal_acc": accuracy_score(
                            self.y_crystal, model_obj.predict(self.x_crystal)
                        )
                    }
                )
                print(
                    "crystal f1 score: {:.2f}".format(
                        f1_score(
                            self.y_crystal,
                            model_obj.predict(self.x_crystal),
                            average="weighted",
                        )
                    )
                )
                wandb.log(
                    {
                        "crystal_f1": f1_score(
                            self.y_crystal,
                            model_obj.predict(self.x_crystal),
                            average="weighted",
                        )
                    }
                )

            if self.test_md:
                print("test md...")
                predictions = model_obj.predict(self.x_md_test)
                print(len(self.x_md_test), len(predictions), len(self.y_md_test))

                acc_score = accuracy_score(self.y_md_test, predictions)
                f1score = f1_score(self.y_md_test, predictions, average="weighted")
                print("md test acc:           {:.2f}".format(acc_score))
                print("md test f1 score:      {:.2f}".format(f1score))
                wandb.log(
                    {
                        "md_test_acc": accuracy_score(
                            self.y_md_test, model_obj.predict(self.x_md_test)
                        )
                    }
                )
                wandb.log(
                    {
                        "md_test_f1": f1_score(
                            self.y_md_test,
                            model_obj.predict(self.x_md_test),
                            average="weighted",
                        )
                    }
                )

        run.finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="options for hyperparam tune")
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
        "--aug", action="store_true", dest="aug", default=False, help="augment data"
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
    if method == "bayes":
        sweep_config["metric"] = {"name": "md_test_f1", "goal": "maximize"}
    # trainer = training(model=model, pca_tf=pca_tf, test_crystal=True, test_md=True)
    # config = config(nestimators=100, max_depth=6)
    # trainer.train(config)

    sweep_id = wandb.sweep(sweep_config, project="HemeML_MD")
    training_obj = training(model=model, pca_tf=pca_tf, test_crystal=False, test_md=True)
    wandb.agent(sweep_id, function=training_obj.train, count=count)
