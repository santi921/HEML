import numpy as np

from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import accuracy_score, f1_score

from HEML.utils.data import pull_mats_from_MD_folder, pull_mats_w_label
from HEML.utils.attrib import *
from HEML.utils.model import *
from HEML.utils.fields import pca, aug_all


class training:
    def __init__(self, model, pca_tf=True, aug=True, test_crystal=False, test_md=False,  joint_pca=False):
        self.aug = aug
        self.pca_tf = pca_tf
        self.model = model
        self.test_crystal = test_crystal
        self.test_md = test_md
        self.joint_pca = joint_pca

        pca_comps = 25

        # df = pd.read_csv("../../data/protein_data.csv")
        x, y, names = pull_mats_from_MD_folder(label_ind=3)

        x_sign = np.sign(x)
        # getting absolute value of every element
        x_abs = np.abs(x)
        # applying log1p
        x_log1p = np.log1p(x_abs)
        # getting sign back
        x = np.multiply(x_log1p, x_sign)
        y = [np.argmax(i) for i in y]
        print(set(y))
        (
            self.X_train,
            self.X_test,
            self.y_train,
            self.y_test,
            self.names_train,
            self.names_test,
        ) = train_test_split(x, y, names, test_size=0.1, random_state=0)
        # print(self.y_test.shape)

        if self.test_md:
            x_md_test, y_md_test, names_test = pull_mats_from_MD_folder(
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

            if self.aug:
                x_md_test, y_md_test = aug_all(x_md_test, y_md_test, xy=True, z=False)
                # create new list with repeated element in names_test
            
            names_test = [item for item in names_test for i in range(4)]
            self.y_md_test = [np.argmax(i) for i in y_md_test]
            self.names_md_test = names_test

        if self.pca_tf:
            self.X_train_untransformed = self.X_train
            self.X_test_untransformed = self.X_test

            if self.aug:
                x_train_temp, y_train_temp = aug_all(
                    np.concatenate((self.X_train, self.X_test)),
                    np.concatenate((self.y_train, self.y_test)),
                    xy=True,
                    z=True,
                )
                if self.joint_pca:
                    _, self.pca_obj = pca(np.concatenate((self.X_train, self.X_test, x_md_test)), verbose=True, pca_comps=pca_comps)

                else:
                    _, self.pca_obj = pca(x_train_temp, verbose=True, pca_comps=pca_comps)
           
            else:    
                _, self.pca_obj = pca(self.X_train, verbose=True, pca_comps=pca_comps)


            self.X_train, self.pca_obj_train = pca(self.X_train, self.pca_obj)
            self.X_test, self.pca_obj_test = pca(self.X_test, self.pca_obj)

        if self.test_md:
            self.x_md_test, _ = pca(x_md_test, self.pca_obj)


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

    def train(self, config):
        model_obj = self.make_model(config)

        kf = KFold(n_splits=5, random_state=0, shuffle=True)
        acc_train, acc_val, f1_val = [], [], []

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
                x_train_temp, y_train_temp = aug_all(
                    self.X_train_untransformed[ind_train], y_train_temp, xy=True, z=True
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
            f1_val.append(f1_score(y_val_temp, y_val_pred, average="weighted"))

        print("final testing...")
        y_test_pred = model_obj.predict(self.X_test)
        print(set(y_test_pred))
        acc_test = accuracy_score(self.y_test, y_test_pred)
        f1_test = f1_score(self.y_test, y_test_pred, average="weighted")

        print("train acc:        {:.2f}".format(np.mean(np.array(acc_train))))
        print("validation acc:   {:.2f}".format(np.mean(np.array(acc_val))))
        print("test acc:         {:.2f}".format(acc_test))
        print("validation f1:    {:.2f}".format(np.mean(np.array(f1_val))))
        print("test f1:          {:.2f}".format(f1_test))

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
                    accuracy_score(self.y_crystal, model_obj.predict(self.x_crystal))
                )
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

        if self.test_md:
            print("test md...")
            predictions = model_obj.predict(self.x_md_test)
            # predictions = []
            # for i in range(len(self.x_md_test)):
            #    prediction_single = model_obj.predict(self.x_md_test[i].reshape(1, -1))
            #    predictions += prediction_single.tolist()

            # predictions = model_obj.predict(self.x_md_test)
            # print(len(self.x_md_test), len(predictions), len(self.y_md_test))

            acc_score = accuracy_score(self.y_md_test, predictions)
            f1score = f1_score(self.y_md_test, predictions, average="weighted")
            print("md test acc:           {:.2f}".format(acc_score))
            print("md test f1 score:      {:.2f}".format(f1score))

            names = self.names_md_test
            result = {}
            print(set(predictions))
            for i in range(len(predictions)):
                name_pro = names[i]

                if name_pro not in result:
                    result[name_pro] = [0, 0, 0]

                pred = predictions[i]
                if pred == 0:
                    result[name_pro][0] += 1
                elif pred == 1:
                    result[name_pro][1] += 1
                else:
                    result[name_pro][2] += 1
            print(result)


class config:
    def __init__(
        self,
        nestimators,
        max_depth,
        eta,
        gamma,
        reg_lambda,
        alpha,
        subsample,
        min_samples_leaf,
        bootstrap,
    ):
        self.nestimators = nestimators
        self.max_depth = max_depth
        self.eta = eta
        self.gamma = gamma
        self.reg_lambda = reg_lambda
        self.alpha = alpha
        self.subsample = subsample
        self.min_samples_leaf = min_samples_leaf
        self.bootstrap = bootstrap


if __name__ == "__main__":
    pca_tf = True
    model = "xgb"  # "xgb" "brfc" "eec"

    trainer = training(
        model=model, 
        pca_tf=pca_tf, 
        test_crystal=False, 
        test_md=True, 
        aug=True, 
        joint_pca=True

    )

    config_best_joint_pca = config(
        nestimators=800,
        #nestimators=400, 
        #max_depth=5,
        max_depth=4,
        eta=0.70399291160943,
        gamma=63.27723608448208,
        #gamma=50.0,
        #reg_lambda=0.0,
        reg_lambda=0.0,
        alpha=0.0009978623441125043,
        subsample=0.8322187403324615,
        min_samples_leaf=3,
        bootstrap=True,
    )

    """config = config(
        nestimators=800,
        #nestimators=400, 
        min_samples_leaf=5,
        bootstrap=False,
        eta=0.7,
        gamma=63,
        max_depth=4,
        subsample=0.63,
        reg_lambda=0.000002194560486561219,
        alpha=0.0009978623441125043,
        #eval_metric="mlogloss",
    )"""
    config = config_best_joint_pca

    trainer.train(config)
