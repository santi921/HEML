from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from imblearn.ensemble import BalancedRandomForestClassifier, EasyEnsembleClassifier


def hyperparameter_dicts():
    dict_hyper = {}

    dict_ada = {
        "nestimators": {"min": 50, "max": 600},
        "lr": {"min": 1e-5, "max": 1e-1, "distribution": "log_uniform_values"},
    }
    dict_xgb = {
        "nestimators": {"min": 50, "max": 600},
        "max_depth": {"min": 2, "max": 8},
        "subsample": {"min": 0.4, "max": 1.0},
        "eta": {"min": 0.1, "max": 0.9},
        "gamma": {"min": 0.00001, "max": 100, "distribution": "log_uniform_values"},
        "reg_lambda": {
            "min": 0.000001,
            "max": 0.0001,
            "distribution": "log_uniform_values",
        },
        "alpha": {"min": 0.000001, "max": 0.001, "distribution": "log_uniform_values"},
    }
    dict_rfc = {
        "nestimators": {"min": 50, "max": 600},
        "max_depth": {"min": 2, "max": 8},
        "bootstrap": {"values": [False, True]},
    }
    dict_log = {
        "C": {"min": 1e-3, "max": 100, "distribution": "log_uniform_values"},
        "maxiter": {"min": 1000, "max": 10000},
        "l1": {"min": 0.000001, "max": 0.01, "distribution": "log_uniform_values"},
    }
    dict_qda = {
        "reg_param": {"min": 1e-5, "max": 0.7, "distribution": "log_uniform_values"},
    }
    dict_knc = {"nneighbors": {"min": 3, "max": 9}, "p": {"min": 1, "max": 5}}
    dict_ridge = {
        "alpha": {"min": 1e-3, "max": 100, "distribution": "log_uniform_values"},
    }
    dict_eec = {"nestimators": {"min": 50, "max": 600}}
    dict_brfc = {
        "nestimators": {"min": 50, "max": 600},
        "max_depth": {"min": 3, "max": 8},
        "bootstrap": {"values": [False, True]},
        "min_samples_leaf": {"min": 1, "max": 3},
    }

    dict_hyper["ada"] = dict_ada
    dict_hyper["xgb"] = dict_xgb
    dict_hyper["rfc"] = dict_rfc
    dict_hyper["log"] = dict_log
    dict_hyper["qda"] = dict_qda
    dict_hyper["knc"] = dict_knc
    dict_hyper["ridge"] = dict_ridge
    dict_hyper["eec"] = dict_eec
    dict_hyper["brfc"] = dict_brfc

    return dict_hyper


def construct_models(config, model="xgb"):

    if model == "ada":
        model_obj = AdaBoostClassifier(
            learning_rate=config.lr, n_estimators=config.nestimators
        )

    elif model == "xgb":
        model_obj = XGBClassifier(
            eta=config.eta,
            gamma=config.gamma,
            max_depth=config.max_depth,
            subsample=config.subsample,
            reg_lambda=config.reg_lambda,
            alpha=config.alpha,
        )

    elif model == "rfc":
        model_obj = RandomForestClassifier(
            n_estimators=config.nestimators,
            max_depth=config.max_depth,
            bootstrap=config.bootstrap,
            n_jobs=-1,
        )

    elif model == "log":
        model_obj = LogisticRegression(
            C=config.C, max_iter=config.maxiter, l1_ratio=config.l1
        )

    elif model == "qda":
        model_obj = QuadraticDiscriminantAnalysis(reg_param=config.reg_param,)

    elif model == "knc":
        model_obj = KNeighborsClassifier(
            n_neighbors=config.nneighbors, p=config.p, n_jobs=-1
        )

    elif model == "ridge":
        model_obj = RidgeClassifier(alpha=config.alpha)

    elif model == "eec":
        model_obj = EasyEnsembleClassifier(n_estimators=config.nestimators, n_jobs=-1)

    elif model == "brfc":
        model_obj = BalancedRandomForestClassifier(
            n_estimators=config.nestimators,
            min_samples_leaf=config.min_samples_leaf,
            max_depth=config.max_depth,
            bootstrap=config.bootstrap,
            n_jobs=-1,
        )

    else:
        model_obj = XGBClassifier()

    return model_obj
