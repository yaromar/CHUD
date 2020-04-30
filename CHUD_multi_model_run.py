#!/usr/bin/env python
# coding: utf-8

# ### Prepare Data

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import plot_precision_recall_curve
import matplotlib.pyplot as plt

from joblib import dump, load

print("Loading data...")

DATA_ROOT = "/ysm-gpfs/project/smz25/CHUD/"
DATA_FILE = "sample_var_phenos_leuk_topmed_v2.feather"

sample_var_phenos_file = DATA_ROOT + DATA_FILE
sample_var_phenos = pd.read_feather(sample_var_phenos_file)

print("Data loaded from feather")

# Drop the pneumonia columns and id

sample_var_phenos.drop(
    columns=['All_Pneumonia', 
            'AgeAt_All_Pneumonia', 
            'Incd_All_Pneumonia', 
            'All_Pneumonia_FollowUp', 
            'Prev_All_Pneumonia',
            'id',
            "AD",
            "CHIP_Gene"], 
    inplace = True)
sample_var_phenos.drop(
    sample_var_phenos.filter(
        regex = 'Prev\_|Incd\_|\_FollowUp', 
        axis="columns"), 
    axis = 1, inplace = True)


# Drop ADs and CHIP_Gene columns
sample_var_phenos.drop(
    sample_var_phenos.filter(
        regex = 'AD\_|CHIP\_Gene\_', 
        axis="columns"), 
    axis = 1, inplace = True)

out = open("perf_metric_all.csv", "w+")
out.write("modelname,accuracy,precision,recall\n")
out.flush()

for PREDICTION_PHENO in ["Coronary_Artery_Disease_SOFT", "hasCHIP", "AML", "MPN"]:
    print("Phenotype: ", PREDICTION_PHENO)
    # Define the Prediction Phenotype
    # PREDICTION_PHENO = 

    # 'Hypercholesterolemia',
    # 'Hypertension',
    # 'hasCHIP',
    # 'AML',
    # 'MPN',
    # 'SmokingStatus',
    # 'ever_smoked',
    # 'Coronary_Artery_Disease_SOFT',
    # 'Heart_Failure',

    # get dummy var columns for the two categorical variables

    # sample_phenos_onehot = pd.get_dummies(sample_var_phenos , columns=["CHIP_Gene", "AD"])

    # print("Data Dummified")
    # Do the imputation before we move forward. First drop rows with NA PCs
    sample_phenos_onehot = sample_var_phenos
    sample_phenos_onehot.dropna(subset=["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"], inplace=True)

    # Fill in an NA vals with 0 when we can

    fillnavals = {"Heart_Failure":0, "AF":0, "AML":0, "MPN":0}
    sample_phenos_onehot.fillna(value=fillnavals)

    print("Cleaned of NAs in data")

    # Create the data and training labels from the pandas dataframe

    all_except_class = sample_phenos_onehot.columns[sample_phenos_onehot.columns != PREDICTION_PHENO]

    X = sample_phenos_onehot.loc[:, all_except_class]

    # Scale data

    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler(copy=False)
    X = scaler.fit_transform(X)

    print("Data scaled")


    Y = np.array(sample_phenos_onehot.loc[:, PREDICTION_PHENO])
    np.nan_to_num(Y, copy=False)
    print("Created Numpy arrays")

    # Impute the rest

    from sklearn.impute import SimpleImputer
    imp = SimpleImputer(strategy='mean')
    X = imp.fit_transform(X)

    print("Imputed remaining NAs")

    # Check no more nans
    # nanargs = np.argwhere(np.isnan(X))

    # print(nanargs)
    # nanargs[:, 1]
    # print(np.unique(nanargs[:,1]).size)
    # sample_phenos_onehot.dtypes.iloc[list(np.unique(nanargs[:,1]))]
    # sample_phenos_onehot.iloc[:, list(np.unique(nanargs[:,1]))].isna().sum()

    # PCA

    # pca = PCA(n_components=3)
    # X = pca.fit_transform(X)

    # fig = plt.figure(1, figsize=(4, 3))
    # plt.clf()
    # ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)

    # plt.cla()
    # pca.fit(X)

    # Create the train test sets

    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.33, random_state=42, stratify=Y)
    print("Created Train/Test sets")

    print(X_train.shape)
    print(X_test.shape)
    print(y_train.shape)
    print(y_test.shape)

    # Set all NaN labels to 0
    np.nan_to_num(y_train, copy=False)
    np.nan_to_num(y_test, copy=False)


    # Oversample data 


    from imblearn.over_sampling import SMOTE, ADASYN, SMOTENC

    from collections import Counter 

    # Clean up a bit

    # del sample_phenos_onehot
#     del sample_var_phenos

    print("Oversampling Minority Class")


    counter = Counter(y_train)
    print("Before resampling: ", counter)

    sampler = ADASYN(sampling_strategy=1.0, n_jobs=-1)
    X_train_resampled, y_train_resampled = sampler.fit_resample(X_train, y_train)

    print("After resampling: ", sorted(Counter(y_train_resampled).items()))

    ### Model Training
    models = {}

    # ## Elastic Net Linear

    print("Starting Elastic Net Linear Regression...")

    from sklearn.linear_model import ElasticNet

    regr = ElasticNet(# cv=5, 
                        random_state=0,
                        # verbose=2,
                        # selection='random',
                        # n_jobs=-1
                       )
    regr.fit(X_train_resampled, y_train_resampled)

    print("Completed Training: Elastic Net Linear Regression")

    y_hat = regr.predict(X_test)

    dump(regr, 'sample_var_phenos_leuk_somatic_topmed_v2_elasticlinear.joblib') 

    yhat_rounded = np.where(y_hat > 0.5, 1, 0)

    importancedf = pd.DataFrame(columns = sample_phenos_onehot.loc[:, all_except_class].columns)
    importancedf.loc[len(importancedf)] = regr.coef_
    importancedf.sort_values(by=0, axis=1, ascending=False, inplace=True)
    importancedf.to_csv(DATA_FILE + "_" + PREDICTION_PHENO +"_linear_coefscaled.csv" , index=False, float_format='%.8f', mode="w", )

    linear_acc = sum(yhat_rounded == y_test)/len(yhat_rounded)
    linear_precision = precision_score(y_test, yhat_rounded)
    linear_recall = recall_score(y_test, yhat_rounded)
    print(f"Elastic Net Linear Regression Accuracy: {linear_acc}")
    print(f"Elastic Net Linear Regression Precision: {linear_precision}")
    print(f"Elastic Net Linear Regression Recall: {linear_recall}")
    out.write(f"{PREDICTION_PHENO}_elasticnetlinear,{linear_acc},{linear_precision},{linear_recall}\n")
    out.flush()

    models['elasticnetlinear'] = [regr, yhat_rounded]


    # Elastic Net LogReg

    print("Starting Elastic Net LogReg...")

    from sklearn.linear_model import SGDClassifier

    sgdregr = SGDClassifier(
        loss="log", 
        penalty="elasticnet",
        l1_ratio=0.5,
        n_jobs=-1,
        random_state=30,
        class_weight="balanced",
        # verbose=100,
        )

    sgdregr.fit(X_train_resampled, y_train_resampled)

    print("Completed Training: Elastic Net LogReg")

    sgd_y_hat = sgdregr.predict(X_test)

    dump(sgdregr, 'sample_var_phenos_leuk_somatic_topmed_v2_elasticlogreg.joblib') 


    importancedf = pd.DataFrame(columns = sample_phenos_onehot.loc[:, all_except_class].columns)
    importancedf.loc[len(importancedf)] = sgdregr.coef_[0]
    importancedf.sort_values(by=0, axis=1, ascending=False, inplace=True)
    importancedf.to_csv(DATA_FILE + "_" + PREDICTION_PHENO +"_logreg_coefscaled.csv" , index=False, float_format='%.8f', mode="w", )


    logreg_acc = sum(sgd_y_hat == y_test)/len(sgd_y_hat)
    logreg_precision = precision_score(y_test, sgd_y_hat)
    logreg_recall = recall_score(y_test, sgd_y_hat)
    print(f"Elastic Net LogReg Accuracy: {logreg_acc}")
    print(f"Elastic Net LogReg Precision: {logreg_precision}")
    print(f"Elastic Net LogReg Recall: {logreg_recall}")
    out.write(f"{PREDICTION_PHENO}_elasticnetlogreg,{logreg_acc},{logreg_precision},{logreg_recall}\n")
    out.flush()

    models['elasticnetlogreg'] = sgdregr


    # ## Random Forests

    print("Starting Random Forests...")


    from sklearn.ensemble import RandomForestClassifier
    rdmforest_clf = RandomForestClassifier(max_depth=2, 
                                 random_state=50,
                                 # verbose=100, 
                                 #class_weight="balanced_subsample",
                                 class_weight={0: 0.3, 1: 0.7},
                                 max_samples=1000,
                                 n_jobs=-1)


    rdmforest_clf.fit(X_train_resampled, y_train_resampled)

    print("Completed Training: Random Forests")

    forest_yhat = rdmforest_clf.predict(X_test)

    importancedf = pd.DataFrame(columns = sample_phenos_onehot.loc[:, all_except_class].columns)
    importancedf.loc[len(importancedf)] = rdmforest_clf.feature_importances_
    importancedf.sort_values(by=0, axis=1, ascending=False, inplace=True)
    importancedf.to_csv(DATA_FILE + "_" + PREDICTION_PHENO +"_randforest_importancescaled.csv" , index=False, float_format='%.8f', mode="w", )

    forest_acc = sum(forest_yhat == y_test)/len(forest_yhat)
    forest_precision = precision_score(y_test, forest_yhat)
    forest_recall = recall_score(y_test, forest_yhat)

    print(f"Random Forest Accuracy: {forest_acc}")
    print(f"Random Forest Precision: {forest_precision}")
    print(f"Random Forest Recall: {forest_recall}")
    out.write(f"{PREDICTION_PHENO}_randomforest,{forest_acc},{forest_precision},{forest_recall}\n")
    out.flush()

    dump(rdmforest_clf, 'sample_var_phenos_leuk_somatic_topmed_v2_elasticrandforest_resampled.joblib') 

    models['randomforest'] = rdmforest_clf

    # ## Support Vectors

    from sklearn.svm import SVC
    from sklearn import datasets, svm
    from sklearn.kernel_approximation import Nystroem


    # #### Approximate RBF Kernel Map and Linear SVM

    print("Starting RBF Kernel SVM...")

    linearsvm_clf = svm.LinearSVC()
    feature_map_nystroem = Nystroem(gamma=.2,
                                    random_state=100,
                                    n_components=300)

    data_transformed = feature_map_nystroem.fit_transform(X_train_resampled)

    linearsvm_clf.fit(data_transformed, y_train_resampled)

    print("Trained: RBF Kernel SVM")

    # evaluation on test set

    X_test_transformed = feature_map_nystroem.fit_transform(X_test)
    linearsvm_yhat = linearsvm_clf.predict(X_test_transformed)

    rbf_acc = sum(linearsvm_yhat == y_test)/len(linearsvm_yhat)
    rbf_precision = precision_score(y_test, linearsvm_yhat)
    rbf_recall = recall_score(y_test, linearsvm_yhat)

    print(f"Nystroem and Linear SVM Accuracy: {rbf_acc}")
    print(f"Nystroem and Linear SVM Precision: {rbf_precision}")
    print(f"Nystroem and Linear SVM Recall: {rbf_recall}")
    out.write(f"{PREDICTION_PHENO}_rbfsvmapprox,{rbf_acc},{rbf_precision},{rbf_recall}\n")
    out.flush()

    models['rbfsvmapprox'] = linearsvm_clf

    # Linear SVM On its own

    print("Staring Linear Kernel SVM...")

    linearsvmalone_clf = svm.LinearSVC()
    linearsvmalone_clf.fit(X_train_resampled, y_train_resampled)

    print("Trained: Linear Kernel SVM")

    linearsvmalone_yhat = linearsvmalone_clf.predict(X_test)


    importancedf = pd.DataFrame(columns = sample_phenos_onehot.loc[:, all_except_class].columns)
    importancedf.loc[len(importancedf)] = linearsvmalone_clf.coef_[0]
    importancedf.sort_values(by=0, axis=1, ascending=False, inplace=True)
    importancedf.to_csv(DATA_FILE + "_" + PREDICTION_PHENO +"_linearsvm_coefscaled.csv" , index=False, float_format='%.8f', mode="w", )

    linearsvm_acc = sum(linearsvmalone_yhat == y_test)/len(linearsvmalone_yhat)
    linearsvm_precision = precision_score(y_test, linearsvmalone_yhat)
    linearsvm_recall = recall_score(y_test, linearsvmalone_yhat)

    print(f"Linear SVM Accuracy: {linearsvm_acc}")
    print(f"Linear SVM Precision: {linearsvm_precision}")
    print(f"Linear SVM Recall: {linearsvm_recall}")
    out.write(f"{PREDICTION_PHENO}_linearsvm,{linearsvm_acc},{linearsvm_precision},{linearsvm_recall}\n")
    out.flush()

    models['linearsvm'] = linearsvmalone_clf

    # Make up some plots
    from sklearn.metrics import plot_roc_curve
    from sklearn.metrics import roc_curve, auc, RocCurveDisplay
    from sklearn.metrics import plot_confusion_matrix
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import ConfusionMatrixDisplay

    for model_name, model in models.items():
        print("Model: ", model_name)
        # ROC Curve
        if model_name == "elasticnetlinear":

            fpr, tpr, thresholds = roc_curve(model[1], y_test)
            roc_auc = auc(fpr, tpr)
            svc_disp = RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc,
                                            estimator_name='Elastic Net Linear')

            svc_disp.plot()
            plt.title(f"ROC Curve of {PREDICTION_PHENO} by {model_name} (Scaled)")
            plt.savefig(f"finalplots/scaledroc_{PREDICTION_PHENO}_{model_name}.png", 
                        dpi=600, transparent=True,
                       bbox_inches="tight",
                       pad_inches=0.3)

        else:
            if model_name == "rbfsvmapprox":
                svc_disp = plot_roc_curve(model, feature_map_nystroem.fit_transform(X_test), y_test)
            else:
                svc_disp = plot_roc_curve(model, X_test, y_test)
            svc_disp.plot()
            plt.title(f"ROC Curve of {PREDICTION_PHENO} by {model_name} (Scaled)")
            plt.savefig(f"finalplots/scaledroc_{PREDICTION_PHENO}_{model_name}.png", 
                        dpi=600, transparent=True,
                       bbox_inches="tight",
                       pad_inches=0.3)


        # Confusion Matrix
        if model_name == "elasticnetlinear":
            conf_mat = confusion_matrix(model[1], y_test)
            disp = ConfusionMatrixDisplay(confusion_matrix=conf_mat, display_labels=[0,1])
            disp.plot()
            plt.title(f"Confusion Matrix of {PREDICTION_PHENO} by {model_name} (Scaled)")        
            plt.savefig(f"finalplots/scaledconfusion_mat_{PREDICTION_PHENO}_{model_name}.png", 
                        dpi=600, transparent=True,
                       bbox_inches="tight",
                       pad_inches=0.3)

        else:
            if model_name == "rbfsvmapprox":
                disp = plot_confusion_matrix(model, feature_map_nystroem.fit_transform(X_test), y_test,
                                             display_labels=[0,1],
                                             cmap=plt.cm.Blues,
                                             normalize=None)
            else:
                disp = plot_confusion_matrix(model, X_test, y_test,
                                             display_labels=[0,1],
                                             cmap=plt.cm.Blues,
                                             normalize=None)
            disp.plot()
            plt.title(f"Confusion Matrix of {PREDICTION_PHENO} by {model_name} (Scaled)")        
            plt.savefig(f"finalplots/scaledconfusion_mat_{PREDICTION_PHENO}_{model_name}.png", 
                        dpi=600, transparent=True,
                       bbox_inches="tight",
                       pad_inches=0.3)

        # Didn't plot PR curve, since they show nothing at all

out.close()
