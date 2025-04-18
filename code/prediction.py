import numpy as np
import math
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
import scipy as sp
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import KFold
import random


def gradient_boosting(X_train,Y_train,X_test,Y_test):
    params={'n_estimators': 10000, 'max_depth': 7, 'min_samples_split': 3,
                'learning_rate': 0.001,'max_features':'sqrt','subsample':0.5,
                }
    regr = GradientBoostingRegressor(**params)
    regr.fit(X_train,Y_train)
    pearson_coorelation = sp.stats.pearsonr(Y_test,regr.predict(X_test))
    mse = mean_squared_error(Y_test, regr.predict(X_test))
    mae = mean_absolute_error(Y_test, regr.predict(X_test))
    return pearson_coorelation[0],pow(mse,0.5),mae
    
    





def cv10():
    
    # protein-rna
    train_X = np.load('feature/pri-train-feature.npy')
    train_Y = np.load('feature/pri-train-label.npy')
    test_X = np.load('feature/pri-test-feature.npy')
    test_Y = np.load('feature/pri-test-label.npy')
    
    feature = np.vstack([train_X,test_X])
    label = np.array( np.array(train_Y.tolist()+test_Y.tolist()) )
    
    
    print(train_X.shape)
    
    
    #print(feature.shape,label.shape,np.max(feature),np.min(feature))
    cv = 10
    seed = 42
    #skf = StratifiedKFold(n_splits=cv, shuffle=True, random_state=seed)
    kf = KFold(n_splits=cv, shuffle=True, random_state=seed)
    
    
    pcc = 0
    rmse = 0
    mae = 0
    print('start')
    for now, (train_index, test_index) in enumerate(kf.split(feature)):
        train_feature, test_feature = feature[train_index], feature[test_index]
        train_label, test_label     = label[train_index], label[test_index]
        
        print('CV:',now)
        
    
        scaler = MinMaxScaler(feature_range=(-1, 1))
        train_feature = scaler.fit_transform(train_feature)
        test_feature = scaler.transform(test_feature)
        
        
        print('range:',np.max(train_feature),np.min(train_feature),np.max(test_feature),np.min(test_feature))
        
        pcc1,rmse1,mae1 = gradient_boosting(train_feature,train_label,test_feature,test_label)
        
        
        
        print('PCC:',np.round(pcc1,3),'RMSE:',np.round(rmse1,3),'MAE:',np.round(mae1,3))
        
        pcc = pcc + pcc1
        rmse = rmse + rmse1
        mae = mae + mae1
    
    
    print('PCC:',np.round(pcc/cv,3),'RMSE:',np.round(rmse/cv,3),'MAE:',np.round(mae/cv,3))


#cv10()