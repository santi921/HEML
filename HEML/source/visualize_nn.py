import tensorflow as tf

import numpy as np
import pandas as pd

from HEML.utils.data_utils import *
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, r2_score
from tensorflow.keras.regularizers import l1 as L1
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt

from HEML.utils.attrib_utils import *



def split_and_filter(mat, cutoff = 95, min_max = True, std_mean = False):

    arr_mean, arr_std, arr_min, arr_max  = np.mean(mat), np.std(mat), np.min(mat), np.max(mat)
    if(min_max):
        mat = (mat - arr_min) / (arr_max - arr_min + 10e-10)

    if(std_mean):
        mat = (mat - arr_mean) / (arr_std)
    
    try:
        u = mat[0][:,:,:,0].flatten()
        v = mat[0][:,:,:,1].flatten()
        w = mat[0][:,:,:,2].flatten()
    except:
        u = mat[:,:,:,0].flatten()
        v = mat[:,:,:,1].flatten()
        w = mat[:,:,:,2].flatten()

    component_distro = [np.sqrt(u[ind]**2 + v[ind]**2 + w[ind]**2) for ind in range(len(u))]
    cutoff = np.percentile(component_distro, cutoff)

    for ind, i in enumerate(component_distro): 
        if (i < cutoff): 
            u[ind], v[ind], w[ind] = 0,0,0  
    u = np.around(u, decimals=1)
    v = np.around(v, decimals=1)
    w = np.around(w, decimals=1)
    return u, v, w


if __name__ == "__main__":


    dense = True
    gradcam = False
    pca_tf = True

    no_classes = 3
    df = pd.read_csv("../../data/protein_data.csv")
    x, y = pull_mats_w_label('./dat')
    arr_min, arr_max,  = np.min(x), np.max(x)
    x = (x - arr_min) / (arr_max - arr_min + 1e-18)
    shape_mat = x.shape
    X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.25, random_state=1)
    #X_train, y_train = aug_all(X_train, y_train)
    #X_test, y_test = aug_all(X_test, y_test)
    
    if (pca_tf):
        print("::::pca::::" * 5)

        mat, pca_obj = pca(np.concatenate((X_train, X_test)))
        X_train = mat[:len(X_train)]
        X_test = mat[len(X_train):]
        img_1 = X_train[0,:,] # y_test label 3 is 0 
        img_2 = X_train[1,:] # y_test label 2 is 1 
        img_3 = X_train[27,:] # y_test label 1 is 27
        print(np.shape(img_1))
        mat_pca_inverse = unwrap_pca(mat, pca_obj, shape_mat)
        print(np.shape(mat_pca_inverse))


    if(dense):
        
        from xgboost import XGBClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.ensemble import RandomForestClassifier

        model_xgb = XGBClassifier(max_depth = 4, colsample_bytree = 0.2, subsample = 0.3, eta = 0.03)
        #model_xgb = RandomForestClassifier(max_depth=2, bootstrap=True, ccp_alpha=0.2, )
        flat_y_train = [np.argmax(i) for i in y_train]
        flat_y_test = [np.argmax(i) for i in y_test]

        if(pca_tf):
            model_xgb.fit(X_train, flat_y_train) 
            y_train_pred = model_xgb.predict(X_train)
            y_test_pred = model_xgb.predict(X_test)

            #model = tf.keras.Sequential([
            #tf.keras.layers.Flatten(input_shape=(10)),
            #tf.keras.layers.Dense(128, activation='relu'),
            #tf.keras.layers.Dense(64, activation='reli'),
            #tf.keras.layers.Dense(32, activation='reli'),
            #tf.keras.layers.Dense(np.shape(y)[1], activation = 'softmax', name = 'visualized_layer')
            #])
            #model.summary()

        else:

            model_xgb.fit(X_train.reshape(X_train.shape[0], X_train.shape[1] * X_train.shape[2] * X_train.shape[3] * X_train.shape[4]
            ), flat_y_train, verbose = True) 
            y_train_pred = model_xgb.predict(X_train.reshape(X_train.shape[0], X_train.shape[1] * X_train.shape[2] * X_train.shape[3] * X_train.shape[4]))
            y_test_pred = model_xgb.predict(X_test.reshape(X_test.shape[0], X_test.shape[1] * X_test.shape[2] * X_test.shape[3] * X_test.shape[4]))
            #model = tf.keras.Sequential([
            #tf.keras.layers.Flatten(input_shape=(30, 30, 30, 3)),
            #tf.keras.layers.Dense(128, activation='sigmoid'),
            #tf.keras.layers.Dense(64, activation='sigmoid'),
            #tf.keras.layers.Dense(32, activation='sigmoid'),
            #tf.keras.layers.Dense(np.shape(y)[1], activation = 'softmax', name = 'visualized_layer')
            #])

        print("test r^2: " + str(accuracy_score(flat_y_train, y_train_pred)))
        print("test r^2: " + str(accuracy_score(flat_y_test, y_test_pred)))
        
       
    else:
        model = tf.keras.Sequential([
            tf.keras.layers.Conv3D(128, (4, 4, 4), strides = (1,1,1), 
            activation='sigmoid',
            padding='same',
            kernel_regularizer=L1(10e-4),activity_regularizer = L1(10e-4),
            input_shape=(30, 30, 30, 3)),
            #tf.keras.layers.Conv3D(32, (3,3,3), strides = 1, activation="relu"),
            tf.keras.layers.AveragePooling3D(pool_size=4),
            tf.keras.layers.BatchNormalization(),
            tf.keras.layers.Conv3D(64, (3,3,3), activation="sigmoid"),
            tf.keras.layers.MaxPooling3D(pool_size=4),
            tf.keras.layers.BatchNormalization(),
            tf.keras.layers.GlobalAveragePooling3D(), #tf.keras.layers.Flatten(),
            tf.keras.layers.Flatten(),
            tf.keras.layers.Dense(128, activation='sigmoid'),
            tf.keras.layers.Dense(np.shape(y)[1], activation = 'softmax', name = 'visualized_layer')
        ])

    model.compile(optimizer='adam',
                loss=tf.keras.losses.CategoricalCrossentropy(),
                #loss = 'kullback_leibler_divergence', 
                metrics = ['accuracy']

                )


    if(pca_tf):
        print("pass")

    else:
        model.fit(X_train, y_train,epochs=10, validation_data = (X_test,y_test), verbose = True)
        y_pred = model.predict(X_test)

        img_1 = X_train[0,:,:,:,:].reshape([1,30,30,30,3])
        img_2 = X_train[1,:,:,:,:].reshape([1,30,30,30,3])
        img_3 = X_train[27,:,:,:,:].reshape([1,30,30,30,3])

        print("checking equality of original matrices")
        print(np.array_equal(img_1, img_2))
        print(np.array_equal(img_2, img_3))
        print(np.array_equal(img_1, img_3))

        if(dense):

            print("\n")
            print("::::saliency map::::" * 5)
            print("\n")

            mat_1 = saliency_map_dense(model, img_1)
            mat_2 = saliency_map_dense(model, img_2)
            mat_3 = saliency_map_dense(model, img_3)        
            
        else: 

            if (gradcam): 
                print("\n")
                print("::::gradcam::::" * 5)
                print("\n")

                mat_1 = make_gradcam_heatmap(img_1, model, model.layers[0].name)
                mat_2 = make_gradcam_heatmap(img_2, model, model.layers[0].name)
                mat_3 = make_gradcam_heatmap(img_3, model, model.layers[0].name)

            else:
                print("\n")
                print("::::integrated gradients::::" * 5)
                print("\n")

                orig_img = np.copy(img_1[0])
                preds = model.predict(img_1)    
                top_pred_idx = tf.argmax(preds[0])
                grads = get_gradients(img_1, model, top_pred_idx=top_pred_idx)
                igrads = random_baseline_integrated_gradients(
                    np.copy(orig_img), model = model, top_pred_idx=top_pred_idx, num_steps=50, num_runs=10 )
                mat_1 = igrads.numpy()
                
                orig_img = np.copy(img_2[0])
                preds = model.predict(img_2)    
                top_pred_idx = tf.argmax(preds[0])
                grads = get_gradients(img_2, model, top_pred_idx=top_pred_idx)
                igrads = random_baseline_integrated_gradients(
                    np.copy(orig_img), model = model, top_pred_idx=top_pred_idx, num_steps=50, num_runs=10 )
                mat_2 = igrads.numpy()
                
                orig_img = np.copy(img_3[0])
                preds = model.predict(img_3)    
                top_pred_idx = tf.argmax(preds[0])
                grads = get_gradients(img_3, model, top_pred_idx=top_pred_idx)
                igrads = random_baseline_integrated_gradients(
                    np.copy(orig_img), model = model, top_pred_idx=top_pred_idx, num_steps=50, num_runs=10 )
                mat_3 = igrads.numpy()

    
    y_pred_max = [np.argmax(i) for i in y_pred]
    y_test_arg = [np.argmax(i) for i in y_test]

    # todo: redo this for different sample types
    x, y, z = np.meshgrid(np.arange(-3, 2.8, 0.2),
                        np.arange(-3, 2.8, 0.2),
                        np.arange(-3, 2.8, 0.2))
    iron_size = 1.26
    resolution = 100
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    
    x_sphere = iron_size * np.cos(u)*np.sin(v) 
    y_sphere = iron_size * np.sin(u)*np.sin(v) 
    z_sphere = iron_size * np.cos(v) 
    #print(x_sphere)
    
    
    fig = make_subplots(
        rows=1, cols=3,
        subplot_titles=("cat 1", "cat 2", "cat 3"),
        specs=[[{'type': 'cone'}, {'type': 'cone'}, {'type': 'surface'}],
            ])

    u_1, v_1, w_1 = split_and_filter(mat_1)
    print(np.max(u_1), np.max(w_1), np.max(v_1))
    fig.add_trace(
        go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), u=u_1, v=v_1, w=w_1),
        row=1, col=1)
    #fig.add_trace(go.Surface(x=x_sphere, y=y_sphere, z=z_sphere, surfacecolor=x**2 + y**2 + z**2), row=1, col=1)

    
    u_2, v_2, w_2 = split_and_filter(mat_2)
    print(np.max(u_2), np.max(w_2), np.max(v_2))
    
    fig.add_trace(
        go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), u=u_2, v=v_2, w=w_2),
        row=1, col=2)
    #fig.add_trace(go.Surface(x=x_sphere, y=y_sphere, z=z_sphere, surfacecolor=x**2 + y**2 + z**2), row=1, col=2)

    u_3, v_3, w_3 = split_and_filter(mat_3)
    print(np.max(u_3), np.max(w_3), np.max(v_3))
        
    fig.add_trace(
        go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), u=u_3, v=v_3, w=w_3),
        row=1, col=3)
    #fig.add_trace(go.Surface(x=x_sphere, y=y_sphere, z=z_sphere), row=1, col=3)


    fig.show()
    fig.write_html("out_final.html")
    print("------Checking Equality Conditions------")
    print("...first and second example...")
    print(np.array_equal(u_1, u_2))
    print(np.array_equal(w_1, w_2))
    print(np.array_equal(v_1, v_2))
    print("...first and third example...")
    print(np.array_equal(u_1, u_3))
    print(np.array_equal(w_1, w_3))
    print(np.array_equal(v_1, v_3))
    print("...second and third example...")
    print(np.array_equal(u_2, u_3))
    print(np.array_equal(w_2, w_3))
    print(np.array_equal(v_2, v_3))    