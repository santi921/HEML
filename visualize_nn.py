import numpy as np
import pandas as pd

from data_utils import *
import tensorflow as tf
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score
from tensorflow.keras.regularizers import L1
#tf.compat.v1.disable_eager_execution()

import plotly.graph_objects as go


import matplotlib.pyplot as plt
import numpy as np


from attrib_utils import *

dense = False
no_classes = 3
df = pd.read_csv("protein_data.csv")
x, y = pull_mats_w_label('./dat')
arr_min, arr_max,  = np.min(x), np.max(x)
x = (x - arr_min) / (arr_max - arr_min + 1e-18)
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42)


model2 = tf.keras.Sequential([
    tf.keras.layers.Flatten(input_shape=(30, 30, 30, 3)),
    tf.keras.layers.Dense(128, activation='relu'),
    tf.keras.layers.Dense(np.shape(y)[1], activation = 'softmax', name = 'visualized_layer')
])

model = tf.keras.Sequential([
    tf.keras.layers.Conv3D(32, (4, 4, 4), strides = (1,1,1), 
    activation="relu",
    padding='same',
    kernel_regularizer=L1(10e-4),activity_regularizer =L1(10e-4),
    input_shape=(30, 30, 30, 3)),
    tf.keras.layers.Conv3D(32, (3,3,3), strides = 1, activation="relu"),
    tf.keras.layers.AveragePooling3D(pool_size=3),
    tf.keras.layers.BatchNormalization(),
    tf.keras.layers.Conv3D(64, (2,2,2), activation="relu"),
    tf.keras.layers.MaxPooling3D(pool_size=2),
    tf.keras.layers.BatchNormalization(),
    tf.keras.layers.GlobalAveragePooling3D(), #tf.keras.layers.Flatten(),
    tf.keras.layers.Flatten(),
    tf.keras.layers.Dense(128, activation="relu"),
    tf.keras.layers.Dense(np.shape(y)[1], activation = 'softmax', name = 'visualized_layer')
])

model2.compile(optimizer='adam',
              loss=tf.keras.losses.CategoricalCrossentropy(),
              #loss = 'kullback_leibler_divergence', 
              metrics = ['accuracy']

              )
model.compile(optimizer='adam',
              loss=tf.keras.losses.CategoricalCrossentropy(),
              #loss = 'kullback_leibler_divergence', 
              metrics = ['accuracy']

              )

model.fit(X_train, y_train,epochs=5, validation_data = (X_test,y_test))

if(dense):
    model2.fit(X_train, y_train,epochs=5, validation_data = (X_test,y_test))
    # works
    mat = saliency_map_dense(model2, X_train[0,:,:,:,:].reshape([1,30,30,30,3]))
    y_pred = model.predict(X_test)

else: 
    y_pred = model2.predict(X_test)
    # works for convnet
    img = X_train[0,:,:,:,:].reshape([1,30,30,30,3])
    orig_img = np.copy(img[0])
    preds = model.predict(img)
    top_pred_idx = tf.argmax(preds[0])
    grads = get_gradients(img, model, top_pred_idx=top_pred_idx)
    igrads = random_baseline_integrated_gradients(
        np.copy(orig_img), model = model, top_pred_idx=top_pred_idx, num_steps=50, num_runs=1
    )
    mat = igrads.numpy()
    print(mat)
    # grad cam
    model.summary()
    make_gradcam_heatmap(img, model, model.layers[0].name)



y_pred_max = [np.argmax(i) for i in y_pred]
y_test_arg = [np.argmax(i) for i in y_test]

x, y, z = np.meshgrid(np.arange(-3, 2.8, 0.2),
                      np.arange(-3, 2.8, 0.2),
                      np.arange(-3, 2.8, 0.2))


arr_mean, arr_std, arr_min, arr_max  = np.mean(mat), np.std(mat), np.min(mat), np.max(mat)
mat = (mat - arr_min) / (arr_max - arr_min)

print("mat shape: " + str(np.shape(mat)))
try:
    u = mat[0][:,:,:,0].flatten()
    v = mat[0][:,:,:,1].flatten()
    w = mat[0][:,:,:,2].flatten()
except:
    u = mat[:,:,:,0].flatten()
    v = mat[:,:,:,1].flatten()
    w = mat[:,:,:,2].flatten()


component_distro = [np.sqrt(u[ind]**2 + v[ind]**2 + w[ind]**2) for ind in range(len(u))]
cutoff = np.percentile(component_distro, 98)

print(cutoff)

for ind, i in enumerate(component_distro): 
    if (i < cutoff): 
        u[ind], v[ind], w[ind] = 0,0,0

fig = go.Figure(data=go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), u=u, v=v, w=w))
fig.show()