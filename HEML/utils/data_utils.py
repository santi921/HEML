import pandas as pd 
import numpy as np 
from glob import glob 
import os 
from sklearn.decomposition import PCA
from plotly.subplots import make_subplots
import plotly.graph_objects as go


# ./cpet/efield_cox_5o4k.dat
def mat_pull(file):

    with open(file) as f: 
        lines = f.readlines()

    steps_x = 2 * int(lines[0].split()[2]) + 1
    steps_y = 2 * int(lines[0].split()[3]) + 1
    steps_z = 2 * int(lines[0].split()[4][:-1]) + 1
    mat = np.zeros((steps_x, steps_y, steps_z, 3))

    #gap_x = round(np.abs(float(lines[steps_x*steps_y + 7].split()[0]) - float(lines[7].split()[0])), 4)
    #gap_y = round(np.abs(float(lines[steps_x+8].split()[1]) - float(lines[7].split()[1])), 4)
    #gap_z = round(np.abs(float(lines[8].split()[2]) - float(lines[7].split()[2])), 4)
    
    for ind, i in enumerate(lines[7:]):
        line_split = i.split()
        #print(i)
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 0] = float(line_split[-3])
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 1] = float(line_split[-2])
        mat[int(ind/(steps_z*steps_y)), int(ind/steps_z % steps_y), ind%steps_z, 2] = float(line_split[-1])
    return mat  

def pull_mats_w_label(dir_dat):

    x, y = [], []
    df = pd.read_csv("../../data/protein_data.csv")
    y_count, h_count, c_count = 0, 0, 0
    for row in df.iterrows():
        print(row[1]['name'])
        cpet_name = "../../data/cpet/efield_cox_" + row[1]['name'] + ".dat"
        
        if(os.path.exists(cpet_name)):
            x.append(mat_pull(cpet_name))
            if(row[1]['label'] == 'Y'):
                y.append([1,0,0])
                y_count += 1
            elif(row[1]['label'] == 'H'):
                y.append([0,1,0])
                h_count += 1
            else:
                y.append([0,0,1])
                c_count += 1
    print(y_count, h_count, c_count)
    return np.array(x), np.array(y)

def aug_all(mat, target, xy = True, z = False, mut = False):
    full_aug  = []
    full_aug_target = []

    for ind, i in enumerate(mat):
        x_aug, y_aug = augment_mat_field(i, target[ind])
        [full_aug.append(j) for j in x_aug]
        [full_aug_target.append(j) for j in y_aug]

    return np.array(full_aug), np.array(full_aug_target)

def augment_mat_field(mat, target, xy = True, z = False, mut = False):
    aug_target = []
    aug_mat = []

    if(xy):
        x_flip = np.array(np.flip(mat, axis = 0), dtype=float)
        y_flip = np.array(np.flip(mat, axis = 1), dtype=float)
        xy_flip = np.array(np.flip(np.flip(mat, axis = 1), axis = 0), dtype=float)

        x_flip[:,:,:,0] = -1*x_flip[:,:,:,0]
        y_flip[:,:,:,1] = -1*y_flip[:,:,:,1]
        xy_flip[:,:,:,0] = -1*xy_flip[:,:,:,0]
        xy_flip[:,:,:,1] = -1*xy_flip[:,:,:,1]
        
        aug_mat.append(mat)
        aug_mat.append(x_flip)
        aug_mat.append(y_flip)
        aug_mat.append(xy_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        
    if(z):
        z_flip = np.array(np.flip(mat, axis = 2), dtype=float)
        xz_flip = np.array(np.flip(np.flip(mat, axis = 0), axis = 2), dtype=float)
        yz_flip = np.array(np.flip(np.flip(mat, axis = 1), axis = 2), dtype=float)
        xyz_flip = np.array(np.flip(np.flip(np.flip(mat, axis = 2), axis = 1), axis = 0), dtype=float)

        z_flip[:,:,:,0] = -1*z_flip[:,:,:,0]
        xz_flip[:,:,:,0] = -1*xz_flip[:,:,:,0]
        xz_flip[:,:,:,2] = -1*xz_flip[:,:,:,2]
        yz_flip[:,:,:,1] = -1*yz_flip[:,:,:,1]
        yz_flip[:,:,:,2] = -1*yz_flip[:,:,:,2]
        xyz_flip[:,:,:,0] = -1*xyz_flip[:,:,:,0]
        xyz_flip[:,:,:,1] = -1*xyz_flip[:,:,:,1]
        xyz_flip[:,:,:,2] = -1*xyz_flip[:,:,:,2]
        
        aug_mat.append(z_flip)
        aug_mat.append(xz_flip)
        aug_mat.append(yz_flip)
        aug_mat.append(xyz_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        
    return aug_mat, aug_target
      
    
# try to sparify image data 
# alternatively use just the pca compression as input
def pca(mat): 
    mat_transform = mat.reshape(mat.shape[0], mat.shape[1] * mat.shape[2] * mat.shape[3] * mat.shape[4])
    pca = PCA(n_components=10)
    mat_transform = pca.fit_transform(mat_transform)
    cum_explained_var = []
    for i in range(0, len(pca.explained_variance_ratio_)):
        if i == 0:
            cum_explained_var.append(pca.explained_variance_ratio_[i])
        else:
            cum_explained_var.append(pca.explained_variance_ratio_[i] + 
                                    cum_explained_var[i-1])

    pc0 = pca.components_[0]
    print(np.shape(pc0))
    pc0 = pc0.reshape(1, mat.shape[1], mat.shape[2], mat.shape[3], mat.shape[4])

    fig = make_subplots(
        rows=1, cols=1,
        specs=[[{'type': 'cone'}]
            ])
    x, y, z = np.meshgrid(np.arange(-3, 2.8, 0.2),
                        np.arange(-3, 2.8, 0.2),
                        np.arange(-3, 2.8, 0.2))
    u = pc0[0][:,:,:,0].flatten()
    v = pc0[0][:,:,:,1].flatten()
    w = pc0[0][:,:,:,2].flatten()
    fig.add_trace(
        go.Cone(x=x.flatten(), y=y.flatten(), z=z.flatten(), u=u, v=v, w=w),
        row=1, col=1)

    fig.write_html("out_pca.html")

    print(cum_explained_var)
    
    return mat_transform, pca

def unwrap_pca(mat, pca, shape): 
    mat = pca.inverse_transform(mat)
    mat = mat.reshape(len(mat), shape[1], shape[2], shape[3], shape[4])
    return mat 

#def unwrap_and_show_comp(mat, pca, shape, comp = 1 ):

    