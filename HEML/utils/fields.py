import numpy as np 
from sklearn.cluster import AffinityPropagation
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from sklearn.decomposition import PCA

from HEML.utils.data import mat_pull

def split_and_filter(mat, cutoff=95, min_max=True, std_mean=False, log1=False):

    arr_mean, arr_std, arr_min, arr_max = (
        np.mean(mat),
        np.std(mat),
        np.min(mat),
        np.max(mat),
    )

    if log1:
        shape = mat.shape
         # get magnitude of vector field
        mags = np.sqrt((mat**2).sum(axis=3))
        #mags = mags.reshape([1, shape[0], shape[1], shape[2]])
        mags_log1p = np.log1p(mags)
        # get ratio magnitude to mags_log1p
        ratio = mags_log1p / mags
        multiply = np.multiply(mat, ratio.reshape([shape[0], shape[1], shape[2], 1])) 
        mat = multiply

    if min_max:
        mat = (mat - arr_min) / (arr_max - arr_min + 10e-10)

    if std_mean:
        mat = (mat - arr_mean) / (arr_std)

    try:
        u = mat[0][:, :, :, 0].flatten()
        v = mat[0][:, :, :, 1].flatten()
        w = mat[0][:, :, :, 2].flatten()
    except:
        u = mat[:, :, :, 0].flatten()
        v = mat[:, :, :, 1].flatten()
        w = mat[:, :, :, 2].flatten()

    component_distro = [
        np.sqrt(u[ind] ** 2 + v[ind] ** 2 + w[ind] ** 2) for ind in range(len(u))
    ]
    cutoff = np.percentile(component_distro, cutoff)

    for ind, i in enumerate(component_distro):
        if i < cutoff:
            u[ind], v[ind], w[ind] = 0, 0, 0

    u = np.around(u, decimals=2)
    v = np.around(v, decimals=2)
    w = np.around(w, decimals=2)

    return u, v, w


def save_numpy_as_dat(dict_meta_data, average_field, name):
    """
    Saves np array in original format from cpet output
    Takes: 
        dict_meta_data: dictionary with meta data
        average_field: np array with average field
        name: name of file to save
    """
    

    first_line = dict_meta_data["first_line"]
    steps_x = dict_meta_data["steps_x"]
    steps_y = dict_meta_data["steps_y"]
    steps_z = dict_meta_data["steps_z"]
    step_size_x = dict_meta_data["step_size_x"]
    step_size_y = dict_meta_data["step_size_y"]
    step_size_z = dict_meta_data["step_size_z"]
    print(dict_meta_data)
    lines = [first_line]
    # add six lines starting with #
    lines = lines + ["#\n"] * 6

    # write as 
    for i in range(steps_x):
        for j in range(steps_y):
            for k in range(steps_z):
                line = "{:.3f} {:.3f} {:.3f} {:.6f} {:.6f} {:.6f}\n".format(
                    step_size_x * (i - (steps_x - 1)/2), 
                    step_size_y * (j - (steps_y - 1)/2), 
                    step_size_z * (k - (steps_z - 1)/2), 
                    average_field[i, j, k, 0], average_field[i, j, k, 1], average_field[i, j, k, 2])
                lines.append(line)
    
    with open(name, "w") as f:
        f.writelines(lines)


def average_fields(mat): 
    """Average the fields in the matrix."""
    return np.mean(mat, axis=0)


def aug_all(mat, target, xy=True, z=False, mut=False):
    full_aug = []
    full_aug_target = []

    for ind, i in enumerate(mat):
        x_aug, y_aug = augment_mat_field(i, target[ind])
        [full_aug.append(j) for j in x_aug]
        [full_aug_target.append(j) for j in y_aug]

    return np.array(full_aug), np.array(full_aug_target)


def augment_mat_field(mat, target, xy=True, z=False):
    aug_target = []
    aug_mat = []

    if xy:
        x_flip = np.array(np.flip(mat, axis=0), dtype=float)
        y_flip = np.array(np.flip(mat, axis=1), dtype=float)
        xy_flip = np.array(np.flip(np.flip(mat, axis=1), axis=0), dtype=float)

        x_flip[:, :, :, 0] = -1 * x_flip[:, :, :, 0]
        y_flip[:, :, :, 1] = -1 * y_flip[:, :, :, 1]
        xy_flip[:, :, :, 0] = -1 * xy_flip[:, :, :, 0]
        xy_flip[:, :, :, 1] = -1 * xy_flip[:, :, :, 1]

        aug_mat.append(mat)
        aug_mat.append(x_flip)
        aug_mat.append(y_flip)
        aug_mat.append(xy_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)

    if z:
        z_flip = np.array(np.flip(mat, axis=2), dtype=float)
        xz_flip = np.array(np.flip(np.flip(mat, axis=0), axis=2), dtype=float)
        yz_flip = np.array(np.flip(np.flip(mat, axis=1), axis=2), dtype=float)
        xyz_flip = np.array(
            np.flip(np.flip(np.flip(mat, axis=2), axis=1), axis=0), dtype=float
        )

        z_flip[:, :, :, 0] = -1 * z_flip[:, :, :, 0]
        xz_flip[:, :, :, 0] = -1 * xz_flip[:, :, :, 0]
        xz_flip[:, :, :, 2] = -1 * xz_flip[:, :, :, 2]
        yz_flip[:, :, :, 1] = -1 * yz_flip[:, :, :, 1]
        yz_flip[:, :, :, 2] = -1 * yz_flip[:, :, :, 2]
        xyz_flip[:, :, :, 0] = -1 * xyz_flip[:, :, :, 0]
        xyz_flip[:, :, :, 1] = -1 * xyz_flip[:, :, :, 1]
        xyz_flip[:, :, :, 2] = -1 * xyz_flip[:, :, :, 2]

        aug_mat.append(z_flip)
        aug_mat.append(xz_flip)
        aug_mat.append(yz_flip)
        aug_mat.append(xyz_flip)

        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)
        aug_target.append(target)

    return aug_mat, aug_target


def pca(
        mat, 
        pca=None, 
        pca_comps=10, 
        verbose=False, 
        write=False,
        bounds={'x': [-3.0, 3.0], 'y': [-3.0, 3.0], 'z': [-3.0, 3.0]} , 
        step_size = {"x": 0.3, "y": 0.3, "z": 0.3}):

    mat_transform = mat.reshape(
        mat.shape[0], mat.shape[1] * mat.shape[2] * mat.shape[3] * mat.shape[4]
    )
    if pca == None:
        pca = PCA(n_components=pca_comps)
        mat_transform = pca.fit_transform(mat_transform)

    else:
        mat_transform = pca.transform(mat_transform)

    cum_explained_var = []
    for i in range(0, len(pca.explained_variance_ratio_)):
        if i == 0:
            cum_explained_var.append(pca.explained_variance_ratio_[i])
        else:
            cum_explained_var.append(
                pca.explained_variance_ratio_[i] + cum_explained_var[i - 1]
            )

    pc0 = pca.components_[0]
    # print(np.shape(pc0))
    pc0 = pc0.reshape(1, mat.shape[1], mat.shape[2], mat.shape[3], mat.shape[4])

    fig = make_subplots(rows=1, cols=1, specs=[[{"type": "cone"}]])
    x, y, z = np.meshgrid(
        np.arange(bounds['x'][0], bounds['x'][1]+step_size['x'], step_size['x']),
        np.arange(bounds['y'][0], bounds['y'][1]+step_size['y'], step_size['y']),
        np.arange(bounds['z'][0], bounds['z'][1]+step_size['z'], step_size['z'])
    )
    
    u = pc0[0][:, :, :, 0].flatten()
    v = pc0[0][:, :, :, 1].flatten()
    w = pc0[0][:, :, :, 2].flatten()

    comp_vect_field = pc0.reshape(
        mat.shape[1], mat.shape[2], mat.shape[3], mat.shape[4]
    )

    u_1, v_1, w_1 = split_and_filter(
        comp_vect_field, cutoff=98, std_mean=True, min_max=False
    )

    vector_scale = 3
    fig.add_trace(
        go.Cone(
            x=x.flatten(),
            y=y.flatten(),
            z=z.flatten(),
            u=u_1,
            v=v_1,
            w=w_1,
            sizeref=vector_scale,
        ),
        row=1,
        col=1,
    )

    if write:
        fig.write_html("./out_pca.html")

    if verbose:
        print("cumulative explained vars ratio: \n" + str(cum_explained_var))
        
    return mat_transform, pca


def unwrap_pca(mat, pca, shape):
    mat = pca.inverse_transform(mat)
    mat = mat.reshape(len(mat), shape[1], shape[2], shape[3], shape[4])
    return mat


def helmholtz_hodge_decomp_approx(
    file="../../data/cpet/efield_cox_1sj21.dat", show=False
):
    Vf = mat_pull(file)
    NX, NY, NZ = Vf[:, :, :, 1].shape
    # print(NX, NY, NZ)
    Vfx = Vf[:, :, :, 0]
    Vfy = Vf[:, :, :, 1]
    Vfz = Vf[:, :, :, 2]

    vx_f = np.fft.fftn(Vfx)
    vy_f = np.fft.fftn(Vfy)
    vz_f = np.fft.fftn(Vfz)

    kx = np.fft.fftfreq(NX).reshape(NX, 1, 1)
    ky = np.fft.fftfreq(NY).reshape(NY, 1)
    kz = np.fft.fftfreq(NZ)
    k2 = kx**2 + ky**2 + kz**2
    k2[0, 0, 0] = 1.0  # to avoid inf. we do not care about the k=0 component

    div_Vf_f = vx_f * kx + vy_f * ky + vz_f * kz  # * 1j
    V_compressive_overk = div_Vf_f / k2
    V_compressive_x = np.fft.ifftn(
        V_compressive_overk * kx
    )  # [:,np.newaxis,np.newaxis])
    V_compressive_y = np.fft.ifftn(V_compressive_overk * ky)
    V_compressive_z = np.fft.ifftn(V_compressive_overk * kz)

    V_solenoidal_x = Vfx - V_compressive_x
    V_solenoidal_y = Vfy - V_compressive_y
    V_solenoidal_z = Vfz - V_compressive_z

    # check if the solenoidal part really divergence-free
    divVs = np.fft.ifftn(
        (
            np.fft.fftn(V_solenoidal_x) * kx
            + np.fft.fftn(V_solenoidal_y) * ky
            + np.fft.fftn(V_solenoidal_z) * kz
        )
        * 1j
        * 2.0
        * np.pi
    )

    # print('div_solenoidal max:', abs(divVs).max())
    # check the power in solenoidal and compressive components
    # print('variance:')
    # print( 'solenoidal x,y,z:', V_solenoidal_x.var(), V_solenoidal_y.var(), V_solenoidal_z.var())
    # print('compressive x,y,z:', V_compressive_x.var(), V_compressive_y.var(), V_compressive_z.var())

    if show:
        # plot one slice of the decomposed field on X-Y plane
        X, Y = np.meshgrid(range(NY), range(NX))
        scale = 1
        plt.figure()
        plt.quiver(
            X,
            Y,
            V_solenoidal_x[:, :, 0] / scale,
            V_solenoidal_y[:, :, 0] / scale,
            norm=True,
        )
        plt.figure()
        plt.quiver(
            X,
            Y,
            V_compressive_x[:, :, 0] / scale,
            V_compressive_y[:, :, 0] / scale,
            norm=True,
        )
        plt.ion()
        plt.show()

    solenoidal = {
        "x": V_solenoidal_x.real,
        "y": V_solenoidal_y.real,
        "z": V_solenoidal_z.real,
    }
    compressize = {
        "x": V_compressive_x.real,
        "y": V_compressive_y.real,
        "z": V_compressive_z.real,
    }
    return solenoidal, compressize


def compress(distance_matrix):
    compressed_dictionary = {}
    affinity = AffinityPropagation(affinity="precomputed")
    affinity.fit(distance_matrix)
    cluster_centers_indices = affinity.cluster_centers_indices_
    labels = list(affinity.labels_)
    n_clusters_ = len(cluster_centers_indices)

    print(f"Estimated number of clusters: {n_clusters_}")
    # get count of a value in a list
    for i in range(n_clusters_):
        compressed_dictionary[i] = {
            "count": str(labels.count(i)),
            "index_center": str(cluster_centers_indices[i]),
        }
    return compressed_dictionary

