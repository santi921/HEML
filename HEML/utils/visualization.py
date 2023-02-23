import numpy as np
from HEML.utils.data import pull_mats_w_label, pca


def shift_and_rotate(
    xyz_list, center=[0, 0, 0], x_axis=[1, 0, 0], y_axis=[0, 1, 0], z_axis=[0, 0, 1]
):
    for i in range(len(xyz_list)):
        xyz_list[i] = xyz_list[i] - center
        xyz_list[i] = np.array(
            [
                np.dot(xyz_list[i], x_axis),
                np.dot(xyz_list[i], y_axis),
                np.dot(xyz_list[i], z_axis),
            ]
        )
    return xyz_list


def ms(x, y, z, radius, resolution=50):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0 : 2 * np.pi : resolution * 2j, 0 : np.pi : resolution * 1j]
    X = radius * np.cos(u) * np.sin(v) + x
    Y = radius * np.sin(u) * np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)
