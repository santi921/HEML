from math import radians
from turtle import color
from HEML.utils.data import *
from HEML.utils.attrib import *
from HEML.utils.model import *
import plotly.io as pio


def ms(x, y, z, radius, resolution=50):
    """Return the coordinates for plotting a sphere centered at (x,y,z)"""
    u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
    X = radius * np.cos(u)*np.sin(v) + x
    Y = radius * np.sin(u)*np.sin(v) + y
    Z = radius * np.cos(v) + z
    return (X, Y, Z)


def main():

    x, y = pull_mats_w_label('./dat')
    arr_min, arr_max,  = np.min(x), np.max(x)
    x = (x - arr_min) / (arr_max - arr_min + 1e-18)

    y = [np.argmax(i) for i in y]


    x_untransformed = x
    x_pca, pca_obj = pca(x, verbose = True) 
    shape_mat = x.shape
    #pca_obj.components_[0]
    
    fig = make_subplots(
        rows=2, cols=3,


        subplot_titles=(
            "PCA0({:.2f})".format(pca_obj.explained_variance_ratio_[0]), 
            "PCA1({:.2f})".format(pca_obj.explained_variance_ratio_[1]), 
            "PCA2({:.2f})".format(pca_obj.explained_variance_ratio_[2]),
            "PCA3({:.2f})".format(pca_obj.explained_variance_ratio_[3]),
            "PCA4({:.2f})".format(pca_obj.explained_variance_ratio_[4])
            ),
        
        vertical_spacing = 0.1,
        horizontal_spacing = 0.1,
        
        #[{'type': 'surface'}, {'type': 'surface'}, {'type': 'surface'}],
        #[{'type': 'surface'}, {'type': 'surface'}, {'type': 'surface'}]
        specs=[
                [{'type': 'surface'}, {'type': 'surface'}, {'type': 'surface'}],
                [{'type': 'surface'}, {'type': 'surface'}, {'type': 'surface'}],
            ])
        

    for ind,pca_comp in enumerate(pca_obj.components_[:-5]):
        comp_vect_field = pca_comp.reshape(shape_mat[1], shape_mat[2], shape_mat[3], shape_mat[4])

        x, y, z = np.meshgrid(
                        np.arange(-3, 3.3, 0.3),
                        np.arange(-3, 3.3, 0.3),
                        np.arange(-3, 3.3, 0.3)
                        )

        u_1, v_1, w_1 = split_and_filter(
            comp_vect_field, 
            cutoff=95, 
            std_mean=True, 
            min_max=False
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
                sizeref=vector_scale),
                row=int(1+np.floor((ind)/3)), 
                col=1+(ind)%3)
            
        #plot atoms
        dict_atom = {
                    "Fe": [ 0.0, 0.0,  0.0],
                    "N1": [-1.9, 0.0,  0.0],
                    "N2": [ 0.0, 2.1,  0.0],
                    "N3": [ 1.9, 0.0,  0.0],
                    "N4": [ 0.0,-2.1,  0.0],
                    "O1": [ 0.0, 0.0,  1.8]
                    }

        for atom, pos in dict_atom.items():
            if(atom[:2] == "Fe"): 
                radius = 1.40
                color = 'rgb(255, 140, 0)'
            elif(atom[:1] == "N"):
                radius = 0.65
                color = 'rgb(39, 116, 174)'
            elif(atom[:1] == "O"):
                radius = 0.6
                color = 'rgb(255, 0, 0)'
            else: 
                radius = 0.25
                color = 'rgb(255,255, 255)'

            draw_list = ms(pos[0], pos[1], pos[2], radius=radius)

            fig.add_trace(
                go.Surface(
                    x=draw_list[0], 
                    y=draw_list[1], 
                    z=draw_list[2], 
                    opacity=0.999,
                colorscale = [color,color]), 
            row=int(1+np.floor((ind)/3)), 
            col=1+(ind)%3)



        fig.update_layout(yaxis_range=[-3,3], xaxis_range=[-3,3])


    fig.show()
    fig.write_html("../../reporting/pca_components2.html")
    pio.write_image(fig, 'filename.pdf', scale=2, width=800, height=800)

main()