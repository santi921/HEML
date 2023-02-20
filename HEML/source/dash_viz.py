import dash
import dash_bio as dashbio
from dash import html
from dash_bio.utils import PdbParser, create_mol3d_style, xyz_reader
from dash.dependencies import Input, Output
import numpy as np 
from dash import html, dcc
import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from HEML.utils.data import  get_nodes_and_edges_from_pdb
from HEML.utils.data import *
from HEML.utils.attrib import *
from HEML.utils.model import *

atom_int_dict = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'F': 9,
    'P': 15,
    'S': 16,
    'Cl': 17,
    'Br': 35,
    'Fe': 26, 
    'FE': 26, 
    'I': 53
}

int_atom_dict = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    15: 'P',
    16: 'S',
    17: 'Cl',
    35: 'Br',
    26: 'Fe',
    53: 'I'
}

atomic_size = {
    'H': 0.5,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.80,
    'S': 1.80,
    'Cl': 1.75,
    'Br': 1.85,
    'Fe': 1.80,
    'I': 1.98
}

atom_colors = {
    'H': 'white',
    'C': 'black',
    'N': 'blue',
    'O': 'red',
    'F': 'orange',
    'P': 'green',
    'S': 'yellow',
    'Cl': 'green',
    'Br': 'brown',
    'Fe': 'orange',
    'I': 'purple'
}

app = dash.Dash(__name__)

int_atom_dict = {
    1: 'H',
    6: 'C',
    7: 'N',
    8: 'O',
    9: 'F',
    15: 'P',
    16: 'S',
    17: 'Cl',
    35: 'Br',
    26: 'Fe',
    53: 'I'
}

atom_mass = {
    'H': 1.00794,
    'C': 12.0107,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984,
    'P': 30.9738,
    'S': 32.065,
    'Cl': 35.453,
    'Br': 79.904,
    'Fe': 55.845,
    'I': 126.904
}

list_data = []
bonds_list = []


def get_cones_viz_from_pca(vector_scale = 3, components = 10, data_file = "../../data/protein_data.csv", dir_fields = "../../data/cpet/"): 

    cones = []

    x, _ = pull_mats_w_label(dir_data = data_file, dir_fields = dir_fields)
    arr_min, arr_max,  = np.min(x), np.max(x)
    #x = (x - arr_min) / np.abs(arr_max - arr_min + 0.1)
    # getting sign of every element
    x_sign = np.sign(x)
    # getting absolute value of every element
    x_abs = np.abs(x)
    # applying log1p
    x_log1p = np.log1p(x_abs)
    # getting sign back
    x = np.multiply(x_log1p, x_sign)
    
    x_untransformed = x
    x_pca, pca_obj = pca(x, verbose = True, pca_comps = components, write = False) 
    shape_mat = x.shape


    for ind,pca_comp in enumerate(pca_obj.components_):
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
        
        cones.append(go.Cone(
            x=x.flatten(), 
            y=y.flatten(), 
            z=z.flatten(), 
            u=u_1,
            v=v_1, 
            w=w_1,
            sizeref=vector_scale,
            opacity=0.4, 
            showscale=False,))
        
    return cones 
        



filtered_atom, bonds, filtered_xyz = get_nodes_and_edges_from_pdb('../../data/pdbs_processed/1a4e.pdb')
vector_field_pca = get_cones_viz_from_pca(vector_scale = 5, components = 10)
print(len(filtered_atom))
me = True

for i in range(len(filtered_atom)):
    list_data.append({
        #"serial": i,
        #"name": int_atom_dict[filtered_atom[i]],
        "symbol": int_atom_dict[filtered_atom[i]],
        #"elem": int_atom_dict[filtered_atom[i]], 
        #"atom": int_atom_dict[filtered_atom[i]], 
        #"positions": list(filtered_xyz[i]), 
        "x": filtered_xyz[i][0],
        "y": filtered_xyz[i][1],
        "z": filtered_xyz[i][2],
        #"residue_name": 'ARG', 
        #"residue_index": 1,
        #"chain": "A",
        #"mass_magnitude": atom_mass[int_atom_dict[filtered_atom[i]]]
        })
#print(list_data)
for i in range(len(bonds)): 
    bonds_dict = {"atom1_index": bonds[i][0], "atom2_index": bonds[i][1], "bond_order": 1.0}
    bonds_list.append(bonds_dict)
data = {}
data["atoms"] = list_data


#data["bonds"] = bonds_list
#data = urlreq.urlopen(
#    'https://git.io/speck_methane.xyz'
#).read().decode('utf-8')
#data = xyz_reader.read_xyz(datapath_or_datastring=data, is_datafile=False)
#print(data)

#styles = create_mol3d_style(
#    data["atoms"], 
#    visualization_type='cartoon', 
#    color_element='atom',)
fig = go.Figure(data=[vector_field_pca[7]],
            layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )

app.layout = html.Div([
    dcc.Dropdown(
        id='default-speck-preset-views',
        options=[
            {'label': 'Default', 'value': 'default'},
            {'label': 'Ball and stick', 'value': 'stickball'}
        ],
        value='default'
    ),
    dcc.Graph(figure=fig),
    dashbio.Speck(
        id='default-speck',
        data=list_data,
        view={
                'resolution': 400,
                'ao': 0.1,
                'outline': 1,
                'atomScale': 0.25,
                'relativeAtomScale': 0.33,
                'bonds': True
            }
    ),
])
@app.callback(
    Output('default-speck', 'presetView'),
    Input('default-speck-preset-views', 'value')
)
def update_preset_view(preset_name):
    return preset_name


if __name__ == '__main__':
    app.run_server(
        debug=True, 
        host='127.0.0.1', 
        port=8032)
    