import numpy as np 
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from glob import glob
from tqdm import tqdm
from HEML.utils.data import mat_pull
from HEML.utils.visualization import mat_to_cones


def show_in_out_plots(folder_plot):
    #fig = make_subplots(
    #    rows=1, cols=2,
    #    specs=[[{'type': 'surface'}, {'type': 'surface'}]],
    #            subplot_titles=("Scalar In", "Scalar Out"),
    #            )

    frames_in = []
    # get all the files in the folder ending in .dat
    files = glob(folder_plot + '/*.dat')        
    shape = mat_pull(files[0]).shape
    # add 1 to the shape to get the number of slices
    shape = (1, shape[0], shape[1], shape[2], shape[3])
    print("Gathering frames....")
    
    #for i in range(len(files)):
    for i in tqdm(range(int(len(files)/10))):
        if i == 0:
            frame_init = mat_to_cones(
                mat_pull(files[i]),
                shape ,
                vector_scale=1
                )
        frames_in.append(
            mat_to_cones(
                mat_pull(files[i]),
                shape ,
                vector_scale=1
            ))

        
    fig = go.Figure(data=frame_init)
    fig.frames = [
        go.Frame(
            data=[frames_in[i]], 
            name=str(i)) 
            for i in range(len(frames_in))
        ]    
 
    sliders = [
                {
                    "pad": {"b": 10, "t": 60},
                    "len": 0.9, "x": 0.1, "y": 0,
                    "steps": [
                        {
                            "args": [[f.name], frame_args(0)],
                            "label": str(k),
                            "method": "animate",
                        }
                        for k, f in enumerate(fig.frames)
                    ],
                }
            ]

    # Layout
    fig.update_layout(
            scene = dict(
                xaxis = dict(nticks=10, range=[-3,3],),
                yaxis = dict(nticks=10, range=[-3,3],),
                zaxis = dict(nticks=10, range=[-3,3],)),

            title='Slices in volumetric data',
            width=1000, height=1000,
            updatemenus = [
                {
                    "buttons": [
                        {
                            "args": [None, frame_args(50)],
                            "label": "&#9654;", # play symbol
                            "method": "animate",
                        },
                        {
                            "args": [[None], frame_args(0)],
                            "label": "&#9724;", # pause symbol
                            "method": "animate",
                        },
                    ],
                    "direction": "left",
                    "pad": {"r": 10, "t": 70},
                    "type": "buttons",
                    "x": 0.1,
                    "y": 0,
                }
            ],
            sliders=sliders
    )
    fig.show()
    


def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True,
            "transition": {"duration": duration, "easing": "linear"},
        }




