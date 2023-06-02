import os, json
import numpy as np 
def main():
    bohr_to_ang_conv = True
    for file in os.listdir(os.path.join(os.path.dirname(__file__), 'cpet')):
        if file.split(".")[-1] == "txt" and "field" in file:
            file_header = file.split(".")[0]
            run_index = file_header.split("-")[0].split("_")[-3]
            run_ind = file_header.split("-")[0].split("_")[-2]
            protein= file_header.split("-")[0].split("_")[-1]
            run = file_header.split("-")[-2][3:]
            frame = file_header.split("-")[-1]
            #print(run_ind, protein, run, frame)
            with open(os.path.join(os.path.dirname(__file__), 'plot_options.json')) as f:
                plot_options = json.load(f)
            # get first line from file 
            with open(os.path.join(os.path.dirname(__file__), 'cpet', file)) as f:
                first_line = f.readline()
            first_line_data = first_line.split()
            
            bohr_to_ang = 1
            if bohr_to_ang_conv:
                bohr_to_ang = 1.88973
            
            center = first_line_data[1]                    
            x_axis = first_line_data[2]
            y_axis = first_line_data[3]
            center = np.array([float(i) for i in  center.split(":")]) * bohr_to_ang
            x_axis = np.array([float(i) for i in  x_axis.split(":")]) * bohr_to_ang
            y_axis = np.array([float(i) for i in  y_axis.split(":")]) * bohr_to_ang
            print("center raw: " + center)
            print("x_axis raw: " + x_axis)
            print("y_axis raw: " + y_axis)
            plot_options["alignment_dict"]["alignment_method"] = "dict"
            plot_options["alignment_dict"]["center"] = list(center)
            plot_options["alignment_dict"]["x_axis"] = list((x_axis-center))
            plot_options["alignment_dict"]["y_axis"] = list((y_axis-center))

            # save as plot_options_{1}_{run_ind}_{protein}_{run}_{frame}.json
            json_name = os.path.join(os.path.dirname(__file__), 'cpet','plot_options_{run_index}_{run_ind}_{protein}_{run}_{frame}.json'.format(run_index=run_index, run_ind=run_ind, protein=protein, run=run, frame=frame))
            #efield_cox_1_01_WT-run3-178.dat
            field_file = os.path.join(os.path.dirname(__file__), "cpet", 'efield_cox_{run_index}_{run_ind}_{protein}-run{run}-{frame}.dat'.format(run_index=run_index, run_ind=run_ind, protein=protein, run=run, frame=frame))
            pdb_file = os.path.join(os.path.dirname(__file__), "frames", '{run_index}_{run_ind}_{protein}-run{run}-{frame}.pdb'.format(run_index=run_index, run_ind=run_ind, protein=protein, run=run, frame=frame))
            save_file = os.path.join(os.path.dirname(__file__), "cpet", "{run_index}_{run_ind}_{protein}-run{run}-{frame}.html".format(run_index=run_index, run_ind=run_ind, protein=protein, run=run, frame=frame))
            plot_options["save_name"]=save_file
            with open(json_name, 'w') as f:
                json.dump(plot_options, f, indent=4)
            # run command to plot 
            os.system("plot_field_on_struct.py -options {} -pdb_file {} -field_file {}".format(json_name, pdb_file, field_file))

main()
