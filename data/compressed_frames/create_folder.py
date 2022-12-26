import os 

def create_folders():
    """Create a folder for each protein in the current directory and move the
        protein into the folder.    
    """
    # Get the current directory
    current_dir = os.getcwd()
    # Get the list of files in the current directory
    files = os.listdir(current_dir)
    # Loop over the files
    for file in files:
        if file[-3:] == "pdb":
            # Get the name of the protein
            protein = file.split(".")[0].split("_")[1]
            frame = file.split(".")[0].split("_")[0]
            # Create a folder for the protein
            if not os.path.exists(frame + "_" + protein):
                os.mkdir(frame + "_" + protein)
            # Move the protein into the folder
            os.rename(file, frame + "_" + protein + "/" + file)

create_folders()

