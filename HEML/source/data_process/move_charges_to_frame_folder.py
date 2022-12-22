import os 
def move_charges_into_folder(charges_root, compressed_frame_folder):
    """Move the charges into the folder of the frame.
    Takes: 
        charges_root: the root directory of the charges
        compressed_frame_folder: the root directory of the compressed frames
    """
    #get list of frame number and protein name from compressed_frame_folder
    # get name of every folder in compressed_frame_folder
    folders = [f for f in os.listdir(compressed_frame_folder) if os.path.isdir(os.path.join(compressed_frame_folder, f))]
    #suffix = "_movie.pqr"
    #suffix = "_todo_process.pqr"
    suffix = ".pqr"
    for i in folders: 
        frame = i.split("_")[0]
        protein = i.split("_")[1]
        # check if the charges file exists
        #print(i)
        if os.path.exists(os.path.join(charges_root, frame + "_" + protein + suffix)):
            # move the charges file into the folder
            os.system("cp " + os.path.join(charges_root, frame + "_" + protein + suffix) + " " + os.path.join(compressed_frame_folder, i + "/" + i +"_charges.pqr"))
        else: print(i + " not found")

move_charges_into_folder(
    "../charges_2/processed_charges/", 
    "./"
)
#move_charges_into_folder(
#    "../charges_processed", 
#    "./"
#)

