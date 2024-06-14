import os
import shutil

def mass_migration(source_folder, destination_folder):
    # Create destination folder if it doesn't exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
    
    # Iterate through files in the source folder
    for filename in os.listdir(source_folder):
        # Check if the file ends with '.dat'
        if filename.endswith('.dat'):
            # Construct source and destination paths
            source_path = os.path.join(source_folder, filename)
            folder_name = source_folder.split('/')[-1]  # Extract folder name
            destination_path = os.path.join(destination_folder, f"{folder_name}_{filename}")
            
            # Move the file
            shutil.move(source_path, destination_path)
            print(f"Moved {source_path} to {destination_path}")

# Example usage:
source_folder = './3vxi'
destination_folder = './'
mass_migration(source_folder, destination_folder)
