import yaml
import os
from astropy.io import fits

with open("/home/archive/Documents/ARIES-archive/config/paths.yaml", "r") as file:
    paths_yaml = yaml.safe_load(file)


adfosc_final = paths_yaml['adfosc_final']
# folders =os.listdir(adfosc_final)

# for folder in folders:

folder ='20260316'
folder_to_share = os.path.join(paths_yaml['adfosc_final'],folder)

for path, dirs, files in os.walk(folder_to_share):
    # print(path)
    for file in files:
        if file.endswith('.fits'):
            # print(file)
            file_to_share = os.path.join(path,file)
            with fits.open(file_to_share) as hdul:
                header = hdul[0].header
                path_to_move =os.path.dirname(header['ORIGFILE'])
                # print(path_to_move)
                # print(file_to_share)
                user_folder_structure = os.path.join(paths_yaml['adfosc_data_share'],path_to_move,file)
                os.makedirs(os.path.dirname(user_folder_structure),exist_ok=True)
                print(user_folder_structure)
                import shutil

                shutil.copy(file_to_share, user_folder_structure)