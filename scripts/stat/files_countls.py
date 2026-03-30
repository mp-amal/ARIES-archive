import os
import yaml 
from astropy.io import fits
import pandas as pd

with open("/home/archive/Documents/ARIES-archive/config/paths.yaml", "r") as file:
    paths_yaml = yaml.safe_load(file)


rows =[]
for path, dirs, files in os.walk(paths_yaml['rel_path']):
    for file in files:
        if file.endswith('.png') or file.endswith('.jpg'):
            continue
        if file.startswith("S-"):
            f = os.path.join(path,file)

            try:
                with fits.open(f) as hdul:
                    # print(f"OK   : {os.path.basename(f)}")
                    header = hdul[0].header
                    rows.append({
                        "filename": file,
                        "OBJECT": header.get("OBJECT", ""),
                        "DATE-OBS": header.get("DATE-OBS", "")
                    })

            except Exception as e:
                print(f"BAD  : {os.path.basename(f)} --> {e}")
















#             with fits.open(f) as hdul:
#                 

#                 rows.append({
#                     "filename": file,
#                     "OBJECT": header.get("OBJECT", ""),
#                     "DATE-OBS": header.get("DATE-OBS", "")
#                 })
df = pd.DataFrame(rows)
df.to_csv('./output.csv', index=False)