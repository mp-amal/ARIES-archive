import os
import pandas as pd
import shutil


filestocopy_csv = pd.read_csv('/home/archive/Documents/ARIES-archive/scripts/extra_features/advo_file_names.csv')
print(filestocopy_csv['filename'])


for filetocopy in filestocopy_csv['filename']:
    for path, dirs, files in os.walk('/data/archived_data/astro_data/final_data'):
        for file in files:
            if filetocopy == file :
                shutil.copy(os.path.join(path,file),os.path.join('/data/d_share/pranshu_advo',file))
                print('copied')
                




