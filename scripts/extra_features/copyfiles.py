import os
import pandas as pd
import shutil


filestocopy_csv = pd.read_csv('/home/archive/Documents/ARIES-archive/scripts/extra_features/file_list_1(1).csv')
# print(filestocopy_csv['filename'])


for filetocopy in filestocopy_csv['filename']:
    for path, dirs, files in os.walk('/data/archived_data/astro_data/final_data'):
        for file in files:
            if filetocopy == file :
                # print(file)
                shutil.copy(os.path.join(path,file),os.path.join('/data/d_share/pranshu_advo',file))
                print('copied')
                




