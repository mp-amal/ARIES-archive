TELESCOPE = "DOT"
# Standard Library
import os
import re
import sys
import glob
import time
import json
import shutil
import random
import subprocess
from collections import defaultdict
from datetime import datetime, timedelta, date
from collections import Counter
import ephem
from pathlib import Path
# Third-Party Libraries
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.visualization import simple_norm
# import ephem
# import aplpy
import paramiko
from scp import SCPClient

def get_current_cycle_name():
    """
    Return cycle name based on today's date.

    A cycle = Feb, Mar, Apr, May
    B cycle = Oct, Nov, Dec, and next year's Jan

    Example:
    if today is 2027-01-10 -> returns '2026B'
    """
    today = datetime.today()
    # today  = datetime.strptime("2026-01-06 10:30:00", "%Y-%m-%d %H:%M:%S")
    # print(today)
    year = today.year
    month = today.month

    if month in [2, 3, 4, 5]:
        return f"{year}-C1"
    elif month in [10, 11, 12]:
        return f"{year}-C2"
    elif month == 1:
        return f"{year - 1}-C2"
    else:
        return None
    
cycle =get_current_cycle_name()
print(cycle)

def check_log_files(log_path,folder_to_process):
    # print("check_folder debug :",folder_to_process)
    folder_to_process = datetime.strptime(folder_to_process, "%Y%m%d")
    this_day_fmt = folder_to_process.strftime("%Y_%m_%d")
    prev_day = folder_to_process - timedelta(days=1)
    prev_day_fmt = prev_day.strftime("%Y_%m_%d")

    dates = [this_day_fmt,prev_day_fmt]
    
    for date in dates:
        # print(date)
        patterns = ["adfosc_ics_log_" + date, "tcs_telnet_log_" + date, "adfosc_telnet_log_" + date]
        files = ["adfosc_ics_log_", "tcs_telnet_log_", "adfosc_telnet_log_"]
        log_count = []
        for pattern, file_prefix in zip(patterns, files):
            # print(pattern)
            matching_files = [f for f in os.listdir(log_path) if pattern in f]
            file_pattern = file_prefix + date + r"_(\d{2})_(\d{2})_(\d{2})\.txt"
            
            files_sorted = sorted(matching_files, key=lambda f: re.search(file_pattern, f).groups())
            hourly_counts = defaultdict(int)
            # print(files_sorted)
            if len(files_sorted):
                log_count.append(len(files_sorted))

        if len(log_count)<3:
            c = 'No log'            
        else:
            c = 'go ahed log is there'
        return c


def process_fits_file(file_path,rel_path):
    try:
        # Get PROPNO from external function
        propid = PROPNO(file_path)
        if propid is None:
            propid = ["PXX"]

        # Add PROPNO and ORIGFILE to header
        with fits.open(file_path, mode='update') as hdulist:
            header = hdulist[0].header
            header["PROPNO"] = propid[0].replace('_', '')
            header["ORIGFILE"] = rel_path

    #     # Open again to check if HISTORY keyword is present
        with fits.open(file_path, ignore_missing_simple=True) as hdul:
            # relative_path = os.path.relpath(file_path, start=local_path)
            header = hdul[0].header

            if "HISTORY" in header:
                # print("\n", relative_path,"File header already upated")
                mode = header.get("MODE", "-") or "-"
                type_ = header.get("TYPE", "-") or "-"
                category = header.get("CATEGORY", "-") or "-"
                print('log updated')
                print(file_path)
                with open("./ARIES-archive/logs/file_track.dat", "a") as f:
                    f.write(f'{os.path.basename(file_path):<30} {mode:<15} {type_:<10} {category:<15} {rel_path}\n')
            else:
                print(f'🚀 Updating : {rel_path}')
                cmd = [
                    'python3',
                    '/home/archive/Documents/program/DOT/ADFOSC/ADFOSC-ARIES-main_V1.2/2024_Head_Manage_ULogs_v01-2025.py',
                    '-i', file_path
                ]
                result = subprocess.run(cmd, capture_output=True, text=True)
                # print(result)
    #             if result.returncode == 0:
    #                 with fits.open(file_path, ignore_missing_simple=True) as hdul:
    #                     header = hdul[0].header
    #                 with open("/data/DOT/ADFOSC/rawdata/file_track.dat", "a") as f:
    #                     f.write(f'{os.path.basename(file_path):<30} {header["MODE"]:<15} {header["TYPE"]:<10} {header["CATEGORY"]:<15} {relative_path}\n')
    #             else:
    #                 print(f"\n 🔴 Error occurred with file: {file_path}")
    #                 print(result.stderr)
    #                 move_file_path = os.path.join("/data/DOT/ADFOSC/rawdata/keyword_problem", relative_path)
    #                 move_file_dir = os.path.dirname(move_file_path)
    #                 os.makedirs(move_file_dir, exist_ok=True)
    #                 shutil.move(file_path, move_file_path)
    #                 print("\n FILE moved !!!!\n")
    #                 pass

    except Exception as e:
        print(f"❗ Unexpected error processing {file_path}: {e}")

def PROPNO(file_name):
    pattern = r'_P\d{1,2}'
    last_folder = os.path.basename(os.path.dirname(file_name))
    # proppath = os.path.join(last_folder,os.path.basename(file_name))
    
    match = re.findall(pattern, os.path.basename(file_name), flags=re.IGNORECASE)
    # print(file_name)
    # print("proposal id   - ",match)
    # Display the result
    if match:
        return match
    else:
        print("No Proposal ID found in the file name.")

def main():
    date = "20260309"
    log_path = f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata/{date}/{date}_adfosc_log/'
    c = check_log_files(log_path,date)
    print(c)
    
    process_folder_path = f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata/{date}/'
    for path,dirs,files in os.walk(process_folder_path):
        for file in files:
            if not "adfosc_log" in path:
                # print(path.upper())
                if not 'TEST' in path.upper() and not "TEST" in file.upper():
                    # print(path)

                    file_path = os.path.join(path,file)
                    # print(file_path)
                    rel_path = Path(file_path).relative_to(Path(f"/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata"))
                    print(rel_path)
                    process_fits_file(file_path,str(rel_path))
    #                 # if not os.path.exists(file_path):
    #                 #         continue
    #                 with fits.open(file_path, ignore_missing_simple=True) as hdul:
    #                     relative_path = os.path.relpath(file_path, start=local_path)
    #                     header = hdul[0].header
                        
    #                     process_dir_path =os.path.join(process_path,folder)
    #                     try:
    #                         type = header["TYPE"]
    #                     except KeyError:
    #                         print(f"Skipping file : 'TYPE' keyword not found in header.")
    #                         # return  # or `continue` if this is inside a loop
    #                         continue
    #                     # type = header["TYPE"]
    #                     # print(type)
    #                     if type  == "LAMP":
    #                         os.makedirs(process_dir_path,exist_ok=True)
    #                         rename_fits_file(file_path, process_dir_path)
    #                     if type  == "BIAS":
    #                         os.makedirs(process_dir_path,exist_ok=True)
    #                         rename_fits_file(file_path, process_dir_path)
    #                     if type  == "FLAT":
    #                         os.makedirs(process_dir_path,exist_ok=True)
    #                         rename_fits_file(file_path, process_dir_path)
    #                     if type  == "OBJECT":
    #                         os.makedirs(process_dir_path,exist_ok=True)
    #                         mode = header["MODE"]
    #                         if mode == "Spectroscopy":
    #                             os.makedirs(process_dir_path,exist_ok=True)
    #                             rename_fits_file(file_path, process_dir_path)
    #                         else:
    #                             print('Mode : ',mode)
    #                             if not "imaging" in file_path:
    #                                 os.makedirs(os.path.join(os.path.dirname(file_path),"imaging"), exist_ok=True)
    #                                 # print( "\n folder created :",os.path.join(os.path.dirname(file_path),"imaging"))

    #                                 image_path = os.path.join(os.path.dirname(file_path),"imaging",os.path.basename(file_path))
    #                                 shutil.move(file_path,image_path)
    #                                 print("\n file moved : ",os.path.join(os.path.dirname(file_path),"imaging",os.path.basename(file_path)))
    #                                 # here do astrometry
    #                                 with fits.open(image_path, ignore_missing_simple=True) as hdul:
    #                                     header = hdul[0].header
    #                                     ra =header['RA']
    #                                     dec = header['DEC']
    #                                     date_obs = header["DATE-OBS"]
    #                                     telescope = header['TELESCOP']
    #                                     instru = header['INSTRUME']
    #                                     print("RA  : ",ra)
    #                                     print("DEC : ",dec)

    #                                     try:
    #                                         subprocess.call('solve-field --continue --downsample 2 --no-plots --config /home/archive/Documents/program/DOT/Astrometry.cfg --ra ' + str(ra) + ' --dec ' + str(dec) + ' --radius 20 ' +image_path, timeout=60, shell=True)
    #                                     except subprocess.TimeoutExpired:
    #                                         print("\n"+"Astrometry Timed Out (60s) For : "+file)
    #                                         failed_path =os.path.join(os.path.dirname(image_path), "timeout_files")
    #                                         # if not os.path.exists(failed_path):
    #                                         os.makedirs(failed_path,exist_ok=True)
    #                                         shutil.move(image_path,os.path.join(failed_path,file))
    #                                             # return False
    #                                     else:
    #                                         print("\n"+"Astrometry Ran Sucessfully For : "+file)
    #                                     # os.system('rm -rf *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs')
    #                                     # return True
    #                                         print("Image PAth ",image_path)
    #                                         # new_filename = os.path.join(os.path.dirname(image_path),file).split('.')[0] + '.new'

    #                                         new_filename = os.path.splitext(image_path)[0] + '.new'
    #                                         print("New path ",new_filename)
    #                                         image_final_path = os.path.join(process_path,folder)
    #                                         os.makedirs(image_final_path,exist_ok=True)
    #                                         if os.path.exists(new_filename):
    #                                             append_astrometryheader(image_path,new_filename)
    #                                             print('\n................exist...................\n')
    #                                             # code = 'S'
    #                                             # new_name=code+'-2023APXX-'+date_obs+'-'+telescope+'-'+instru+'.fits'
                                                
    #                                             # final_file = os.path.join()
    #                                             with fits.open(new_filename) as hdul:
    #                                                 header = hdul[0].header
    #                                                 propid =header["PROPNO"] 
    #                                                 date_obs = header['DATE-OBS']
    #                                                 telescope = header['TELESCOP']
    #                                                 instru = header['INSTRUME']
    #                                                 new_name = f"S-{year}AP{propid}-{date_obs}-{telescope}-{instru}.fits"
    #                                                 new_path = os.path.join(os.path.join(process_path,folder), new_name)
    #                                                 os.rename(new_filename, new_path)
    #                                                 print(f"Renamed:  {new_path}")
    #                                             # print(header.keys())
    #                                             os.makedirs(image_final_path+"/thumbnails",exist_ok=True)
    #                                             (prefix, sep, suffix) = new_name.rpartition('.')
    #                                             thumb_name=prefix
    #                                             # print(os.path.join(db_path,new_name,'/{}.png'.format(thumb_name)))
    #                                             # import os
    #                                             # print(os.path.isfile(final_path+"/"+dir+'/Imaging'+'/'+new_name))
    #                                             with fits.open(new_path) as hdul:
    #                                                 data = hdul[0].data 
    #                                             # Normalize data for better visualization
    #                                             norm = simple_norm(data, 'sqrt', percent=99)

    #                                             # Create a grayscale image
    #                                             plt.figure(figsize=(8, 8))
    #                                             plt.imshow(data, norm=norm, cmap='gray', origin='lower')
    #                                             plt.xlabel('RA (J2000)')
    #                                             plt.ylabel('Dec (J2000)')
    #                                             plt.colorbar(label='Pixel Value')

    #                                             # plt.savefig(final_path+"/"+dir+'/Imaging/thumbnails'+'/{}.png'.format(thumb_name),dpi=60)
    #                                             plt.savefig(image_final_path+"/thumbnails"+'/{}.png'.format(thumb_name), dpi=60)
    #                                             # plt.savefig(os.path.join(db_path,'thumbnails',thumb_name+'.png'))
    #                                             print('Thumbnail saved!!')
    #                                             print("File saved to : ===============================================>>>>>>>>>>>> ",new_path)
    #                                         else:
    #                                             print("\nAstrometry faild in '"+os.path.basename(image_path)+"'")
    #                                             if not os.path.exists(os.path.join(os.path.dirname(image_path),'failed_astrometry')):
    #                                                 os.makedirs(os.path.join(os.path.dirname(image_path),'failed_astrometry'))
    #                                                 shutil.move(image_path,os.path.join(os.path.dirname(image_path),'failed_astrometry'))
                        
 






main()