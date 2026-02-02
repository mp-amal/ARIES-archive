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

year = "2025"
cycle = "C2"
server_ip = "000.00.0.00"  # Replace with your server's IP
username = "username"    # Replace with your username
password = "password"    # Replace with your password

# -----------------------------------------------------------------------------
# NOTE:
# The paths below are DUMMY placeholders for sharing / public repos.
# Replace them with the real deployment paths on your archive server.
# -----------------------------------------------------------------------------

DIR_DOC = '/path/to/ADA_PROGRAM/config/'    # Path to the main config file directory (dummy)

file_keywords = DIR_DOC + 'MasterHeaderList.dat'  # Contains header details (dummy path used)
keywords_df = pd.read_csv(file_keywords, sep=r'\s+', comment='#').set_index('Header')  # MasterHeader file
keywords_df = keywords_df.fillna('NULL')

# Logs (dummy paths)
log_path = "/path/to/DOT_LOG/adfosc_log"
log_list = "/path/to/ADA_PROGRAM/output/DOT/ADFOSC/log_list"

# Data paths (dummy base paths; structure preserved)
remote_path  = f"/path/to/observation_Data/ADFOSC/{year}-{cycle}"
local_path   = f"/path/to/DOT/ADFOSC/rawdata/{year}-{cycle}"
process_path = f"/path/to/DOT/ADFOSC/Processed_Data/{year}-{cycle}"

file_telescopes = DIR_DOC + 'TelescopeList.csv'   # Contains telescope details (dummy path used)

telescope_df = pd.read_csv(file_telescopes, sep=',', comment='#').set_index('ShortName')    # DataFrame of telsecope details
#  Helping funtions 

def process_fits_file(file_path, local_path):
    try:
        # Get PROPNO from external function
        propid = PROPNO(file_path)
        if propid is None:
            propid = ["PXX"]

        # Add PROPNO and ORIGFILE to header
        with fits.open(file_path, mode='update') as hdulist:
            header = hdulist[0].header
            header["PROPNO"] = propid[0].replace('_', '')
            header["ORIGFILE"] = os.path.basename(file_path)

    #     # Open again to check if HISTORY keyword is present
        with fits.open(file_path, ignore_missing_simple=True) as hdul:
            relative_path = os.path.relpath(file_path, start=local_path)
            header = hdul[0].header

            if "HISTORY" in header:
                # print("\n", relative_path,"File header already upated")
                mode = header.get("MODE", "-") or "-"
                type_ = header.get("TYPE", "-") or "-"
                category = header.get("CATEGORY", "-") or "-"
                print('log updated')
                print(file_path)
                with open("/data/DOT/ADFOSC/rawdata/file_track.dat", "a") as f:
                    f.write(f'{os.path.basename(file_path):<30} {mode:<15} {type_:<10} {category:<15} {relative_path}\n')
            else:
                print(f'🚀 Updating : {relative_path}')
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

def check_log_files(log_path,folder_to_process):
    print("check_folder debug :",folder_to_process)
    folder_to_process = datetime.strptime(folder_to_process, "%Y%m%d")
    this_day_fmt = folder_to_process.strftime("%Y_%m_%d")
    prev_day = folder_to_process - timedelta(days=1)
    prev_day_fmt = prev_day.strftime("%Y_%m_%d")

    dates = [this_day_fmt,prev_day_fmt]
    
    for date in dates:
        print(date)
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
            c = 'leave'            
        else:
            c = 'copy'
        return c



        #     if len(files_sorted):
        #         log_count.append(len(files_sorted))
        # print(log_count)
        # if not log_count:
        #     return False
        
            # for f in files_sorted:
            #     match = re.search(file_pattern, f)
            #     if match:
            #         hour = match.group(1)
            #         hourly_counts[hour] += 1
            
            # if date == this_day_fmt:
            #     hours_range = range(9)
            # else:
            #     hours_range = range(17, 24)
            
            # expected_count = 12 if "adfosc_telnet_log_" in file_prefix else 6
            
            # for hour in hours_range:
            #     hour_str = f"{hour:02d}"
            #     count = hourly_counts.get(hour_str, 0)
                
            #     if count < expected_count:
            #         status = f"Hour {hour_str}: Missing {expected_count - count} files"

            #     elif count > expected_count:
            #         status = f"Hour {hour_str}: Extra {count - expected_count} files"
            #     else:
            #         status = f"Hour {hour_str}: OK"
                # print(status)
def copy_subfolders_from_server(remote_path, processing_path, server_ip, username, password):

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname=server_ip, username=username, password=password)

    # List all subdirectories inside remote_path
    # command = f"ls -d {remote_path}/"
    command = f"cd {remote_path} && ls"
    stdin, stdout, stderr = ssh.exec_command(command)
    folder_list = stdout.read().decode().strip().split('\n')
    folder_list = [folder.strip() for folder in folder_list if folder.strip()]



    with SCPClient(ssh.get_transport()) as scp:
        # folder_list = stdout.read().decode().strip().split('\n')
        # folder_list = [folder.strip() for folder in folder_list if folder.strip()]
        copied_folder = []
        for i, folder in enumerate(folder_list, start=1):
            folder_name = os.path.basename(folder.rstrip('/'))
            print(i, folder_name)

            try:
                c = check_log_files(log_path, folder_name)
                if c == 'copy':
                    print('Log available, copying')
                    dest_path = os.path.join(processing_path, folder_name)
                    if not os.path.exists(dest_path):
                        print(f"Copying {folder_name} to {dest_path} ...")
                        scp.get(os.path.join(remote_path, folder_name),
                                dest_path, recursive=True)
                        print(f"{folder_name} 📁✅.")
                        copied_folder.append(folder_name)
                    else:
                        print(f"📁 {folder_name} already exists ✅.")
                else:
                    print('no log')
            except Exception as e:
                print(f"[ERROR] While processing {folder_name}: {e}")
                # continue with next folder
                continue

    return copied_folder










    # try:
    #     # Establish SSH connection
    #     ssh = paramiko.SSHClient()
    #     ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    #     ssh.connect(hostname=server_ip, username=username, password=password)

    #     # List all subdirectories inside remote_path
    #     # command = f"ls -d {remote_path}/"
    #     command = f"cd {remote_path} && ls"
    #     stdin, stdout, stderr = ssh.exec_command(command)
    #     folder_list = stdout.read().decode().strip().split('\n')
    #     folder_list = [folder.strip() for folder in folder_list if folder.strip()]
    #     print(folder_list)
    #     # print(folder_list)
    #     # Create SCP client and copy folders
    #     with SCPClient(ssh.get_transport()) as scp:
    #         i = 0
    #         copied_folder =[]
    #         for folder in folder_list:
    #             print("debug in copy func - folder :", folder)
    #             i = i+1
    #             # if i == 6:
    #             #     break
    #             folder_name = os.path.basename(folder.rstrip('/'))
    #             print(folder_name)
    #             c= check_log_files(log_path,folder_name)
    #             if c == 'copy':
    #                 print('Log available, copying')
    #                 dest_path = os.path.join(processing_path, folder_name)
    #                 if not os.path.exists(dest_path):
    #                     print(f"Copying {folder_name} to {dest_path} ...")
    #                     print(f"🔄 {folder_name} ...")
    #                     scp.get(os.path.join(remote_path,folder_name), dest_path, recursive=True)
    #                     # print(f"Copied {folder_name} successfully.")
    #                     print(f"{folder_name} 📁✅.")
    #                     copied_folder.append(folder_name)
    #                 else:
    #                     print(f"📁 {folder_name}  ✅.")
    #                     # continue
    #             else:
    #                 print('no log')
    #     return copied_folder
    # except Exception as e:
    #     print(f"Error: {e}")
    # finally:
    #     ssh.close()
def is_valid_yyyymmdd(date_str):
    if date_str.isdigit() and len(date_str) == 8:
        try:
            datetime.strptime(date_str, "%Y%m%d")
            return True
        except ValueError:
            return False
    return False
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
def rename_fits_file(filename, db_path):
    """
    Rename a FITS file based on its TYPE, DATE-OBS, TELESCOP, and INSTRUME header keywords.

    Parameters:
    - filename : str : Path to the FITS file.
    - db_path  : str : Directory where the renamed file should be moved.
    """
    with fits.open(filename) as hdul:
        header = hdul[0].header
        type_ = header['TYPE']
        date_obs = header['DATE-OBS']
        telescope = header['TELESCOP']
        instru = header['INSTRUME']
        # propid =header["PROPNO"]
        propid = header.get("PROPNO", "PXX")
        type_code_map = {
            'LAMP': 'L',
            'FLAT': 'F',
            'BIAS': 'B',
            'OBJECT': 'S'
        }

        if type_ in type_code_map:
            code = type_code_map[type_]
            print("TYPE CODE : ",code)
            new_name = f"{code}-{year}AP{propid[0].replace('_','').upper()}-{date_obs}-{telescope}-{instru}.fits"
            new_path = os.path.join(db_path, new_name)
            os.rename(filename, new_path)
            # shutil.move(filename, new_path)
            # print(f"{file} --------------> :  {new_path}")
        # else:
        #     print(f"Unrecognized TYPE '{type_}' in {filename}, skipping.")
def append_astrometryheader(ast_path,new_filename):
    """
    Append Astrometry Keywords to the header of the file 'filename'. Also, compute central RA and DEC
    using those keywords and append them to the header
    Args:
        filename : FITS file to which Astrometry header details have to be appended
    Returns:
        None
    """
    # new_filename = filename.split('.')[0] + '.new'
    headernew = fits.getheader(new_filename, ext=0)
    w = wcs.WCS(new_filename, naxis=2)
    radec = w.wcs_pix2world(headernew['NAXIS1'] / 2, headernew['NAXIS2'] / 2, 1)

    astrometrykeys = ['WCSAXES', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                      'CUNIT1', 'CUNIT2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

    with fits.open(ast_path, mode='update') as hdulist:
        header = hdulist[0].header

        for keyword in astrometrykeys:
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header.append(card=(keyword, headernew[keyword]))

        for idx, keyword in enumerate(['RA', 'DEC']):
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header. append(card=(keyword, str(radec[idx])))
    print("\n"+"<----- .  Astrometry header updated  ----->")
def append_nullheader(filename):
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header
        header.remove('OBJECT', ignore_missing=True, remove_all=True) #neha
        for keyword in keywords_df.index:
            if keyword not in list(header.keys()):
                header.append(card=(keyword, keywords_df.loc[keyword, 'Value']))
    print("<----- 1.  Null Header updated  ----->")

def init_telescope(telescopename):
    """
    Defines a Telescope object in Ephem.
    Args:
        telescopename : Name of the Telescope from which the data was observed
    Returns:
        telescope     : Ephem.Observer object containing Telescope site details
    """
    _, OBS_LONG, OBS_LAT, OBS_ALT, _, _, _ = telescope_df.loc[telescopename].values

    telescope = ephem.Observer()
    telescope.lon = OBS_LONG
    telescope.lat = OBS_LAT
    telescope.elevation = OBS_ALT
    telescope.pressure = 0
    telescope.epoch = ephem.J2000

    return telescope

def get_LST(telescopename, time_obs):
    """
    Calculates Local Sidereal Time for a 'telescope' at time specified in utc by 'time_obs'.
    Args:
        telescopename : Name of the Telescope from which the data was observed
        time_obs      : Time of observation in UTC
    Returns:
        LST_deg       : LST in degrees
    """
    date_obs, time_utc = time_obs.split('T')
    datetime_utc = str(date_obs) + ' ' + str(time_utc)

    telescope = init_telescope(telescopename)
    telescope.date = datetime_utc

    LST_hms = telescope.sidereal_time()
    LST_hms = str(LST_hms).split(':')
    LST_deg = (float(LST_hms[0]) * 1.0 + float(LST_hms[1]) / 60.0 + float(LST_hms[2]) / 3600.0) * 360.0 / 24.0

    return LST_deg
def do_astrometry(TELESCOPE, filename, pixtolerance=0.02):
    """
    Performs Astrometry on images with constraints specified by RA and DEC if UT is specified.
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file on which Astrometry has to be performed
        pixtolerance  : Platescale tolerance in percentage
    Returns:
        None
    """
    print('\n'*3+'\n')
    _, _, OBS_LAT, _, _, _, _ = telescope_df.loc[TELESCOPE].values

    header = fits.getheader(filename, ext=0)
    date_obs = header['DATE-OBS']
    print("\n"+'Doing astrometry...'+"\n")
    try:
        # PIXSCALE = header['PIXSCALE']
        PIXSCALE = 0.200
        PIX_l = str(PIXSCALE - (PIXSCALE * pixtolerance))
        PIX_u = str(PIXSCALE + (PIXSCALE * pixtolerance))

        LST_deg = get_LST(TELESCOPE, date_obs)
        print(LST_deg,date_obs,'LST_deg & date_obs')
        # if telescopename == 'DFOT':
        try:
            subprocess.call('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIX_l + ' --scale-high '
                            + PIX_u + ' --scale-units app --config '+DIR_DOC+'Astrometry.cfg --ra ' + str(LST_deg) + ' --dec ' +
                            str(OBS_LAT) + ' --radius 100 ' + filename, timeout=120, shell=True)
        except subprocess.TimeoutExpired:
            print("\n"+"Astrometry Timed Out (60s) For : "+os.path.basename(filename))
            return False
        else:
            print("\n"+"Astrometry Ran Sucessfully For : "+os.path.basename(filename))
            os.system('rm -rf *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs')
            return True

        # elif telescopename == 'ST':
        #    os.system('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIXSCALE_l + ' --scale-high '
        #              + PIXSCALE_u + ' --scale-units app --config Astrometry.cfg  ' + filename)
        # else:
        #    sys.exit(1)
    except KeyError:
        print("\n"+"ERROR: Header Keyword 'PIXSCALE' not found")
        return False 


def main():
    # check log file and copy if no log will not copy
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(hostname=server_ip, username=username, password=password)

    # List all subdirectories inside remote_path
    # command = f"ls -d {remote_path}/"
    command = f"cd {remote_path} && ls"
    stdin, stdout, stderr = ssh.exec_command(command)
    folder_list = stdout.read().decode().strip().split('\n')
    folder_list = [folder.strip() for folder in folder_list if folder.strip()]

    with SCPClient(ssh.get_transport()) as scp:
            # folder_list = stdout.read().decode().strip().split('\n')
            # folder_list = [folder.strip() for folder in folder_list if folder.strip()]
            copied_folder = []
            for i, folder in enumerate(folder_list, start=1):
                folder_name = os.path.basename(folder.rstrip('/'))
                print(i, folder_name)
                dest_path = os.path.join(local_path, folder_name)
                if os.path.exists(os.path.join(process_path,folder_name)):
                    print("folder exist")
                    continue

                try:
                    c = check_log_files(log_path, folder_name)
                    if c == 'copy':
                        print('Log available, copying')
                        dest_path = os.path.join(local_path, folder_name)
                        if not os.path.exists(dest_path):
                            print(f"Copying {folder_name} to {dest_path} ...")
                            scp.get(os.path.join(remote_path, folder_name),
                                    dest_path, recursive=True)
                            print(f"{folder_name} 📁✅.")
                        else:
                            continue
                            # copied_folder.append(folder_name)

                        if not is_valid_yyyymmdd(folder):
                            if os.path.exists(os.path.join(process_path,folder)):
                                print("Folder processed : ",folder)
                            continue
                        for path,dir,files in os.walk(os.path.join(local_path,folder)):
                            print(os.path.join(local_path,folder))
                            for file in files:
                                # print(file)
                                # if ".fit" in file:
                                if file.lower().endswith((".fit", ".fits")):
                                    file_path = os.path.join(path,file)
                                    # print(file_path)
                                    if not "test" in file_path.lower() and not "ivt" in file_path.lower()and not "point" in file_path.lower() and not "ts_" in file_path.lower():
                                        process_fits_file(file_path, local_path)
                                        # if not os.path.exists(file_path):
                                        #         continue
                                        with fits.open(file_path, ignore_missing_simple=True) as hdul:
                                            relative_path = os.path.relpath(file_path, start=local_path)
                                            header = hdul[0].header
                                            
                                            process_dir_path =os.path.join(process_path,folder)
                                            try:
                                                type = header["TYPE"]
                                            except KeyError:
                                                print(f"Skipping file : 'TYPE' keyword not found in header.")
                                                # return  # or `continue` if this is inside a loop
                                                continue
                                            # type = header["TYPE"]
                                            # print(type)
                                            if type  == "LAMP":
                                                os.makedirs(process_dir_path,exist_ok=True)
                                                rename_fits_file(file_path, process_dir_path)
                                            if type  == "BIAS":
                                                os.makedirs(process_dir_path,exist_ok=True)
                                                rename_fits_file(file_path, process_dir_path)
                                            if type  == "FLAT":
                                                os.makedirs(process_dir_path,exist_ok=True)
                                                rename_fits_file(file_path, process_dir_path)
                                            if type  == "OBJECT":
                                                os.makedirs(process_dir_path,exist_ok=True)
                                                mode = header["MODE"]
                                                if mode == "Spectroscopy":
                                                    os.makedirs(process_dir_path,exist_ok=True)
                                                    rename_fits_file(file_path, process_dir_path)
                                                else:
                                                    print('Mode : ',mode)
                                                    if not "imaging" in file_path:
                                                        os.makedirs(os.path.join(os.path.dirname(file_path),"imaging"), exist_ok=True)
                                                        # print( "\n folder created :",os.path.join(os.path.dirname(file_path),"imaging"))

                                                        image_path = os.path.join(os.path.dirname(file_path),"imaging",os.path.basename(file_path))
                                                        shutil.move(file_path,image_path)
                                                        print("\n file moved : ",os.path.join(os.path.dirname(file_path),"imaging",os.path.basename(file_path)))
                                                        # here do astrometry
                                                        with fits.open(image_path, ignore_missing_simple=True) as hdul:
                                                            header = hdul[0].header
                                                            ra =header['RA']
                                                            dec = header['DEC']
                                                            date_obs = header["DATE-OBS"]
                                                            telescope = header['TELESCOP']
                                                            instru = header['INSTRUME']
                                                            print("RA  : ",ra)
                                                            print("DEC : ",dec)

                                                            try:
                                                                subprocess.call('solve-field --continue --downsample 2 --no-plots --config /home/archive/Documents/program/DOT/Astrometry.cfg --ra ' + str(ra) + ' --dec ' + str(dec) + ' --radius 20 ' +image_path, timeout=60, shell=True)
                                                            except subprocess.TimeoutExpired:
                                                                print("\n"+"Astrometry Timed Out (60s) For : "+file)
                                                                failed_path =os.path.join(os.path.dirname(image_path), "timeout_files")
                                                                # if not os.path.exists(failed_path):
                                                                os.makedirs(failed_path,exist_ok=True)
                                                                shutil.move(image_path,os.path.join(failed_path,file))
                                                                    # return False
                                                            else:
                                                                print("\n"+"Astrometry Ran Sucessfully For : "+file)
                                                            # os.system('rm -rf *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs')
                                                            # return True
                                                                print("Image PAth ",image_path)
                                                                # new_filename = os.path.join(os.path.dirname(image_path),file).split('.')[0] + '.new'

                                                                new_filename = os.path.splitext(image_path)[0] + '.new'
                                                                print("New path ",new_filename)
                                                                image_final_path = os.path.join(process_path,folder)
                                                                os.makedirs(image_final_path,exist_ok=True)
                                                                if os.path.exists(new_filename):
                                                                    append_astrometryheader(image_path,new_filename)
                                                                    print('\n................exist...................\n')
                                                                    # code = 'S'
                                                                    # new_name=code+'-2023APXX-'+date_obs+'-'+telescope+'-'+instru+'.fits'
                                                                    
                                                                    # final_file = os.path.join()
                                                                    with fits.open(new_filename) as hdul:
                                                                        header = hdul[0].header
                                                                        propid =header["PROPNO"] 
                                                                        date_obs = header['DATE-OBS']
                                                                        telescope = header['TELESCOP']
                                                                        instru = header['INSTRUME']
                                                                        new_name = f"S-{year}AP{propid}-{date_obs}-{telescope}-{instru}.fits"
                                                                        new_path = os.path.join(os.path.join(process_path,folder), new_name)
                                                                        os.rename(new_filename, new_path)
                                                                        print(f"Renamed:  {new_path}")
                                                                    # print(header.keys())
                                                                    os.makedirs(image_final_path+"/thumbnails",exist_ok=True)
                                                                    (prefix, sep, suffix) = new_name.rpartition('.')
                                                                    thumb_name=prefix
                                                                    # print(os.path.join(db_path,new_name,'/{}.png'.format(thumb_name)))
                                                                    # import os
                                                                    # print(os.path.isfile(final_path+"/"+dir+'/Imaging'+'/'+new_name))
                                                                    with fits.open(new_path) as hdul:
                                                                        data = hdul[0].data 
                                                                    # Normalize data for better visualization
                                                                    norm = simple_norm(data, 'sqrt', percent=99)

                                                                    # Create a grayscale image
                                                                    plt.figure(figsize=(8, 8))
                                                                    plt.imshow(data, norm=norm, cmap='gray', origin='lower')
                                                                    plt.xlabel('RA (J2000)')
                                                                    plt.ylabel('Dec (J2000)')
                                                                    plt.colorbar(label='Pixel Value')

                                                                    # plt.savefig(final_path+"/"+dir+'/Imaging/thumbnails'+'/{}.png'.format(thumb_name),dpi=60)
                                                                    plt.savefig(image_final_path+"/thumbnails"+'/{}.png'.format(thumb_name), dpi=60)
                                                                    # plt.savefig(os.path.join(db_path,'thumbnails',thumb_name+'.png'))
                                                                    print('Thumbnail saved!!')
                                                                    print("File saved to : ===============================================>>>>>>>>>>>> ",new_path)
                                                                else:
                                                                    print("\nAstrometry faild in '"+os.path.basename(image_path)+"'")
                                                                    if not os.path.exists(os.path.join(os.path.dirname(image_path),'failed_astrometry')):
                                                                        os.makedirs(os.path.join(os.path.dirname(image_path),'failed_astrometry'))
                                                                        shutil.move(image_path,os.path.join(os.path.dirname(image_path),'failed_astrometry'))
                        
                        
                        
                        
                        
                        
                        else:
                            print(f"📁 {folder_name} already exists ✅.")
                    
                    
                    else:
                        print('no log')
                except Exception as e:
                    print(f"[ERROR] While processing {folder_name}: {e}")
                    # continue with next folder
                    continue

        # return copied_folder

main()