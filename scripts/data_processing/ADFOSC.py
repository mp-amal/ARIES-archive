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
import yaml
from datetime import datetime, timedelta, timezone, date
# Combined legacy import line kept as-is intentionally
import os,sys,psycopg2,glob
import shutil

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

with open("/home/archive/Documents/ARIES-archive/config/paths.yaml", "r") as file:
    paths_yaml = yaml.safe_load(file)

with open("/home/archive/Documents/ARIES-archive/config/credentials/db_archive.yaml", "r") as cred:
    cred_yaml = yaml.safe_load(cred)
with open("/home/archive/Documents/ARIES-archive/config/credentials/1_9_archive.yaml", "r") as cred:
    arch_cred = yaml.safe_load(cred)
    
# print(cred_yaml['user'])      

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


def process_fits_file(file_path,rel_path,log_path):
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
                    '-i', file_path, '-l', log_path
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

def rename_fits_file(filename, db_path):

    """
    Rename a FITS file based on its TYPE, DATE-OBS, TELESCOP, and INSTRUME header keywords.

    Parameters:
    - filename : str : Path to the FITS file.
    - db_path  : str : Directory where the renamed file should be moved.
    """
    if "TEST" in filename.upper():
            new_name = f"T-{year}AP{propid[0].replace('_','').upper()}-{date_obs}-{telescope}-{instru}.fits"
            new_path = os.path.join(db_path, new_name)
            os.rename(filename, new_path)
    else:
        pass

    with fits.open(filename) as hdul:
        header    = hdul[0].header
        type_     = header['TYPE']
        date_obs  = header['DATE-OBS']
        year      =  datetime.fromisoformat(date_obs).year
        telescope = header['TELESCOP']
        instru    = header['INSTRUME']
        mode      = header['MODE']
        # propid =header["PROPNO"]
        propid    = header.get("PROPNO", "PXX")
        type_code_map = {
            'LAMP': 'L',
            'FLAT': 'F',
            'BIAS': 'B',
            'OBJECT': 'S'
        }

        if type_ in type_code_map:
            code = type_code_map[type_]
            print("TYPE CODE : ",code)
            if type ==['OBJECT'] and mode =='Spectroscopy':
                code = 'SP'
                new_name = f"{code}-{year}AP{propid[0].replace('_','').upper()}-{date_obs}-{telescope}-{instru}.fits"
            else:
                new_name = f"{code}-{year}AP{propid[0].replace('_','').upper()}-{date_obs}-{telescope}-{instru}.fits"
            new_path = os.path.join(db_path, new_name)
            os.rename(filename, new_path)

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

def get_month_abbreviation(month_value):
    # Dictionary mapping numerical month values to their abbreviations
    month_abbr = {
        "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
        "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
        "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec"
    }

    # Fetch the abbreviation from the dictionary
    return month_abbr.get(month_value, "Invalid month")


# Create file with header if not exists
ad ='/home/archive/Documents/ARIES-archive/logs/adfosc_processed.dat'
if not os.path.exists(ad):
    with open(ad, "w") as f:
        f.write("date\n")

# Read existing dates
with open(ad, "r") as f:
    lines = [line.strip() for line in f.readlines()]

# print(lines)




def main():
    # rawlist = os.listdir(f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata')

    rawpath = f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata'
    rawlist = [ d for d in os.listdir(rawpath) if re.fullmatch(r"\d{8}", d)]
    final_list = os.listdir(f'/data/{TELESCOPE}/ADFOSC/{cycle}/Processed_Data/Final_data')
    for_process = list(set(rawlist) - set(final_list))

    print(for_process)
    for date in for_process:
        print('Start :',date)

        log_path = f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata/{date}/{date}_adfosc_log/'
        c = check_log_files(log_path,date)
        print(c)
        
        source_folder_path = f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata/{date}/'
        process_folder_path =f'/data/{TELESCOPE}/ADFOSC/{cycle}/Processed_Data/raw_processing/{date}/'

        if not os.path.exists(str(process_fits_file)):
            shutil.copytree(source_folder_path, process_folder_path, dirs_exist_ok=True)
        # else:
            pass
        process_dir_path =os.path.join(Path(f"/data/{TELESCOPE}/ADFOSC/{cycle}/Processed_Data/Final_data"),date)
        for path,dirs,files in os.walk(process_folder_path):
            for file in files:
                # print(file)
                if'.fit' in file :
                    file_path = os.path.join(path,file)
                    # if not 'TEST' in file_path.upper():
                    process_path = Path(f"/data/{TELESCOPE}/ADFOSC/{cycle}/Processed_Data/raw_processing")
                    rel_path = Path(file_path).relative_to(process_path)
                    print(rel_path)
                    process_fits_file(file_path,str(rel_path),log_path)
    # # #                 # if not os.path.exists(file_path):
    # # #                 #         continue
                    with fits.open(file_path, ignore_missing_simple=True) as hdul:
                        # relative_path = os.path.relpath(file_path, start=local_path)
                        header = hdul[0].header
                        # print(header['TYPE'])
                        
                        process_dir_path =os.path.join(Path(f"/data/{TELESCOPE}/ADFOSC/{cycle}/Processed_Data/Final_data"),date)
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
                                    img_prs_path = os.path.join(process_folder_path,"imaging")
                                    os.makedirs(img_prs_path, exist_ok=True)
                                    # print( "\n folder created :",os.path.join(os.path.dirname(file_path),"imaging"))

                                    image_path = os.path.join(img_prs_path,os.path.basename(file_path))
                                    shutil.move(file_path,image_path)
                                    print("\n file moved : ",os.path.join(os.path.dirname(file_path),"imaging",os.path.basename(file_path)))
    #                                 # here do astrometry
                                    with fits.open(image_path, ignore_missing_simple=True) as hdul:
                                        header = hdul[0].header
                                        ra =header['RA']
                                        dec = header['DEC']
                                        date_obs = header["DATE-OBS"]
                                        year      =  datetime.fromisoformat(date_obs).year
                                        telescope = header['TELESCOP']
                                        instru = header['INSTRUME']
                                        print("RA  : ",ra)
                                        print("DEC : ",dec)
                                        propid    = header.get("PROPNO", "PXX")

                                        try:
                                            subprocess.call('solve-field --continue --downsample 2 --no-plots --config /home/archive/Documents/ARIES-archive/config/Astrometry.cfg --ra ' + str(ra) + ' --dec ' + str(dec) + ' --radius 20 ' +image_path, timeout=60, shell=True)
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
                                            # image_final_path = os.path.join(process_dir_path,folder)
                                            os.makedirs(process_dir_path,exist_ok=True)
                                            if os.path.exists(new_filename):
                                                append_astrometryheader(image_path,new_filename)
                                                print('\n................exist...................\n')
                                                code = 'S'
                                                # new_name=code+f'-{year}AP{propid}-'+date_obs+'-'+telescope+'-'+instru+'.fits'
                                                
                                                # final_file = os.path.join()
                                                with fits.open(new_filename) as hdul:
                                                    header = hdul[0].header
                                                    propid =header["PROPNO"] 
                                                    date_obs = header['DATE-OBS']
                                                    telescope = header['TELESCOP']
                                                    instru = header['INSTRUME']
                                                    if "TEST" in new_filename.upper():
                                                        new_name = f"T-{year}AP{propid}-{date_obs}-{telescope}-{instru}.fits"
                                                        new_path = os.path.join( process_dir_path, new_name)
                                                        os.rename(new_filename, new_path)
                                                    else:
                                                        new_name = f"S-{year}AP{propid}-{date_obs}-{telescope}-{instru}.fits"
                                                        new_path = os.path.join( process_dir_path, new_name)
                                                        os.rename(new_filename, new_path)
                                                    print(f"Renamed:  {new_path}")
    #                                             # print(header.keys())
                                                os.makedirs(process_dir_path+"/thumbnails",exist_ok=True)
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
                                                plt.savefig(process_dir_path+"/thumbnails"+'/{}.png'.format(thumb_name), dpi=60)
                                                # plt.savefig(os.path.join(db_path,'thumbnails',thumb_name+'.png'))
                                                print('Thumbnail saved!!')
                                                print("File saved to : ===============================================>>>>>>>>>>>> ",new_path)
                                            else:
                                                print("\nAstrometry faild in '"+os.path.basename(image_path)+"'")
                                                if not os.path.exists(os.path.join(os.path.dirname(image_path),'failed_astrometry')):
                                                    os.makedirs(os.path.join(os.path.dirname(image_path),'failed_astrometry'))
                                                    shutil.move(image_path,os.path.join(os.path.dirname(image_path),'failed_astrometry'))
        if os.path.exists(process_dir_path):                 
            fit_files = [f for f in os.listdir(process_dir_path) if f.lower().endswith('.fits')]
            print(fit_files)

            random_file = random.choice(fit_files)
            # print(random_file)
            random_split = random_file.split('-')

            year = random_split[2]
            mon =random_split[3]
            month =get_month_abbreviation(mon)
            day =random_split[4].split('T')[0]
            date_folder=os.path.basename(str(date))
            telescope =random_split[5] 
            instrument =random_split[6].split('.')[0] 
            db_dest_path = paths_yaml['rel_path']+'/'+str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)
            # db_paths.append(db_dest_path)
            thumpnail_path = paths_yaml['thump_path'] +str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)+'/thumbnails'
            if os.path.exists(db_dest_path):
                print('Destination exist with path : '+db_dest_path)
            else:
                shutil.copytree(process_dir_path,db_dest_path)
                print(process_dir_path+ "--------------> copied to db")
                if os.path.exists(process_dir_path+'/thumbnails'):
                    shutil.copytree(process_dir_path+'/thumbnails',thumpnail_path)               
# remove test files from archvie folder
            for file in os.listdir(db_dest_path):
                file_path = os.path.join(db_dest_path, file)

                if os.path.isfile(file_path) and file.startswith("T"):
                    os.remove(file_path)
                    print(f"Deleted: {file}")

            for file in os.listdir(thumpnail_path):
                file_path = os.path.join(thumpnail_path, file)

                if os.path.isfile(file_path) and file.startswith("T"):
                    os.remove(file_path)
                    print(f"Deleted: {file}")
# -------------------------------------------------------

            if date not in lines[1:]:
                with open(ad, "a") as f:
                    f.write(f"{date}\n")

            folder_to_share = os.path.join(paths_yaml['adfosc_final'],date)

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
                            

                            shutil.copy(file_to_share, user_folder_structure)
            
            import paramiko
            from scp import SCPClient

            server_ip = arch_cred['ip']
            username = arch_cred['user']
            password = arch_cred['password']

            remote_folder = "/home/archive/data/Data_Share/DOT"
            local_folder = os.path.join('/data/DOT/ADFOSC/2026-C1/Processed_Data/data_share/DOT',date)

            # =========================
            # CREATE SSH CONNECTION
            # =========================
            def create_ssh_client(ip, user, passwd):
                ssh = paramiko.SSHClient()
                ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                ssh.connect(hostname=ip, username=user, password=passwd)
                return ssh

            # =========================
            # MAIN
            # =========================
            try:
                ssh = create_ssh_client(server_ip, username, password)

                # create remote folder if not exists
                ssh.exec_command(f'mkdir -p "{remote_folder}"')

                # copy folder
                with SCPClient(ssh.get_transport()) as scp:
                    scp.put(local_folder, remote_path=remote_folder, recursive=True)

                print("Folder copied successfully.")

                ssh.close()

            except Exception as e:
                print("Error:", e)

# uncommet form here to above to copy files to 1.9

        # file_path = (os.path.relpath(db_dest_path, paths_yaml['rel_path'])) # Relative file path to be inserted in the database

        # nofp=0 # No of files processed
        # dataarray=[]
        # if os.path.exists(file_path + "/"+"log.txt"):
        #     print ("Directory seems to be processed. If you want to reprocess the directory, then delete log.txt file and then run this program ")
        #     exit()
        # infile=open(str(db_dest_path) + "/"+"log.txt","w")

        # #Total Number of Fits file in the directory
        # totfits = len(glob.glob1(db_dest_path+"/","*.fits"))
        # #Master Log file located just outside of DataArchive folder structure
        # masterlog = open(os.path.join("/data/archived_data/astro_data/final_data/", "masterlog.txt"), 'a')

        # #infile.write('  {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t \n'.format("OBSTimestamp","Observer","Object","RA","DEC","Filter","Filter1","Filter2","Instrument","FileName","Full Path","PROPID"))

        # #Database connectivity
        # connection = None
        # try:
        #     connection = psycopg2.connect(user=cred_yaml['user'],
        #                                 password=cred_yaml['password'],
        #                                 host=cred_yaml['host'],
        #                                 port=cred_yaml['port'],
        #                                 database=cred_yaml['database'])
        #     cursor = connection.cursor()
        # except Exception as error:
        #         infile.write(str(error))
        #         print (error)
        # for fname in os.listdir(db_dest_path):
        #     if os.path.splitext(fname)[1]!=".fits":continue
        #     fpath=os.path.join(db_dest_path,fname)
        #     #if os.path.isdir(fpath):continue
        #     hdulist=fits.open("%s"%fpath)
            
        #     #hdulist.verify('silentfix')
        #     try:
                
        #         OBSDate,OBSTime = str(hdulist[0].header['DATE-OBS']).split("T")
        #         DATEOBS = OBSDate +" " + OBSTime
                
        #         if "Undefined" in str(hdulist[0].header['OBJECT']):
        #             hdulist[0].header["OBJECT"]="NULL"
        #         OBJECT = hdulist[0].header["OBJECT"]
                
        #         if ("Undefined" in str(hdulist[0].header['RA'])) or (hdulist[0].header["RA"]=="NULL"):
        #             hdulist[0].header["RA"]= 0.000
        #         RA = hdulist[0].header["RA"]
                
        #         if "Undefined" in str(hdulist[0].header['DEC']) or (hdulist[0].header["DEC"]== "NULL"):
        #             hdulist[0].header["DEC"]=0.000
        #         DEC = hdulist[0].header["DEC"]
                    
        #         if "Undefined" in str(hdulist[0].header['PROPNO']):
        #             hdulist[0].header["PROPNO"]="NULL"	
        #         PROPNO = 	hdulist[0].header["PROPNO"]
                
        #         if "Undefined" in str(hdulist[0].header['PROG']):
        #             hdulist[0].header["PROG"]="NULL"	
        #         PROG = 	hdulist[0].header["PROG"]
                
        #         if "Undefined" in str(hdulist[0].header['OBSERVER']):
        #             hdulist[0].header["OBSERVER"]="NULL"	
        #         OBSERVER = 	hdulist[0].header["OBSERVER"]
                
        #         if "Undefined" in str(hdulist[0].header['TELESCOP']):
        #             hdulist[0].header["TELESCOP"]="NULL"	
        #         TELESCOP = 	hdulist[0].header["TELESCOP"]
        #         if "Undefined" in str(hdulist[0].header['INSTRUME']):
        #             hdulist[0].header["INSTRUME"]="NULL"	
        #         INSTRUMENT = 	hdulist[0].header["INSTRUME"]
                
        #         if "Undefined" in str(hdulist[0].header['CATEGORY']):
        #             hdulist[0].header["CATEGORY"]="NULL"	
        #         CATEGORY = 	hdulist[0].header["CATEGORY"]
                
        #         if "Undefined" in str(hdulist[0].header['TYPE']):
        #             hdulist[0].header["TYPE"]="NULL"	
        #         TYPE = 	hdulist[0].header["TYPE"]
                
        #         if "Undefined" in str(hdulist[0].header['MODE']):
        #             hdulist[0].header["MODE"]="NULL"	
        #         MODE = 	hdulist[0].header["MODE"]
                
        #         DATAID = fname
                
        #         if "Undefined" in str(hdulist[0].header['ORIGFILE']):
        #             hdulist[0].header["ORIGFILE"]="NULL"	
        #         ORIGFILE = 	hdulist[0].header["ORIGFILE"]
                
                
        #         RELEASE = "NULL"
                
                
        #         if "Undefined" in str(hdulist[0].header['EXPTIME']):
        #             hdulist[0].header["EXPTIME"]="NULL"	
        #         EXPOSURE = 	hdulist[0].header["EXPTIME"]
                
        #         if "Undefined" in str(hdulist[0].header['FILTER1']):
        #             hdulist[0].header["FILTER1"]="NULL"	
        #         FILTER1 = 	hdulist[0].header["FILTER1"]
                
        #         if "Undefined" in str(hdulist[0].header['FILTER2']):
        #             hdulist[0].header["FILTER2"]="NULL"	
        #         FILTER2 = 	hdulist[0].header["FILTER2"]
                
        #         if "Undefined" in str(hdulist[0].header['GRISM']):
        #             hdulist[0].header["GRISM"]="NULL"
        #         GRISM = hdulist[0].header["GRISM"]
                
        #         if "Undefined" in str(hdulist[0].header['SLIT']):
        #             hdulist[0].header["SLIT"]="NULL"	
        #         SLIT = 	hdulist[0].header["SLIT"]
                
        #         if "Undefined" in str(hdulist[0].header['AIRMASS']):
        #             hdulist[0].header["AIRMASS"]="NULL"	
        #         AIRMASS = 	round(float(hdulist[0].header["AIRMASS"]),3)
                
        #         FILESIZE = os.path.getsize(fpath)
                
        #         FILEPATH = file_path
                
        #         TIMEINSERTED = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%s")				
        #         # print(TIMEINSERTED)
                        
        #         datarow=(DATEOBS,OBJECT,RA,DEC,PROPNO,PROG,OBSERVER,TELESCOP,INSTRUMENT,CATEGORY,TYPE,MODE,DATAID,ORIGFILE,EXPOSURE,FILTER1,FILTER2,GRISM,SLIT,AIRMASS,FILESIZE,FILEPATH,TIMEINSERTED)
        #         # print(datarow)
                        
        #         sql_insert='''INSERT INTO astronomy(DATEOBS,OBJECT,RA,DEC,PROPNO,PROG,OBSERVER,TELESCOP,INSTRUMENT,CATEGORY,TYPE,MODE,DATAID,ORIGFILE,EXPOSURE,FILTER1,FILTER2,GRISM,SLIT,AIRMASS,FILESIZE,FILEPATH,TIMEINSERTED) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);'''

        #         cursor.execute(sql_insert, (datarow))
                
        #         dataarray.append(datarow)
        #         infile.write(str(DATAID) + "\n")
        #     except Exception as error:
        #         infile.write(str(error))
        #         print (error)

        #         connection.rollback()
        #     else:
        #         connection.commit()

        #     hdulist.close()
                
        #     nofp=nofp+1
        # #print (dataarray)
        # infile.write("Total No of fits files  = %d \n" %totfits)
        # infile.write("Total No of files processed = %d" %nofp)  
        # print ("Total number of fits files present in the directory = %d" %totfits)
        # print ("Total No of files processed = %d" %nofp)
        # TIMEINSERTED = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%s")
        # masterlog.write('{}\t {}\t {}\t {}\t \n'.format(TIMEINSERTED,db_dest_path,nofp,totfits))
        # if connection:
        #     cursor.close()
        #     connection.close


            
        # infile.close()




main()

