import os
import re
import shutil
import random
import yaml
from astropy.io import fits
from datetime import datetime, timedelta, date
from pathlib import Path

TELESCOPE = 'DOT'

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

    if month in [2, 3, 4, 5, 6]:
        return f"{year}-C1"
    elif month in [10, 11, 12]:
        return f"{year}-C2"
    elif month == 1:
        return f"{year - 1}-C2"
    else:
        return None
    
cycle =get_current_cycle_name()
# print(cycle) 

def extract_p_folder(path):
    match = re.search(r'/(P\d+)(/|$)', path)
    if match:
        return match.group(1)
    return None

def get_month_abbreviation(month_value):
    # Dictionary mapping numerical month values to their abbreviations
    month_abbr = {
        "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
        "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
        "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec"
    }

    # Fetch the abbreviation from the dictionary
    return month_abbr.get(month_value, "Invalid month")

with open("/home/archive/Documents/ARIES-archive/config/paths.yaml", "r") as file:
    paths_yaml = yaml.safe_load(file)

with open("/home/archive/Documents/ARIES-archive/config/credentials/db_archive.yaml", "r") as cred:
    cred_yaml = yaml.safe_load(cred)
with open("/home/archive/Documents/ARIES-archive/config/credentials/1_9_archive.yaml", "r") as cred:
    arch_cred = yaml.safe_load(cred)


def dest_path(name,date,f_path):
    random_split = name.split('-')
    # print(random_split)
    # print(date)
    year = random_split[1]
    mon =random_split[3]
    month = get_month_abbreviation(mon)
    # day =random_split[4].split('T')[0]
    # date_folder=os.path.basename(str(date))
    date_folder = date
    telescope =random_split[5] 
    instrument =random_split[6].split('.')[0] 
    db_dest_path = paths_yaml['rel_path']+'/'+str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)
    # print(db_dest_path)
    os.makedirs(db_dest_path, exist_ok=True)
    shutil.copy2(f_path,db_dest_path)



def main():
    # rawlist = os.listdir(f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata')

    rawpath = f'/data/{TELESCOPE}/TANSPEC/{cycle}/rawdata'
    rawlist = [ d for d in os.listdir(rawpath) if re.fullmatch(r"\d{8}", d)]
    final_path = f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/Final_data'

    os.makedirs(final_path, exist_ok=True)
    final_list = os.listdir(final_path)
    for_process = list(set(rawlist) - set(final_list))

    print(for_process)
    for date in for_process:
        print('Start :',date)

        
        source_folder_path = f'/data/{TELESCOPE}/TANSPEC/{cycle}/rawdata/{date}/'
        process_folder_path =f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/raw_processing/{date}/'
        os.makedirs(process_folder_path, exist_ok=True)
        for path,dirs,files in os.walk(source_folder_path):
            for file in files:
                # print(file)
                if file.endswith(".log"):
                    if not os.path.exists(os.path.join(process_folder_path, file)):
                        rel_log_path = path.split("rawdata/")[-1]
                        # os.makedirs(os.path.join(final_path,rel_log_path),exist_ok=True)
                        log_path = Path(os.path.join(final_path,rel_log_path,file))
                        log_path = log_path.parent.parent / log_path.name
                        os.makedirs(os.path.dirname(log_path), exist_ok=True)
                        # print(log_path)
                        shutil.copy(os.path.join(path, file), log_path)
                if file.endswith("Z.fits") and file.upper().startswith("SLOPE-") :
                    # print(path,file)
                    if not os.path.exists(os.path.join(process_folder_path, file)):
                        rel_file_path = os.path.join((path.split("rawdata/")[-1]),file)
                        # print(rel_file_path)
                        os.makedirs(os.path.join(f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/raw_processing/',path.split("rawdata/")[-1]),exist_ok=True)

                        shutil.copy(os.path.join(path, file), os.path.join(f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/raw_processing/', rel_file_path))
                        with fits.open(os.path.join(os.path.join(f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/raw_processing/', rel_file_path)), mode='update') as hdul:
                            header = hdul[0].header
                            # header['ORIGFILE'] = file
                            header.set('ORGFILE', rel_file_path, 'FIle path form Archive', after='FNAME')
                            hdul.flush()

# Copying process done -----------------------------------------------------------------------------------
        for path,dirs,files in os.walk(process_folder_path):
            for file in files:
                print(file)
                file_to_process = os.path.join(path, file)
                if "SLOPE-" in file.upper():
                #     # with fits.open(file_to_process) as hdul:
                    with fits.open(file_to_process, mode='update') as hdul:
                        
                        header = hdul[0].header
                #         # header['ORIGFILE'] = file
                #         # header.set('ORGFILE', file, 'File name received in archive for processing', after='FNAME')
                        Date_obs = header["DATE_OBS"] + "T" + header["TIME_OBS"]
                        dt = datetime.fromisoformat(Date_obs)
                        rounded_ms = round(dt.microsecond / 1000)
                        dt_rounded = dt.replace(microsecond=0) + timedelta(milliseconds=rounded_ms)
                        Date_obs = dt_rounded.isoformat(timespec='milliseconds')
                        dt_object = datetime.fromisoformat(header["DATE_OBS"])
                        hdul.flush()
                #         # 2. Print the year
                        year = str(dt_object.year)
                        object = header["OBJECT"].upper() 
                        file_instrument_pc_path = os.path.join(path,file)
                #         # print(cycle)
                        if cycle in file_instrument_pc_path:
                            file_rel_path = file_instrument_pc_path.split("raw_processing")[1][1:]
                            # print(file_rel_path)
                        pi = extract_p_folder(os.path.join(file_instrument_pc_path,header['FNAME']))
                        # print(pi)
                        if pi == None:
                            # print(file,file_instrument_pc_path)
                            pi = "PXX"

                        if "TEST" in file_to_process.upper()   :
                            name = "TEST-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            print(f_path)
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
                #                 # print(f_path)
                            continue
                        if object == "LAMP":
                #             # print(file)
                            name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # print(lamp_name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            # print(f_path)
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
                                continue
                        if "CONT" in file.upper():
                            name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
                                print("renamed")  
                                continue                          
                        if "AR" in file.upper():
                            name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
                                print("renamed")
                                continue
                        if "NE" in file.upper()   :
                            name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
                                print("renamed")
                                continue
                        
                        # elif "TEST" in final_path.upper()   :
                        #     name = "TEST-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                        #     # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                        #     f_path = os.path.join(final_path, os.path.dirname(file_rel_path),name)
                        #     os.makedirs(os.path.dirname(f_path), exist_ok=True)
                        #     if not os.path.exists(f_path):
                        #         os.rename(file_to_process, f_path)
                        #         dest_path(name,date,f_path)
                        #         # print(f_path)
                        # elif "FLAT" in file.upper() :
                        #     name = "F-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                        #     # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                        #     f_path = os.path.join(final_path, os.path.dirname(file_rel_path),name)
                        #     os.makedirs(os.path.dirname(f_path), exist_ok=True)
                        #     if not os.path.exists(f_path):
                        #         os.rename(file_to_process, f_path)
                        #         dest_path(name,date,f_path)
                        # elif header ["EPADU"]==1.07:
                        #     name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                        #     # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                        #     f_path = os.path.join(final_path, os.path.dirname(file_rel_path),name)
                        #     os.makedirs(os.path.dirname(f_path), exist_ok=True)
                        #     if not os.path.exists(f_path):
                        #         os.rename(file_to_process, f_path)
                        #         dest_path(name,date,f_path)
                        else:
                            name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                            # f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",name)
                            f_path = Path(os.path.join(final_path, os.path.dirname(file_rel_path),name))
                            f_path = f_path.parent.parent / f_path.name
                            os.makedirs(os.path.dirname(f_path), exist_ok=True)
                            if not os.path.exists(f_path):
                                os.rename(file_to_process, f_path)
                                dest_path(name,date,f_path)
        import paramiko
        from scp import SCPClient
        server_ip = arch_cred['ip']
        username = arch_cred['user']
        password = arch_cred['password']

        remote_folder = "/home/archive/data/Data_Share/DOT/TANSPEC"
        local_folder = os.path.join(final_path,date)

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
main()
