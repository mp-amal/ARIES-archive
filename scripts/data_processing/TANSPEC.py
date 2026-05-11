import os
import re
import shutil
from astropy.io import fits
from datetime import datetime, timedelta, date


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

    if month in [2, 3, 4, 5]:
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


def main():
    # rawlist = os.listdir(f'/data/{TELESCOPE}/ADFOSC/{cycle}/rawdata')

    rawpath = f'/data/{TELESCOPE}/TANSPEC/{cycle}/rawdata'
    rawlist = [ d for d in os.listdir(rawpath) if re.fullmatch(r"\d{8}", d)]
    final_path = f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/Final_data'

    os.makedirs(final_path, exist_ok=True)
    final_list = os.listdir(final_path)
    for_process = list(set(rawlist) - set(final_list))

    print(for_process)
    for date in rawlist:
        print('Start :',date)

        
        source_folder_path = f'/data/{TELESCOPE}/TANSPEC/{cycle}/rawdata/{date}/'
        process_folder_path =f'/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/raw_processing/{date}/'
        os.makedirs(process_folder_path, exist_ok=True)
        for path,dirs,files in os.walk(source_folder_path):
            for file in files:
                # print(file)
                if file.endswith("Z.fits") :
                    # print(path,file)
                    
                    if not os.path.exists(os.path.join(process_folder_path, file)):
                        # shutil.copy(os.path.join(path, file), os.path.join(process_folder_path, file))
                        # print("copied")
                        pass
                        # print(file_to_process)
# Copying process done -----------------------------------------------------------------------------------
        for path,dirs,files in os.walk(process_folder_path):
            for file in files:
                    # print(file)
                    file_to_process = os.path.join(path, file)
                    if "SLOPE-" in file.upper():
                        with fits.open(file_to_process) as hdul:
                            header = hdul[0].header
                            Date_obs = header["DATE_OBS"] + "T" + header["TIME_OBS"]
                            dt = datetime.fromisoformat(Date_obs)
                            rounded_ms = round(dt.microsecond / 1000)
                            dt_rounded = dt.replace(microsecond=0) + timedelta(milliseconds=rounded_ms)
                            Date_obs = dt_rounded.isoformat(timespec='milliseconds')
                            dt_object = datetime.fromisoformat(header["DATE_OBS"])

                            # 2. Print the year
                            year = str(dt_object.year)
                            object = header["OBJECT"].upper() 
                            file_instrument_pc_path = header["PATH"]
                            # print(cycle)
                            if cycle in file_instrument_pc_path:
                                file_rel_path = file_instrument_pc_path.split(cycle)[1][1:]
                            pi = extract_p_folder(os.path.join(file_instrument_pc_path,header['FNAME']))
                            # print(pi)
                            if pi == None:
                                # print(file,file_instrument_pc_path)
                                pi = "PXX"
                            if object == "LAMP":
                                # print(file)
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            elif "CONT" in file.upper():
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")                            
                            elif "AR" in file.upper():
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")
                            elif "NE" in file.upper()   :
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")
                            elif "TEST" in file.upper()   :
                                lamp_name = "TEST-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    # print(f_path)
                            elif "FLAT" in file.upper() :
                                flat_name = "F-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",flat_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            elif header ["EPADU"]==1.07:
                                obj_name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",obj_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            else:
                                obj_name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,"hxrgproc_reprocessed_slope_images",obj_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                    else:
                        with fits.open(file_to_process) as hdul:
                            header = hdul[0].header
                            Date_obs = header["DATE_OBS"] + "T" + header["TIME_OBS"]
                            dt = datetime.fromisoformat(Date_obs)
                            rounded_ms = round(dt.microsecond / 1000)
                            dt_rounded = dt.replace(microsecond=0) + timedelta(milliseconds=rounded_ms)
                            Date_obs = dt_rounded.isoformat(timespec='milliseconds')
                            dt_object = datetime.fromisoformat(header["DATE_OBS"])

                            # 2. Print the year
                            year = str(dt_object.year)
                            object = header["OBJECT"].upper() 
                            file_instrument_pc_path = header["PATH"]
                            # print(cycle)
                            if cycle in file_instrument_pc_path:
                                file_rel_path = file_instrument_pc_path.split(cycle)[1][1:]
                            pi = extract_p_folder(os.path.join(file_instrument_pc_path,header['FNAME']))
                            # print(pi)
                            if pi == None:
                                # print(file,file_instrument_pc_path)
                                pi = "PXX"
                            if object == "LAMP":
                                # print(file)
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            elif "CONT" in file.upper():
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")                            
                            elif "AR" in file.upper():
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")
                            elif "NE" in file.upper()   :
                                lamp_name = "L-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    print("renamed")
                            elif "TEST" in file.upper()   :
                                lamp_name = "TEST-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,lamp_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                                    # print(f_path)
                            elif "FLAT" in file.upper() :
                                flat_name = "F-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,flat_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            elif header ["SIMPLE"]=="T":
                                obj_name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,obj_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)
                            else:
                                obj_name = "S-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                                f_path = os.path.join(final_path, file_rel_path,obj_name)
                                os.makedirs(os.path.dirname(f_path), exist_ok=True)
                                if not os.path.exists(f_path):
                                    os.rename(file_to_process, f_path)




                    # if not "SLOPE-" in file.upper():
                    #     print(file.upper(),object)
                           # Handle other object types
                           # = "O-"+cycle+pi+Date_obs+"-"+TELESCOPE+"-TANSPEC.fits"
                           # f_path = os.path.join(final_path, file_rel_path,other_name)
                           # os.makedirs(os.path.dirname(f_path), exist_ok=True)
                           # os.rename(file_to_process, f_path)









        #         if'.fit' in file :
        #             if file.endswith("Z.fits"):
        #                 print(path,file)

                        # shutil.copy(os.path.join(path, file), os.path.join(process_folder_path, file))
                        # with fits.open(os.path.join(process_folder_path, file)) as hdul:
                        #     # Process the FITS file
                        #     pass

        # process_dir_path =os.path.join(f"/data/{TELESCOPE}/TANSPEC/{cycle}/Processed_Data/Final_data",date)
        # for path,dirs,files in os.walk(process_folder_path):
        #     for file in files:
        #         # print(file)
        #         if'.fit' in file :
        #             print(file)
























main()