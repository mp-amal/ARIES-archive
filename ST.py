import os
import re
import sys
import glob
import time
import json
import ephem
# import aplpy
import shutil
import random
import smtplib
import calendar
import schedule
import subprocess
import numpy as np
import pandas as pd
import easygui as eg
from astropy import wcs
from astropy.io import fits
import os,sys,psycopg2,glob
from datetime import datetime,date
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from email.mime.multipart import MIMEMultipart
from astropy.coordinates import Angle, SkyCoord


from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import paramiko
from scp import SCPClient


from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import simple_norm




# -----------------------------------------------------------------------------
# NOTE:
# The paths below are DUMMY placeholders for sharing / public repos.
# Replace them with the real deployment paths on your archive server.
# -----------------------------------------------------------------------------
# Documentry folder (dummy paths)
database_finaldata = '/path/to/archived_data/astro_data/final_data'
thump_path         = '/path/to/archived_data/astro_data/astro_thumbnails/'
rel_path           = '/path/to/archived_data/astro_data/final_data'

DIR_DOC            = '/path/to/ADA_PROGRAM/config/'
final              = '/path/to/ST/2025B'

output_file        = DIR_DOC + "Process_folders.dat"
source_directory   = final + '/rawdata'
Base_folder        = final + '/Processed_data'
final_dest_folder  = final + '/Final_data/'

file_telescopes    = DIR_DOC + 'TelescopeList.csv'       # Contains telescope details
file_instruments   = DIR_DOC + 'InstrumentList.csv'      # Contains instrument details
file_keywords      = DIR_DOC + 'MasterHeaderList.dat'    # Contains header keyword details

today_date = date.today().strftime("%Y%m%d")

# JSON output path (dummy)
json_path = '/path/to/ADA_PROGRAM/output/JSON/daily_json.json'

remote_path         = './2025B'
processing_path     = '/cycle/rawdata'

# with open(json_path, "r") as json_file: 
#   jsondata=json.load(json_file)
# file_template = DIR_DOC+'/2021_temp/TemplateFile__'+{folder_name}+'.dat'
telescope_df = pd.read_csv(file_telescopes, sep=',', comment='#').set_index('ShortName') # DataFrame of telsecope details
instrument_df = pd.read_csv(file_instruments, sep=',', comment='#').set_index('ShortName') # DataFrame of instrument details
keywords_df = pd.read_csv(file_keywords, sep='\s+', comment='#').set_index('Header') # Details of MaseterHeder file
keywords_df = keywords_df.fillna('NULL')

DATE_keyword = 'DATE-OBS'
RA_keyword = 'RA'
DEC_keyword = 'DEC'
UT_keyword = 'UT'
FILTER_keyword = 'FILTER1'
OBJECT_keyword = 'OBJECT'


def copy_subfolders_from_server(remote_path, processing_path, server_ip, username, password):
    try:
        # Establish SSH connection
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(hostname=server_ip, username=username, password=password)

        # List all subdirectories inside remote_path
        # command = f"ls -d {remote_path}/"
        command = f"cd {remote_path} && ls"
        stdin, stdout, stderr = ssh.exec_command(command)
        folder_list = stdout.read().decode().strip().split('\n')
        folder_list = [folder.strip() for folder in folder_list if folder.strip()]
        # print(folder_list)
        # print(folder_list)
        # Create SCP client and copy folders
        with SCPClient(ssh.get_transport()) as scp:
            i = 0
            copied_folder =[]
            for folder in folder_list:
                i = i+1
                # if i == 6:
                #     break
                folder_name = os.path.basename(folder.rstrip('/'))
                if len(folder_name) == 8 and folder_name.isdigit():
                    dest_path = os.path.join(processing_path, folder_name)
                    if not os.path.exists(dest_path):
                        
                        print(f"Copying {folder_name} to {dest_path} ...")
                        print(f"🔄 {folder_name} ...")
                        scp.get(os.path.join(remote_path,folder_name), dest_path, recursive=True)
                        # print(f"Copied {folder_name} successfully.")
                        print(f"{folder_name} 📁✅.")
                        copied_folder.append(folder_name)
                    else:
                        print('folder is there')
        # return copied_folder
    except Exception as e:
        print(f"Error: {e}")
    finally:
        ssh.close()




def generate_template(subfolder_path, file_template):
    # Create the necessary directories for the output file
    os.makedirs(os.path.dirname(file_template), exist_ok=True)
    # print(output_file_path)
    # Open the output file for writing
    with open(file_template, 'w') as output_file:
        output_file.write("FileName\tObject\tFilter1\tpath\n")
        print("FileName\tObject\tFilter1\tpath\n")
        # Walk through the directory structure
        for path, subdirs, files in os.walk(subfolder_path):
            for file_nam in files:
                # print(file_name)
                path =os.path.basename(path)
                if "_" in file_nam[:4]:
                    # print(file_nam)
                    first_underscore_index = file_nam.find("_")
                    second_underscore_index = file_nam.find("_", first_underscore_index + 1)
                    col2 = file_nam[:second_underscore_index]
                    # print(col2)
                    col1 =file_nam
                    re=file_nam[second_underscore_index + 1:]
                    filter1= re.partition("_")[0]
                    output_file.write(f"{file_nam}\t{col2}\t{filter1}\t{path}\n")
                    print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                else:
                        col1= file_nam
                        if 'bais' in col1 or'bias' in col1 or 'Bais' in col1 or 'BIAS' in col1 or 'Bias' in col1 or 'BAIS'in col1:
                            col2 = 'BIAS'
                            filter1 = 'FREE'
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif 'flat' in col1:
                            col2 = 'FLAT'
                        if '_u' in col1 or '_U' in col1:
                            filter1 = 'U'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_b' in col1 or '_B' in col1:
                            filter1 = 'B'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_v' in col1 or '_V' in col1:
                            filter1 = 'V'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_r' in col1 or '_R' in col1:
                            filter1 = 'R'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_i' in col1 or '_I' in col1:
                            filter1 = 'I'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_Ha' in col1:
                            filter1 = 'Ha'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif '_OIII' in col1:
                            filter1 = 'OIII'
                            col2 = col1.partition('_')[0]
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        elif 'test' in col1 or 'TEST' in col1 or 'TST' in col1:
                        
                            filter1 = 'none'
                            col2 = "TEST"
                            output_file.write(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                            print(f"{col1}\t{col2}\t{filter1}\t{path}\n")
                        else:
                            name = os.path.join(path, file_nam.replace(' ','-'))
                            col1 = os.path.basename(name)
                            filter = 'FREE'
                            # col2 = col1.partition('_')[0]
                            col2 ="None"
                            # print(col2)
                            output_file.write(f"{col1}\t{col2}\t{filter}\n")
                            # print(f"{col1}\t{col2}\t{filter}\n")

def plot_fits_wcs(fits_file,temp_path,thumb_name):
    """
    Load a FITS file, extract WCS coordinates, and plot the image with WCS projection.
    
    Parameters:
    fits_file (str): Path to the FITS file.
    output_file (str): Path to save the output figure (default: 'fig.png').
    """
    
    # Load the FITS file
    hdul = fits.open(fits_file)
    header = hdul[0].header  # Use the primary header (or the correct extension)
    if not 'NAXIS3' in header.keys():
        image = hdul[0].data
        hdul.close()
        
        # Normalize the image for better visualization
        norm = simple_norm(image, 'log')
        
        # Create a WCS object
        wcs = WCS(header, naxis=2)
        
        # Get image dimensions
        xaxis= header['NAXIS1']
        yaxis= header['NAXIS2']
        # yaxis, xaxis = image.shape  
        
        # Plot the image with WCS projection
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': wcs})
        ax.imshow(image, origin='lower', cmap='gray', norm=norm)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        
        # Save the figure
        plt.savefig(temp_path+'/{}.png'.format(thumb_name),dpi=60)
        # plt.show()
        
        # print(f"Plot saved as {output_file}")
    else:
        print('KINETIC')

def rename__Calc_files(filename,subfolder_path):
    final_path =subfolder_path+"/Final_data"
    if not os.path.exists(final_path) :
        os.makedirs(final_path)
    data=fits.open(filename)
    header=data[0].header
    image=data[0].data
    obj = header['OBJECT']
    date_obs = header['DATE-OBS']
    telescope = header['TELESCOP']
    instru = header['INSTRUME']
    if (obj == 'test' or obj == 'TEST' or obj == 'Test'):
       print('test file')
    else:
        if (obj == 'flat' or obj == 'FLAT' or obj == 'Flat'):
            code = 'F'
        elif (obj == 'bias' or obj == 'BIAS' or obj == 'Bias'):
            code = 'B'

        else:
            code = 'S'
        new_name=code+'-2023APXX-'+date_obs+'-'+telescope+'-'+instru+'.fits'
        # print(new_name)
        try:
            os.rename(filename,final_path+'/'+new_name)
        except FileNotFoundError:
            pass
        print("\n"+"Calculation file renamed from ::: "+filename+" -----> "+new_name)
def move_fits_files(source_folder, destination_folder):
    # Ensure the destination folder exists
    # if not os.path.exists(destination_folder):
    #     os.makedirs(destination_folder)
    
    # Walk through the source folder, including subdirectories
    for root, dirs, files in os.walk(source_folder):

        # os.makedirs(destination_folder)
        for filename in files:
            if filename.endswith('.fits'):
                filepath = os.path.join(root, filename)
                # try:
                    # Open the FITS file and check NAXIS1 and NAXIS2
                # print(filename)
                with fits.open(filepath) as hdul:
                    header = hdul[0].header
                    naxis1 = header.get('NAXIS1')
                    naxis2 = header.get('NAXIS2')
                    if naxis1 == 512 and naxis2 == 512:
                        if not os.path.exists(destination_folder):
                        # If not, create the directory (and any intermediate directories if necessary)
                            os.makedirs(destination_folder)
                            print(f"polarimetry directory '{destination_folder}' created.")
                        else:
                            # If the directory already exists, just pass
                            print(f"Directory '{destination_folder}' already exists.")
                        print("axis :"+str(naxis1))
                        shutil.move(filepath, os.path.join(destination_folder, filename))
                        print(f"Moved {filename} to {destination_folder}")
                        # if not os.listdir(destination_folder):  # os.listdir() returns a list of files in the directory
                        #     os.rmdir(destination_folder)
                        #     print('Polarimetry folder named--------> '+os.path.basename(destination_folder) +' reomved because of it empty')

                        # send_email(subject, body, to_email,cc_emails)
                            # with open(json_path, "w") as json_file: 
                            #     jsondata[today_date]["ST"]["status"] = 'Not-completed'
                            #     jsondata[today_date]["ST"]["instrument"] = ' -'
                            #     jsondata[today_date]["ST"]["file count"] = '-'
                            #     jsondata[today_date]["ST"]["polarimetry"] = 'Yes'
                            #     jsondata[today_date]["ST"]["manual-check"] = 'check polarimetry'
                            #     json.dump(jsondata, json_file, indent=4)
def send_email(subject, body, to_email,cc_emails):
    # Email settings
    sender_email = "ariesdataarchive@gmail.com"
    sender_password = "password"
    smtp_server = "smtp.gmail.com"
    smtp_port = 000

    # Create a multipart email
    msg = MIMEMultipart()
    msg['From'] = sender_email
    msg['To'] = to_email
    msg['Subject'] = subject
    # msg['Cc'] = ', '.join(cc_emails)
    if cc_emails:
        # If multiple CCs are provided, join them into a single string
        msg['Cc'] = ', '.join(cc_emails)  # Add multiple CCs

    # Attach the email body to the message
    msg.attach(MIMEText(body, 'plain'))

    # Combine To and Cc for sending
    recipients = [to_email]  # To list
    if cc_emails:
        recipients += cc_emails  # Add all CCs to the recipients list

    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()  # Secure the connection
        server.login(sender_email, sender_password)
        text = msg.as_string()
        server.sendmail(sender_email, recipients, text)
        server.quit()

        print("Email sent successfully!")
    except Exception as e:
        print(f"Failed to send email: {str(e)}")
def extract_and_check_filter1(file_path,folder_name,file_template,f):
    # List of valid filters
    valid_filters = {'U', 'B', 'V', 'R', 'I', 'Ha', 'OIII', 'none', 'FREE', 'FLAT'}
    
    # Initialize an empty list to store FILTER1 values
    filter1_array = []

    # Open the .dat file
    with open(file_path, 'r') as file:
        # print(file_path)
        # Read the header to find the index of 'FILTER1'
        header = file.readline().strip().split()
        # print(header)
        filter1_index = header.index('Filter1')
        # print(filter1_index)
        # Loop through each line of the file after the header
        for line in file:
    #         # Split each line by whitespace
            row = line.strip().split()
    #         # Append the value from the 'FILTER1' column to the array
            filter1_array.append(row[filter1_index])

    # # Check if all values in filter1_array are in valid_filters
    if all(filter_value in valid_filters for filter_value in filter1_array):
        print("Correct file")
        f.write(f"{folder_name}\n")
        print("Okey, added to list")
    else:
        print("Not correct file")
        if os.path.exists(Base_folder+'/'+folder_name):
            shutil.rmtree(Base_folder+'/'+folder_name)
        if os.path.exists(file_template):
            os.remove(file_template)
        print("\n deleted folder and template file")
def rewrite_file_name(folder_path):
    # Iterate over all files in the specified folder
    for path, subdirs, files in os.walk(folder_path):
            for filename in files:
        # Create the old file path
                old_file_path = os.path.join(path, filename)
                # Check if it's a file (not a directory)
                if os.path.isfile(old_file_path):
                    # Replace spaces with hyphens in the filename
                    new_filename = filename.replace(" ", "-")
                    
                    # Create the new file path
                    new_file_path = os.path.join(path, new_filename)
                    
                    # Rename the file
                    os.rename(old_file_path, new_file_path)
                    # print(f'Renamed: {old_file_path} --------> {os.path.basename(new_file_path)}')
def get_month_abbreviation(month_value):
    # Dictionary mapping numerical month values to their abbreviations
    month_abbr = {
        "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
        "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
        "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec"
    }

    # Fetch the abbreviation from the dictionary
    return month_abbr.get(month_value, "Invalid month")
def find_final_folders(base_folder):
    final_folders = []
    for root, dirs, files in os.walk(base_folder):
        for dir_name in dirs:
            # print(dir_name)
            if dir_name == "Final_data":
                final_folders.append(os.path.join(root, dir_name))
    return final_folders
def copy_folder(source_dir, dest_dir, folder_name,output_file):
    """
    Copies a folder from source directory to destination directory with a timestamp.

    :param source_dir: The source directory containing the folder to copy.
    :param dest_dir: The destination directory to copy the folder to.
    :param folder_name: The name of the folder to copy.
    """
    try:
        src_folder_path = os.path.join(source_dir, folder_name)
        if not os.path.exists(src_folder_path):
            print(f"Folder '{folder_name}' does not exist in source directory.")
            return
        dest_folder_path = os.path.join(dest_dir, folder_name)
        if not os.path.exists(dest_folder_path):
            
            shutil.copytree(src_folder_path, dest_folder_path) 
            print(f"Folder '{folder_name}' copied to '{dest_dir}' successfully.")
        else:
            print(folder_name+ " Folder exits in working directry")
        
    except Exception as e:
        print(f"Error copying folder: {e}")
def sort_date(source_dir):
    folders = [f for f in os.listdir(source_dir) if os.path.isdir(os.path.join(source_dir, f))]
    # Regular expression pattern for date format YYYYMMDD
    date_pattern = re.compile(r"^\d{8}$")

    # List to store folders with valid date format
    valid_date_folders = []

    for folder in folders:
        if date_pattern.match(folder):
            try:
                # Convert to datetime object to ensure it's a valid date
                datetime.strptime(folder, "%Y%m%d")
                valid_date_folders.append(folder)
                # print(f"{folder} matches the date format YYYYMMDD")
            except ValueError:
                print(f"{folder} matches the pattern but is not a valid date")
        else:
            print(f"{folder} does not match the date format YYYYMMDD")

    # Sort the valid date folders in ascending order
    sorted_valid_date_folders = sorted(valid_date_folders)
    return sorted_valid_date_folders
def identify_instrument(instrument_df, filename):
    """
    Identify the Instrument details using the header of a given file 'filename' in directory to be run.
    Args:
        instrument_df : Pandas DataFrame containing Master-list of Instruments of a specific telescope
        filename      : FITS file from which Instrument details have to be extracted.
    Returns:
        instrument    : Name of the Instrument from which the given file 'filename' was observed.
    """
    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header

        if 'BAXIS1' in header.keys():
            header['XBINNING'] = header['BAXIS1']
            header['YBINNING'] = header['BAXIS2']
        elif 'HBIN' in header.keys():
            header['XBINNING'] = header['HBIN']
            header['YBINNING'] = header['VBIN']

        ysize = int(header['NAXIS1']) * int(header['XBINNING'])
        xsize = int(header['NAXIS2']) * int(header['YBINNING'])
        print(xsize)
        print(ysize)
        query = instrument_df[(instrument_df['NAXIS1'] == ysize) & (instrument_df['NAXIS2'] == xsize)]
        # print(query)
        instrument = query['Instrument'].values[0]

        return instrument
def format_dateobs(dateobs):
    """
    Formats the value of DATE_keyword in the FITS header to account for time in milliseconds.
    Args:
        dateobs : Value of the DATE_keyword that is to be modified
    Returns:
        datenew : Modified value of the DATE_keyword which accounts for time in milliseconds.
    """
    #datetime_master = datetime.strptime(dateobs, '%Y-%m-%dT%H:%M:%S.%f')
    #datenew = datetime_master.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
    ##neha
    # print("dateobs:", dateobs)
    try:
        datetime_master = datetime.strptime(dateobs, '%Y-%m-%d"T"%H:%M:%S.%f')
    except:
        try:
            datetime_master = datetime.strptime(dateobs, '%Y-%m-%dT%H:%M:%S.%f')
        except:
            datetime_master = datetime.strptime(dateobs, '%Y-%m-%dT%H:%M:%S')
    #if 12<datetime_master.hour<18: hour_new=datetime_master.hour-12  #time is wrongly written
    #datetime_master=datetime_master.replace(hour=hour_new)
    #print(datetime_master,hour_new,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
    datenew = datetime_master.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
    return datenew
def extract_extndata(hdulist, extn=0):
    """
    Extracts the data pertaining to a particular extension from a multi-dimensional array stored in a FITS file.
    Args:
        hdulist : HDU which contains the multi-dimensional FITS array
        extn    : Extension number for which the 2-D array is to be extracted
    Returns:
        datanew : 2-D array corresponding to the extension number (extn) specified
    """
    # print("working upto here")
    # print("::::::::::::::::::::::::::::::::::::::::: ------------------------>>  ",hdulist)
    data = hdulist[0].data
    # print("-----print form extract_extend--------------------->>>>>> ",data)
    datanew = np.reshape(data[extn, :, :], (data.shape[1], data.shape[2]))
    # print("-------------------------->>>>>> "+str(datanew))
    return datanew
def modify_dateobs(header, extn=0):
    """
    Modifies the value of DATE_keyword (i.e. time of observation) in the FITS header of multi-extension FITS
    files (Kinetic Mode Images).
    Args:
        header : FITS header which has to be modified for the time of observation
    Returns:
        header : Modified FITS header with the updated value of DATE_keyword
    """
    dateobs = header[DATE_keyword]
    datetime_master = datetime.strptime(dateobs, '%Y-%m-%dT%H:%M:%S.%f')
    datetime_new = datetime_master + extn * timedelta(seconds=header['ACT'])
    header[DATE_keyword] = datetime_new.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]

    return header
def identify_imagetype(filename):

    """
    Append Filter details onto the header of the file 'filename'.
    Args:
        filename : FITS file to which filter details have to be appended.
    Returns:
        True     : When the file is an object
        False    : When the file is a bias or a flat
    """
    
    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header

        if header[OBJECT_keyword] in ['BIAS', 'FLAT']:
            header['CATEGORY'] = 'Calibration'
            header['TYPE'] = header[OBJECT_keyword]
            return False
        else:
            header['CATEGORY'] = 'Science'
            header['TYPE'] = 'OBJ'
            # print(header['TYPE'])
            return True
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
def rename_files(filename,subfolder_path):
    
    final_path =subfolder_path+"/Final_data"
    temp_path =final_path + "/thumbnails"
    if not os.path.exists(final_path):
        os.makedirs(final_path)
    else:
        pass
    if not os.path.exists(temp_path) :
        os.makedirs(temp_path)
    else:
        pass
        
    # print(temp_path)
    # print(final_path)
    data=fits.open(filename)
    header=data[0].header
    image=data[0].data
    obj = header['OBJECT']
    date_obs = header['DATE-OBS']
    telescope = header['TELESCOP']
    instru = header['INSTRUME']
    if (obj == 'flat' or obj == 'FLAT' or obj == 'Flat'):
        code = 'F'
    elif (obj == 'bias' or obj == 'BIAS' or obj == 'Bias'):
        code = 'B'
    elif (obj == 'test' or obj == 'TEST' or obj == 'Test'):
        code = 'T' 
    else:
        code = 'S'
    new_name=code+'-2023APXX-'+date_obs+'-'+telescope+'-'+instru+'.fits'
    # print(new_name)
    try:
        os.rename(filename,final_path+'/'+new_name)
    except FileNotFoundError:
        pass
    #print(header.keys())

    (prefix, sep, suffix) = new_name.rpartition('.')
    thumb_name=prefix
    # print(thumb_path+'/{}.png'.format(thumb_name))

    if ('CRVAL1' in header.keys()):
        print(filename,'\n\n....................................\n\nastrometry done')
        plot_fits_wcs(final_path+'/'+new_name,temp_path,thumb_name)
def append_nullheader(filename):
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    with fits.open(filename, mode ='update') as hdulist:
        header = hdulist[0].header
        header.remove('OBJECT', ignore_missing=True, remove_all=True)
        for keyword in keywords_df.index:
            if keyword not in list(header.keys()):
                header.append(card=(keyword, keywords_df.loc[keyword, 'Value']))
        print(header)
        # headerdata = header
        hdulist.writeto(filename, overwrite=True)     
    print("<----- 1.  Null Header updated  ----->")
def append_templatefiledetails(filename,header_df,file_template):
    """
    Append Filter details onto the header of the file 'filename'.
    Args:
        filename : FITS file to which filter details have to be appended.
    Returns:
        None
    """
    with fits.open(filename, mode ='update') as hdulist :
            header = hdulist[0].header
            file =os.path.basename(filename.replace(' ','-'))
            columns = header_df.columns
            # print(columns)
            if file in header_df.index :
                # print(file)
                # print(filename)     
                for column in columns:
                    val = header_df.loc[header_df.index == file, column].values[0]
                    # print(val)
                    # print(column.upper()+" : "+ str(val).upper())            
                    if column.upper() in header.keys():
                            header[column.upper()] = str(val).upper()
                        #  print(column.upper()+" : "+ str(val).upper())
            else:
                print("ERROR: File '"+file+"' not logged in '~/config/"+os.path.basename(file_template)+"'")
                
    print("\n"+"<----- 2.  Tempalte Header updated  ----->")
def append_telescopedetails(TELESCOPE, instrument_df, filename):

    print('\n')
    print(filename)
    print('\n')
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        instrument_df : Pandas DataFrame containing Master-list of Instruments
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    _, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, HORIZON, ZENITH = telescope_df.loc[TELESCOPE].values
    instrument_df = instrument_df.loc[TELESCOPE]

    instrument = identify_instrument(instrument_df, filename)
    if instrument in instrument_df['Instrument'].values:
        INSTRUME, PIXSCALE, RNOISE, GAIN, _, _ = instrument_df.loc[instrument_df['Instrument'] == instrument].values[0]

        dict_append = {'TELESCOP': TELESCOPE, 'GEOLONG': OBS_LONG, 'GEOLAT': OBS_LAT, 'GEOELEV': OBS_ALT,
                       'HORIZON': HORIZON, 'ZENITH': ZENITH, 'TIMEZONE': OBS_TIMEZONE, 'INSTRUME': INSTRUME,
                       'PIXSCALE': PIXSCALE, 'DET-READ': RNOISE, 'DET-GAIN': GAIN}

        with fits.open(filename, mode='update') as hdulist:
            header = hdulist[0].header
            for keyword, value in dict_append.items():
                header[keyword] = value
    else:
        pass
    print("\n"+"<----- 3.  Telescope Header updated  ----->")
def format_header(TELESCOPE, filename):
    """
    Formats and homogenizes the header of the file 'filename'. Kinetic mode files are sliced into different FITS
    files and DATE_keyword is updated accordingly.
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    instrument = identify_instrument(instrument_df, filename)
    hdulist = fits.open(filename, mode='update')
    header = hdulist[0].header
    header['ORIGFILE'] = os.path.basename(filename)

    # print('step111111111111111',telescopename)
    if TELESCOPE == 'ST':
        # #if instrument == '1k x 1k Imager':   #Neha
        # if instrument == '1K_IMG2POL':        #Neha
        #     date = header['DATE_OBS'].replace('/', '-')
        #     print(date,'++++++++++++++++++++++++')
        #     time_utc = header['TIME']
        #     dateobs = format_dateobs(date + 'T' + time_utc)
        #     print(dateobs,'$$$$$$$$$$$$$$$$$$$$$')
        #     header[DATE_keyword] = dateobs 
        #     #header.extend(cards=(DATE_keyword, dateobs), update=True)
        #     #header.remove('OBJECT', ignore_missing=True, remove_all=True)
        #     print('test++++++++++++++++++++++++')

        #if instrument == '1k x 1k Imager':   #Neha
        if instrument == '1K_IMG2POL':        #Neha
            date = header['DATE_OBS'].replace('/', '-')
            print(date,'++++++++++++++++++++++++')
            time_utc = header['TIME']
            dateobs = format_dateobs(date + 'T' + time_utc)
            print(dateobs,'$$$$$$$$$$$$$$$$$$$$$')
            header[DATE_keyword] = dateobs 
            #header.extend(cards=(DATE_keyword, dateobs), update=True)
            #header.remove('OBJECT', ignore_missing=True, remove_all=True)
            print('test++++++++++++++++++++++++')
    
        if instrument == '1K_IMG1':     
            date = header['DATE_OBS'].replace("/", '-')
            print(date,'++++++++++++++++++++++++')
            time_utc = header['TIME']
            dateobs = format_dateobs(date + 'T' + time_utc)
            print(dateobs,'$$$$$$$$$$$$$$$$$$$$$')
            header[DATE_keyword] = dateobs
            #header.extend(cards=(DATE_keyword, dateobs), update=True)
            #header.remove('OBJECT', ignore_missing=True, remove_all=True)
            print('test++++++++++++++++++++++++')
    
        if instrument == '2K_IMG1':
            #date = header['DATE-OBS'].replace("/", '-')
            #date = header['UTDATE'].replace("/", '-')
            date = header['UTDATE']
            print(date,'++++++++++++++++++++++++')
            time_utc = header['UTSTART']
            h, m, s = time_utc.split(":")  # second should be less than 60!
            if s == '60':
                m = str(int(m)+1)
                s = '00'
            time_utc = h+':'+m+':'+s
            dateobs = format_dateobs(date + 'T' + time_utc)
            print(dateobs,'$$$$$$$$$$$$$$$$$$$$$')
            header[DATE_keyword] = dateobs
            #header.extend(cards=(DATE_keyword, dateobs), update=True)
            #header.remove('OBJECT', ignore_missing=True, remove_all=True)
            print('test++++++++++++++++++++++++')
    
        if instrument == '1K_IMG2':
            dateobs = header['DATE-OBS']
            dateobs = format_dateobs(dateobs)
            header[DATE_keyword] = dateobs
            print('dataobs :  ',dateobs)

        if instrument == '4K_IMG2':
            datetime = header['DATE-OBS']
            datetime2= str.split(datetime, "T")
            date = datetime2[0]
            time_utc=datetime2[1]
            dateobs = date + 'T' + time_utc
            dateobs = format_dateobs(dateobs)
            header[DATE_keyword] = dateobs
            print('^^^^^^^^^dateobs:^^^^^^^^^^^^^',dateobs)
    elif TELESCOPE == 'DFOT':
        if 'FRAME' in header.keys():
            header[DATE_keyword] = format_dateobs(header['FRAME'])
        elif 'DATE' in header.keys():
            header[DATE_keyword] = format_dateobs(header['DATE'])
        else:
            print("ERROR: DATE keyword not found in the Header")
            pass
            # sys.exit(1)

    if 'EXPOSURE' in header.keys():
        header[EXPOSURE_keyword] = header['EXPOSURE']
        header.remove('EXPOSURE', remove_all=True)

    if int(header['NAXIS']) == 3:
        if int(header['NAXIS3']) == 1:
            hdulist[0].data = extract_extndata(hdulist)
        else:
            for extn in range(0, int(header['NAXIS3'])):
                datanew = extract_extndata(hdulist, extn=extn)
                headernew = modify_dateobs(header.copy(), extn=extn)
                hdunew = fits.PrimaryHDU(data=datanew, header=headernew)
                hdunew.writeto("{0}_{1}.fits".format(filename.split('.')[0], extn + 1), overwrite=True)

    hdulist.flush()
    hdulist.close()
    print("\n"+"<----- 4.  Header formating done  ----->")
def do_astrometry(telescopename, filename, pixtolerance=0.02):
    """
    Performs Astrometry on images with constraints specified by RA and DEC if UT is specified.
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file on which Astrometry has to be performed
        pixtolerance  : Platescale tolerance in percentage
    Returns:
        None
    """
    print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    _, _, OBS_LAT, _, _, _, _ = telescope_df.loc[telescopename].values

    header = fits.getheader(filename, ext=0)
    date_obs = header[DATE_keyword]
    #date_obs = '2020-05-18T16:33:57.373'
    print("date_obs:", date_obs)
    #print(header[OBJECT_keyword],'**************')
    #print(header['CATEGORY'],'*****************')
    #if header[OBJECT_keyword] in ['BIAS', 'FLAT']:
        #header['CATEGORY'] = 'Calibration':
    #    print('calibration file ******')
    #    pass

    #else:
    print('doing astrometry ******')
    try:
        PIXSCALE = header['PIXSCALE']
        PIX_l = str(PIXSCALE - (PIXSCALE * pixtolerance))
        PIX_u = str(PIXSCALE + (PIXSCALE * pixtolerance))
        PIX_l = str(0.15)
        PIX_u = str(2.0)
        print('pixel_lower=',PIX_l,'pixel_upper=',PIX_u,'+++++%%%%%')
        print(date_obs, 'date_obs') 
        try: LST_deg = get_LST(telescopename, date_obs)
        except: 
            datetime2= str.split(date_obs, "T")
            date = datetime2[0]
            time_utc=datetime2[1]
            dateobs = date + 'T' + time_utc
            dateobs = format_dateobs(dateobs)
            with fits.open(filename, mode='update') as hdul:
                hdr = hdul[0].header
                hdr[DATE_keyword] = dateobs
            LST_deg = get_LST(telescopename, dateobs)
        print(LST_deg,date_obs,'LST_deg & date_obs')
        #-----------
        #LST_deg = 185.0167
        #OBS_LAT = 5.3428
        #------------
        # if telescopename == 'DFOT':
        try:
            #subprocess.call('solve-field --continue --downsample 2 --no-plots --scale-units app --config '+DIR_DOC+'Astrometry.cfg --ra ' + str(LST_deg) + ' --dec ' + str(OBS_LAT) + ' --radius 50 ' + filename, timeout=2000, shell=True)
            subprocess.call('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIX_l + ' --scale-high ' + PIX_u  + ' --scale-units app --config '+DIR_DOC+'Astrometry.cfg --ra ' + str(LST_deg) + ' --dec ' + str(OBS_LAT) + ' --radius 20 ' + filename, timeout=300, shell=True)
        except subprocess.TimeoutExpired:
            print("Astrometry Timed Out For "+str(filename))
            return False
        else:
            print("Astrometry Ran Sucessfully For "+str(filename))
            os.system('rm -rf *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs')
            return True
        
        print('Astrometry successfull---------------------------111--') 
        # elif telescopename == 'ST':
        #    os.system('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIXSCALE_l + ' --scale-high '
        #              + PIXSCALE_u + ' --scale-units app --config Astrometry.cfg  ' + filename)
        # else:
        #    sys.exit(1)
    except KeyError:
        print("ERROR: Header Keyword 'PIXSCALE' not found")
        return False
    print('Astrometry successfull-----------------------------')
def calculate_airmass(TELESCOPE, filename):
    """
    Calculates AIRMASS for the FITS file and appends respective details in the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file whose header has to be edited
    Returns:
        None
    """
    hdulist = fits.open(filename, mode='update')
    header = hdulist[0].header

    date_obs = header[DATE_keyword]
    object_ra = header[RA_keyword]
    object_dec = header[DEC_keyword]
    
    date_obs, time_utc = date_obs.split('T')
    datetime_utc = str(date_obs) + ' ' + str(time_utc)
    julian_day = ephem.julian_date(str(datetime_utc))

    telescope = init_telescope(TELESCOPE)
    telescope.date = datetime_utc
    time_sidereal = telescope.sidereal_time()

    object_obs = ephem.FixedBody()
    object_obs._ra = object_ra
    object_obs._dec = object_dec
    object_obs._epoch = ephem.J2000
    object_obs.compute(telescope)

    lon= header['GEOLONG']
    lat= header['GEOLAT']
    c = SkyCoord(lon, lat, frame = 'icrs', unit = 'deg')
    lon = c.ra.value
    lat = c.dec.value
    lst = str(time_sidereal)
    LST = Angle(lst + ' hours').deg
    # Alt & Az calculated based on equations of http://www.stargazing.net/kepler/altaz.html
    # Also see the calculator: http://jukaukor.mbnet.fi/star_altitude.html
    ha = LST - float(object_ra)
    if ha < 0.0: ha = ha + 360
    sin_dec = np.sin(np.radians(float(object_dec)))
    sin_lat = np.sin(np.radians(lat))
    cos_lat = np.cos(np.radians(lat))
    cos_dec = np.cos(np.radians(float(object_dec)))
    # Altitude:
    alt0 = sin_dec * sin_lat + cos_dec * cos_lat * np.cos(np.radians(ha))
    alt1 = np.arcsin(alt0)     #in radians
    alt  = np.degrees(alt1)    #in degrees
    sin_alt = alt0
    cos_alt = np.cos(np.radians(alt))
    #Azimuth:
    az0 = np.arccos( (sin_dec - sin_alt*sin_lat) / (cos_alt*cos_lat ))
    if np.sin(np.radians(ha))<0: az = np.degrees(az0)
    else: az = 360 - np.degrees(az0)
    airmass = 1 / np.cos(np.radians(90 - alt))
    alt_angle = Angle(alt, u.deg).to_string(unit = u.degree, sep=':')
    az_angle = Angle(az, u.deg).to_string(unit = u.degree, sep=':')
    #dict_header = {'JD': str(julian_day), 'LST': str(time_sidereal), 'ALTITUDE': str(object_obs.alt),
                  # 'AZIMUTH': str(object_obs.az), 'AIRMASS': str(airmass)}
    dict_header = {'JD': str(julian_day), 'LST': str(time_sidereal), 'ALTITUDE': alt_angle,
                   'AZIMUTH': az_angle, 'AIRMASS': str(airmass)}
    for keyword, value in dict_header.items():
        if keyword in header.keys():
            header.remove(keyword, ignore_missing=True, remove_all=True)
        header.append(card=(keyword, value))

    hdulist.flush()
    hdulist.close()
    print("\n"+"<----- .  Airmass calculation updated  ----->")
def append_astrometryheader(filename,new_filename):
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
    hdulist = []
    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header

        for keyword in astrometrykeys:
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header.append(card=(keyword, headernew[keyword]))

        for idx, keyword in enumerate([RA_keyword, DEC_keyword]):
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header. append(card=(keyword, str(radec[idx])))
    print("\n"+"<----- .  Astrometry header updated  ----->")
def check_astrometry(filename,subfolder_path,TELESCOPE):

    """ 
    Checks whether astrometry is done on the files or not. If astrometry failed on a
    source, it moves the file to the failed_astrometry/ directory.
    """
    if ('fits' in filename):
            data=fits.open(filename)
            header=data[0].header
            try:
                ra,dec=header['RA'],header['DEC']
                filters=header['FILTER1']
                category=header['CATEGORY']
                if (category=='Science'):
                    if (ra=='NULL'):
                        failed_path =subfolder_path + "/failed_astrometry"
                        failed_log_path = subfolder_path + "/failed_astrometry/failed_log.txt"
                        if not os.path.exists(failed_path):
                            os.makedirs(failed_path)
                        # print(temp_path)
                        print("\n"+filename,ra,dec,filters,'\t','astrometry failed')
                        with open(failed_log_path, 'a') as log_file:
                            log_file.write(os.path.basename(filename) + '\n')
                        shutil.move(filename,failed_path)
                    else:
                        print("\n"+filename,ra,dec,filters, '\t','astrometry successful')
                        rename_files(filename,subfolder_path)  
            except KeyError:
                print("\n"+filename,filters,'\t','astrometry failed')
                pass
    print("\n"+"<----- .  Astrometry failed file seperated  ----->")
def extract_and_check_filter1(file_path,folder_name,file_template,f):
    # List of valid filters
    valid_filters = {'U', 'B', 'V', 'R', 'I', 'Ha', 'OIII', 'none', 'FREE', 'FLAT'}
    
    # Initialize an empty list to store FILTER1 values
    filter1_array = []

    # Open the .dat file
    with open(file_path, 'r') as file:
        # print(file_path)
        # Read the header to find the index of 'FILTER1'
        header = file.readline().strip().split()
        # print(header)
        filter1_index = header.index('Filter1')
        print(filter1_index)
        # Loop through each line of the file after the header
        for line in file:
    #         # Split each line by whitespace
            row = line.strip().split()
    #         # Append the value from the 'FILTER1' column to the array
            filter1_array.append(row[filter1_index])
    print(filter1_array)
    # # Check if all values in filter1_array are in valid_filters
    if all(filter_value in valid_filters for filter_value in filter1_array):
        print("Correct file")
        f.write(f"{folder_name}\n")
        print("Okey, added to list")
    else:
        print("Not correct file")
        if os.path.exists(Base_folder+'/'+folder_name):
            shutil.rmtree(Base_folder+'/'+folder_name)
        if os.path.exists(file_template):
            os.remove(file_template)
        print("\n deleted folder and template file")
def astrometry_header(filename,subfolder_path,TELESCOPE):
    if "flat" in filename or "Flat" in filename or "FLAT" in filename or "bais" in filename or "Bais" in filename or "Bias" in filename or"BIAS" in filename or "BAIS" in filename or "bias" in filename or 'test' in filename or 'Test' in filename or 'TEST' in filename:
        rename__Calc_files(filename,subfolder_path)
    else:
        if identify_imagetype(filename):
            do_astrometry(TELESCOPE, filename)
            new_filename = filename.split('.')[0] + '.new'
            if os.path.exists(new_filename):
                append_astrometryheader(filename,new_filename)
                calculate_airmass(TELESCOPE,filename)
                print("updated astrometry header of :"+os.path.basename(filename))
            else:
                print("Astrometry faild in '"+os.path.basename(filename)+"'")
            check_astrometry(filename,subfolder_path,TELESCOPE)  

def datedisplay(source_directory, Base_folder, output_file,ydatefolder):
    folder_list = sort_date(source_directory)
    process_list = sort_date(Base_folder)
    polarimetry_list = sort_date(final +'/Polarimetry/')
    
    print("\n" * 3)
    folder_in_action =[]
    for folder in folder_list:
        if folder not in process_list:
            if folder not in polarimetry_list:
                print(folder)
                folder_in_action.append(folder)
            # process_fits_files(source_directory, folder, instrument_df)
    # Get user input for the folders to copy
    print(folder_in_action)
    if not  ydatefolder in folder_list:
        # send_email(subject, body, to_email,cc_emails)
        print('NO folder' )
        # with open(json_path, "w") as json_file: 
        #     jsondata[today_date]["ST"]["status"] = 'Not-completed'
        #     jsondata[today_date]["ST"]["instrument"] = ' -'
        #     jsondata[today_date]["ST"]["file count"] = '-'
        #     jsondata[today_date]["ST"]["polarimetry"] = '-'
        #     jsondata[today_date]["ST"]["manual-check"] = 'No Folder'
        #     json.dump(jsondata, json_file, indent=4)
    else:
        print('\nFolder is there...')
    folders_to_copy= folder_in_action
    
    with open(output_file, 'w') as f:
        f.write("Folder name\n")
        for folder_name in folders_to_copy:
            if folder_name in folder_list:
                # f.write(f"{folder_name}\n")
                copy_folder(source_directory, Base_folder, folder_name, output_file)
            else:
                print(f"Folder '{folder_name}' is not available in the source directory.")
            rewrite_file_name(Base_folder+'/'+folder_name)            
            desitnation_folder = final +'/Polarimetry/'+folder_name
            move_fits_files(Base_folder+'/'+folder_name, desitnation_folder)
            file_template = DIR_DOC+'2025B_temp/TemplateFile_ST_'+str(folder_name)+'.dat'

            if os.listdir(Base_folder+'/'+folder_name):
                generate_template(Base_folder+'/'+folder_name, file_template)
                extract_and_check_filter1(file_template,folder_name,file_template,f)
            else:
                print('Folder is empty ')
                shutil.rmtree(Base_folder+'/'+folder_name)
                print('Folder is empty and removed')
def header_update(filename,header_df,file_template,TELESCOPE,instrument_df):
    append_nullheader(filename)
    append_templatefiledetails(filename,header_df,file_template)
    append_telescopedetails(TELESCOPE, instrument_df, filename)
    format_header(TELESCOPE, filename)


def plot_fits_wcs(fits_file,temp_path,thumb_name):
    """
    Load a FITS file, extract WCS coordinates, and plot the image with WCS projection.
    
    Parameters:
    fits_file (str): Path to the FITS file.
    output_file (str): Path to save the output figure (default: 'fig.png').
    """
    
    # Load the FITS file
    hdul = fits.open(fits_file)
    header = hdul[0].header  # Use the primary header (or the correct extension)
    if not 'NAXIS3' in header.keys():
        image = hdul[0].data
        hdul.close()
        
        # Normalize the image for better visualization
        norm = simple_norm(image, 'log')
        
        # Create a WCS object
        wcs = WCS(header, naxis=2)
        
        # Get image dimensions
        xaxis= header['NAXIS1']
        yaxis= header['NAXIS2']
        # yaxis, xaxis = image.shape  
        
        # Plot the image with WCS projection
        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': wcs})
        ax.imshow(image, origin='lower', cmap='gray', norm=norm)
        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        
        # Save the figure
        plt.savefig(temp_path+'/{}.png'.format(thumb_name),dpi=60)
        # plt.show()
        
        # print(f"Plot saved as {output_file}")
    else:
        print('KINETIC')

server_ip = "000.00.0.00"  # Replace with your server's IP
username = "servername"    # Replace with your username
password = "password"    # Replace with your password

def main():
    # print('mainfunction')
    TELESCOPE= 'ST'
    today = datetime.today()

    copy_subfolders_from_server(remote_path, processing_path, server_ip, username, password)


    # Calculate yesterday's date
    yesterday = today - timedelta(days=1)
    ydatefolder = str(yesterday.date()).replace('-','')
    datedisplay(source_directory, Base_folder, output_file,ydatefolder)
    process_folders =DIR_DOC+'Process_folders.dat'
    process_file_df = pd.read_csv(process_folders, comment='#').set_index('Folder name')
    for subfolder in process_file_df.index: 
        subfolder_path = os.path.join(Base_folder,str(subfolder))
        file_template = DIR_DOC+'/'+cycle+'_temp/TemplateFile_ST_'+str(subfolder)+'.dat'
        header_df = pd.read_csv(file_template, sep='\s+', comment='#', dtype='string').set_index('FileName')
        # print(header_df)
        for path, subdirs, files in os.walk(subfolder_path):
            for file in files:
                if not 'test' in file and not 'TEST' in file:
                    if not file.startswith('.'):
                        if 'fits' in file :
                            if not 'test' in file or not 'TEST' in file or not 'Test' in file:
                                filename = path+'/'+file 
                                # print(filename)
                                print("\n"+"file name is :--------------   : "+filename+"\n"+"\n")
                                header_update(filename,header_df,file_template,TELESCOPE,instrument_df)
                            print("\n"+":::::::::: Header upadated for '"+os.path.basename(file)+ "' ::::::::::")
            print("\n"+"\n"+"--------------   : header updation completed for : "+ str(subfolder))
        for path, subdirs, files in os.walk(subfolder_path):
                for file in files:
                    if not 'test' in file and not 'TEST' in file:
                        if not file.startswith('.'):
                            if 'fits' in file:
                                filename = path+'/'+file       
                                if "flat" in file or "Flat" in file or "FLAT" in file or "bais" in file or "Bais" in file or "Bias" in file or"BIAS" in file or "BAIS" in file or "bias" in file :
                                    # print(os.path.basename(filename))
                                    rename__Calc_files(filename,subfolder_path)
                                else:
                                    if identify_imagetype(filename):
                                        do_astrometry(TELESCOPE, filename)
                                        new_filename = filename.split('.')[0] + '.new'
                                        if os.path.exists(new_filename):
                                            append_astrometryheader(filename,new_filename)
                                            calculate_airmass(TELESCOPE,filename)
                                            print("\nupdated astrometry header of :"+os.path.basename(filename))
                                        else:
                                            print("\nAstrometry faild in '"+os.path.basename(filename)+"'")
                                        check_astrometry(filename,subfolder_path,TELESCOPE)
        succes_log_path = subfolder_path+"/Final_data/succes_log.txt"
        for filename in glob.glob(os.path.join(subfolder_path+"/Final_data",'*.fits')):
            with open(succes_log_path, 'a') as log_file:
                log_file.write(os.path.basename(filename) + '\n')
        print('\n final data log added')

    final_folders = find_final_folders(Base_folder)
    print("List of 'Final_folder' directories:")
    for subfolder in process_file_df.index:
        subfolder_path = os.path.join(Base_folder,str(subfolder))+"/Final_data"
        # print(subfolder_path)
        for folder in final_folders:
            if folder == subfolder_path :
                    # print("empty")
    #                 # print(folder)
                fits_files = [f for f in os.listdir(folder) if f.endswith('.fits') and os.path.isfile(os.path.join(folder, f))]
                print(len(fits_files))
                if len(fits_files) !=0 :
                    random_file = random.choice(fits_files)
                    # print(random_file)
                    random_split = random_file.split('-')
                    
                    year = random_split[2]
                    mon =random_split[3]
                    month =get_month_abbreviation(mon)
                    day =random_split[4].split('T')[0]
                    date_folder=os.path.basename(str(subfolder))
                    telescope =random_split[5] 
                    instrument =random_split[6].split('.')[0] 
                    db_dest_path = database_finaldata+'/'+str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)
                    dest_path = final_dest_folder+str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)
                    thumpnail_path = thump_path +str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)+'/thumbnails'
                    if os.path.exists(dest_path):
                        print('Destination exist with path : '+dest_path)
                    else:
                        shutil.copytree(folder,db_dest_path)
                        print(folder+ "--------------> copied to db")
                        # print(folder+ "--------------> copied to final data")
                        # shutil.copytree(folder,dest_path)
                        if os.path.exists(thumpnail_path):
                            shutil.copytree(folder+'/thumbnails',thumpnail_path)

                        
                    
    # ---------------Fits to db prgoram added here.................................................

                    # file_path='/data/archived_data/astro_data/final_data/2024/Oct/DFOT/2K_IMG1/20241009'
                    file_path = (os.path.relpath(db_dest_path, rel_path)) # Relative file path to be inserted in the database

                    nofp=0 # No of files processed
                    dataarray=[]
                    if os.path.exists(file_path + "/"+"log.txt"):
                        print ("Directory seems to be processed. If you want to reprocess the directory, then delete log.txt file and then run this program ")
                        exit()
                    infile=open(str(db_dest_path) + "/"+"log.txt","w")

                    #Total Number of Fits file in the directory
                    totfits = len(glob.glob1(db_dest_path+"/","*.fits"))
                    #Master Log file located just outside of DataArchive folder structure
                    masterlog = open(os.path.join("/data/archived_data/astro_data/final_data/", "masterlog.txt"), 'a')

                    #infile.write('  {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t \n'.format("OBSTimestamp","Observer","Object","RA","DEC","Filter","Filter1","Filter2","Instrument","FileName","Full Path","PROPID"))

                    #Database connectivity
                    try:
                        connection = psycopg2.connect(user="username",
                                                    password="password",
                                                    host="db host ip",
                                                    port="5432",
                                                    database="darchive")
                        cursor = connection.cursor()
                    except Exception as error:
                            infile.write(str(error))
                            print (error)
                            


                    for fname in os.listdir(db_dest_path):
                        if os.path.splitext(fname)[1]!=".fits":continue
                        fpath=os.path.join(db_dest_path,fname)
                        #if os.path.isdir(fpath):continue
                        hdulist=fits.open("%s"%fpath)
                        
                        #hdulist.verify('silentfix')
                        try:
                            
                            OBSDate,OBSTime = str(hdulist[0].header['DATE-OBS']).split("T")
                            DATEOBS = OBSDate +" " + OBSTime
                            
                            if "Undefined" in str(hdulist[0].header['OBJECT']):
                                hdulist[0].header["OBJECT"]="NULL"
                            OBJECT = hdulist[0].header["OBJECT"]
                            
                            if ("Undefined" in str(hdulist[0].header['RA'])) or (hdulist[0].header["RA"]=="NULL"):
                                hdulist[0].header["RA"]= 0.000
                            RA = hdulist[0].header["RA"]
                            
                            if "Undefined" in str(hdulist[0].header['DEC']) or (hdulist[0].header["DEC"]== "NULL"):
                                hdulist[0].header["DEC"]=0.000
                            DEC = hdulist[0].header["DEC"]
                                
                            if "Undefined" in str(hdulist[0].header['PROPNO']):
                                hdulist[0].header["PROPNO"]="NULL"	
                            PROPNO = 	hdulist[0].header["PROPNO"]
                            
                            if "Undefined" in str(hdulist[0].header['PROG']):
                                hdulist[0].header["PROG"]="NULL"	
                            PROG = 	hdulist[0].header["PROG"]
                            
                            if "Undefined" in str(hdulist[0].header['OBSERVER']):
                                hdulist[0].header["OBSERVER"]="NULL"	
                            OBSERVER = 	hdulist[0].header["OBSERVER"]
                            
                            if "Undefined" in str(hdulist[0].header['TELESCOP']):
                                hdulist[0].header["TELESCOP"]="NULL"	
                            TELESCOP = 	hdulist[0].header["TELESCOP"]
                            if "Undefined" in str(hdulist[0].header['INSTRUME']):
                                hdulist[0].header["INSTRUME"]="NULL"	
                            INSTRUMENT = 	hdulist[0].header["INSTRUME"]
                            
                            if "Undefined" in str(hdulist[0].header['CATEGORY']):
                                hdulist[0].header["CATEGORY"]="NULL"	
                            CATEGORY = 	hdulist[0].header["CATEGORY"]
                            
                            if "Undefined" in str(hdulist[0].header['TYPE']):
                                hdulist[0].header["TYPE"]="NULL"	
                            TYPE = 	hdulist[0].header["TYPE"]
                            
                            if "Undefined" in str(hdulist[0].header['MODE']):
                                hdulist[0].header["MODE"]="NULL"	
                            MODE = 	hdulist[0].header["MODE"]
                            
                            DATAID = fname
                            
                            if "Undefined" in str(hdulist[0].header['ORIGFILE']):
                                hdulist[0].header["ORIGFILE"]="NULL"	
                            ORIGFILE = 	hdulist[0].header["ORIGFILE"]
                            
                            
                            RELEASE = "NULL"
                            
                            
                            if "Undefined" in str(hdulist[0].header['EXPTIME']):
                                hdulist[0].header["EXPTIME"]="NULL"	
                            EXPOSURE = 	hdulist[0].header["EXPTIME"]
                            
                            if "Undefined" in str(hdulist[0].header['FILTER1']):
                                hdulist[0].header["FILTER1"]="NULL"	
                            FILTER1 = 	hdulist[0].header["FILTER1"]
                            
                            if "Undefined" in str(hdulist[0].header['FILTER2']):
                                hdulist[0].header["FILTER2"]="NULL"	
                            FILTER2 = 	hdulist[0].header["FILTER2"]
                            
                            if "Undefined" in str(hdulist[0].header['GRISM']):
                                hdulist[0].header["GRISM"]="NULL"
                            GRISM = hdulist[0].header["GRISM"]
                            
                            if "Undefined" in str(hdulist[0].header['SLIT']):
                                hdulist[0].header["SLIT"]="NULL"	
                            SLIT = 	hdulist[0].header["SLIT"]
                            
                            if "Undefined" in str(hdulist[0].header['AIRMASS']):
                                hdulist[0].header["AIRMASS"]="NULL"	
                            AIRMASS = 	round(float(hdulist[0].header["AIRMASS"]),3)
                            
                            FILESIZE = os.path.getsize(fpath)
                            
                            FILEPATH = file_path
                            
                            TIMEINSERTED = datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S.%s")				
                            
                                    
                            datarow=(DATEOBS,OBJECT,RA,DEC,PROPNO,PROG,OBSERVER,TELESCOP,INSTRUMENT,CATEGORY,TYPE,MODE,DATAID,ORIGFILE,EXPOSURE,FILTER1,FILTER2,GRISM,SLIT,AIRMASS,FILESIZE,FILEPATH,TIMEINSERTED)
                            
                                    
                            sql_insert='''INSERT INTO astronomy(DATEOBS,OBJECT,RA,DEC,PROPNO,PROG,OBSERVER,TELESCOP,INSTRUMENT,CATEGORY,TYPE,MODE,DATAID,ORIGFILE,EXPOSURE,FILTER1,FILTER2,GRISM,SLIT,AIRMASS,FILESIZE,FILEPATH,TIMEINSERTED) VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);'''

                            cursor.execute(sql_insert, (datarow))
                            
                            dataarray.append(datarow)
                            infile.write(str(DATAID) + "\n")
                        except Exception as error:
                            infile.write(str(error))
                            print (error)
                            connection.rollback()
                        else:
                            connection.commit()

                        hdulist.close()
                            
                        nofp=nofp+1
                    #print (dataarray)
                    infile.write("Total No of fits files  = %d \n" %totfits)
                    infile.write("Total No of files processed = %d" %nofp)  
                    print ("Total number of fits files present in the directory = %d" %totfits)
                    print ("Total No of files processed = %d" %nofp)
                    masterlog.write('{}\t {}\t {}\t {}\t \n'.format(TIMEINSERTED,db_dest_path,nofp,totfits))
                    if connection:
                        cursor.close()
                        connection.close


                        
                    infile.close()
                    subfolder= str(subfolder)
                    # send_email(subject, body, to_email,cc_emails)
                    # with open(json_path, "w") as json_file: 
                    #     jsondata[subfolder]["ST"]["status"] = 'Completed'
                    #     jsondata[subfolder]["ST"]["instrument"] = INSTRUMENT
                    #     jsondata[subfolder]["ST"]["file count"] = nofp
                    #     jsondata[subfolder]["ST"]["polarimetry"] = 'No'
                    #     jsondata[subfolder]["ST"]["manual-check"] = 'No'
                    #     json.dump(jsondata, json_file, indent=4)
                else:
                    print('Folder is empty')

                    # send_email(subject, body, to_email,cc_emails)
                    # with open(json_path, "w") as json_file: 
                    #     jsondata[subfolder]["ST"]["status"] = 'Not-completed'
                    #     jsondata[subfolder]["ST"]["instrument"] = '-'
                    #     jsondata[subfolder]["ST"]["file count"] = '-'
                    #     jsondata[subfolder]["ST"]["polarimetry"] = '-'
                    #     jsondata[subfolder]["ST"]["manual-check"] = 'Folder is empty'
                    #     json.dump(jsondata, json_file, indent=4)

    print('clear')
main()
