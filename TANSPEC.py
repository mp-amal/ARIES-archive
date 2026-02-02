
import re
import os
import ephem
import paramiko
from scp import SCPClient
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import Angle, SkyCoord
from datetime import datetime, timedelta,date
import json
import random
import shutil 
import os,sys,psycopg2,glob
from datetime import datetime, timedelta,timezone,date

# # Path to the log file
TELESCOPE = 'DOT'
year = "2025"
cycle = "C2"

# -----------------------------------------------------------------------------
# NOTE:
# The paths below are DUMMY placeholders for sharing / public repos.
# Replace them with the real deployment paths on your archive server.
# -----------------------------------------------------------------------------

file_keywords   = '/path/to/TANSPEC/config/MasterHeaderList.dat'
file_telescopes = '/path/to/TANSPEC/config/TelescopeList.csv'   # Contains telescope details (dummy)

folder_path = '/path/to/DOT_LOG/tcs_aws_log'

rel_path           = "/path/to/archived_data/astro_data/final_data"
database_finaldata = "/path/to/archived_data/astro_data/final_data"   # Destination folder (dummy)

local_path = f'/path/to/DOT/TANSPEC/{year}-{cycle}/rawdata'            # Local destination (dummy)
final_path = f'/path/to/DOT/TANSPEC/{year}-{cycle}/ProcessedData'      # Local destination (dummy)

os.makedirs(local_path, exist_ok=True)

remote_path = f"/path/to/observation_Data/tanspec/{year}-{cycle}/"     # Remote folder (dummy)


# print(remote_path)
keywords_df = pd.read_csv(file_keywords, sep='\s+', comment='#').set_index('Header')        # Details of MaseterHeder file
keywords_df = keywords_df.fillna('NULL')
telescope_df = pd.read_csv(file_telescopes, sep=',', comment='#').set_index('ShortName')   

today_date = date.today().strftime("%Y%m%d")
yesterday = datetime.now() - timedelta(days=1)
yesterday =  yesterday.strftime("%Y%m%d")
print("\n",yesterday,'\n')
# yesterday = "20250424"
# input_date ='20241025'
server_ip = "000.00.0.00"  # Replace with your server's IP
username = "userid"    # Replace with your username
password = "password"    # Replace with your password



# propid ='P17'
# Extract year, month, and day
lyear = yesterday[:4]  # First 4 characters for the year
lmonth = yesterday[4:6]  # Characters 5 and 6 for the month
lday = yesterday[6:]  # Characters 7 and 8 for the day
# print(lyear)
# helper functions 

def PROPNO(file_name):
    pattern = r'P\d{1,2}'

    # Search for Proposal ID in the file name
    match = re.search(pattern, file_name)

    # Display the result
    if match:
        return match.group()
    else:
        print("No Proposal ID found in the file name.")

def copy_folder_from_server(server_ip, username, password):
    try:
        # Establish SSH connection
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(hostname=server_ip, username=username, password=password)
        # print('connected')
        stdin, stdout, stderr = ssh.exec_command('ls '+remote_path)
        folders = stdout.read().decode().splitlines()
        # print("Folders:", folders)
        if yesterday in folders:
            # print('YEs')
            if not os.path.exists(os.path.join(local_path,yesterday)):
                with SCPClient(ssh.get_transport()) as scp:
                    scp.get(remote_path+yesterday, local_path, recursive=True)
                    print(f"Folder copied successfully from {remote_path} to {local_path}")
            else:
                print("folder already copied")
    except Exception as e:
        print(f"Error: {e}")
    finally:
        ssh.close()

def parse_log_file_to_df(lfile_path):
    """
    Parses a log file and returns a DataFrame with structured data.

    Parameters:
        lfile_path (str): Path to the log file to be parsed.

    Returns:
        pd.DataFrame: A DataFrame containing the structured log data.
    """
    # Initialize an empty list for rows
    rows = []

    # Open the log file with the 'latin-1' encoding
    with open(lfile_path, "r", encoding="latin-1") as lfile:
        for line in lfile:
            if ":" in line:  # Ensure the line contains a timestamp
                # Split the line on ': ' to separate the timestamp from the rest
                try:
                    timestamp, details = line.split(": ", 1)
                    # Split the rest of the details into individual values
                    details_list = details.strip().split(",")
                    # Add the timestamp and details as a row
                    rows.append([timestamp] + details_list)
                except ValueError:
                    print(f"Skipping malformed line: {line.strip()}")

    # Define the columns based on the provided structure
    columns = [
        "Timestamp",
        "MsgNo",
        "Field1",
        "Field2",
        "AirTemp(degC)",
        "RH(%)",
        "Pressure(hPa)",
        "WindSpeed(m/s)",
        "WindDirection(deg)",
        "DewPoint(degC)",
        "Checksum"
    ]

    # Create a DataFrame
    df = pd.DataFrame(rows, columns=columns)

    return df

def merge_files_by_date(folder_path, input_date):
    """
    Merges all files with the same date into a single DataFrame and saves it as a CSV file.
    
    Args:
        folder_path (str): Path to the folder containing log files.
        input_date (str): Date in 'YYYY_MM_DD' format to filter files.
        output_suffix (str): Suffix for the merged output file name. Defaults to '_merged_data.csv'.
    
    Returns:
        str: Path to the merged output file or a message if no files found.
    """
    # Get a list of all log files in the folder
    files = [f for f in os.listdir(folder_path) if f.endswith('.log')]
    # print(files)
    # Filter files that match the input date
    matching_files = [f for f in files if '_'.join(f.split('_')[:3]) == lyear+'_'+lmonth+'_'+lday]
    print("\n",matching_files,"\n")
    if not len(matching_files) == 0:
        dfs =[]
        # # i=1
        for file in matching_files:
            file_path = os.path.join(folder_path, file)
            dfs.append(parse_log_file_to_df(file_path))

        log_data = pd.concat(dfs, ignore_index=True)
        # print(type(log_data))
    
        return log_data
    else:
        print("\n log files is missing \n")

def rename_files(filename, TFOLDER,final_path):
    final_path = os.path.join(final_path, TFOLDER)
    os.makedirs(final_path,exist_ok=True)
    
    data = fits.open(filename)
    header = data[0].header
    # image = data[0].data
    obj = header['OBJECT']
    date_obs = header['DATE-OBS']
    telescope = header['TELESCOP']
    instru = header['INSTRUME']
    # print(filename)
    if 'test' in filename.lower() or  'ts_' in os.path.basename(filename).lower():
        # print(filename)
        print("test")
    else:
        if any(substring in os.path.basename(filename).lower() for substring in ['ar', 'arg', 'ne', 'neo', 'cont']):
            # print(filename)
            code = 'L'
            new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
            print(new_name + ' -------------------------- ' + os.path.basename(filename))
            os.rename(filename, final_path + '/' + new_name)
        elif 'flat' in obj.lower() or 'flat' in filename.lower():
            code = 'F'
            new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
            print(new_name + ' -------------------------- ' + os.path.basename(filename))
            os.rename(filename, final_path + '/' + new_name)
        elif 'bias' in obj.lower() or 'bias' in filename.lower():
            code = 'F'
            new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
            print(new_name + ' -------------------------- ' + os.path.basename(filename))
        else:
            # header['CATEGORY']= 'Science'
            code = 'S'
            new_name=code+'-'+lyear+propid+'-'+date_obs+'-'+telescope+'-'+instru+'.fits'
            os.rename(filename, final_path + '/' + new_name)
            print(new_name+' -------------------------- '+os.path.basename(filename))

    # print(header['OBJECT'])
















    # if 'test' in obj.lower() or 'test' in filename.lower():
    #     # code = 'F'
    #     new_name = 'testfile'
    #     print("test file")
    #     # os.rename(filename, final_path + '/' + new_name)
    # if 'flat' in obj.lower() or 'flat' in filename.lower():
    #     code = 'F'
    #     new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
    #     print(new_name + ' -------------------------- ' + os.path.basename(filename))
    #     # os.rename(filename, final_path + '/' + new_name)
    # if 'bias' in obj.lower() or 'bias' in filename.lower():
    #     code = 'F'
    #     new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
    #     print(new_name + ' -------------------------- ' + os.path.basename(filename))
    #     # os.rename(filename, final_path + '/' + new_name)
    # for lamp in ['ar', 'arg','ne','neo','cont']:
    #     if lamp in filename.lower():
    # # any(lamp in obj.lower() for lamp in ['lamp', 'tanspec target']) or \
    #         code = 'L'
    #         new_name = f"{code}-{lyear}{propid}-{date_obs}-{telescope}-{instru}.fits"
    #         print(new_name + ' -------------------------- ' + os.path.basename(filename))
    #         # os.rename(filename, final_path + '/' + new_name)

    # else:
    #     # header['CATEGORY']= 'Science'
    #     code = 'S'
    #     new_name=code+'-'+lyear+propid+'-'+date_obs+'-'+telescope+'-'+instru+'.fits'
    #     print(new_name+' -------------------------- '+os.path.basename(filename))
    #     # os.rename(filename,final_path+'/'+new_name)
    # #     try:
    # #         os.rename(filename,final_path+'/'+new_name)
    # #     except FileNotFoundError:
    #         pass
    #     # print("\n"+"Calculation file renamed from ::: "+filename+" -----> "+new_name)
    # return new_name

def extract_environmental_data(file_path, log_data):
    # Access the primary header
    header = hdul[0].header
    fdate = header.get('DATE_OBS')  # Common FITS header key for date
    ftime = header.get('TIME_OBS')  # Sometimes TIME may be included in other headers (e.g., DATE-OBS)
    # print(f"Date: {fdate}")
    # print(f"Time: {ftime}")
    f =fdate+" "+ftime
    target_timestamp = datetime.strptime(f, "%Y-%m-%d %H:%M:%S.%f")

    # Define the tolerance in seconds
    tolerance = timedelta(seconds=1)

    # Find rows that match within the ±5-second tolerance
    matches = log_data[(log_data["Timestamp"] >= target_timestamp - tolerance) & 
                (log_data["Timestamp"] <= target_timestamp + tolerance)]
    # Check if there are any matches
    # print(matches)
    if not matches.empty:
    #     # print("Working")
    #     # print(matches)  # Optional: Display the matching rows
    #     # keywords = [
    #     #     "AirTemp(degC)", "RH(%)", "Pressure(hPa)",
    #     #     "WindSpeed(m/s)", "WindDirection(deg)", "DewPoint(degC)", "Checksum"]
    #     hdulist = fits.open(file_path, mode='update')
    #     header = hdulist[0].header

        AIRTEMP =str( matches["AirTemp(degC)"].values[0])
        RELHUM =matches["RH(%)"].values[0]
        PRESSURE = matches[ "Pressure(hPa)"].values[0]
        WINDSPD = matches["WindSpeed(m/s)"].values[0]
        WINDDIR = matches["WindDirection(deg)"].values[0]
        DEWPOINT = matches["DewPoint(degC)"].values[0]
        # CHECKSUM = matches["Checksum"].values[0]
        header['TC-AMBI'] = AIRTEMP
        header['RHUM'] = RELHUM

        header['WIND'] = WINDSPD
        header['WIND-DIR'] = WINDDIR
        # header['DEWPOINT'] = DEWPOINT
        # header['PRESSURE'] = PRESSURE
        # print("<----- Envionmental Header updated  ----->")
    #     # print("Air temprature : "+AIRTEMP)
    #     # print("Humidity: "+RELHUM)
    #     # print("Pressure : "+PRESSURE)
    #     # print("Wind speed : "+WINDSPD)
    #     # print("Wind Dir : "+WINDDIR)
    #     # print('DewPoint : '+DEWPOINT)
    #     # print("Checksum: "+CHECKSUM)                       
    # else:
    #     print("No matches found")
    # else:
        # print("<----- Envionmental details unavilable  ----->")
def append_nullheader(keywords_df,header):
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    # with fits.open(filename, mode='update') as hdulist:
    #     header = hdulist[0].header
        # header.remove('OBJECT', ignore_missing=True, remove_all=True) #neha
    for keyword in keywords_df.index:
        if keyword not in list(header.keys()):
            header.append(card=(keyword, keywords_df.loc[keyword, 'Value']))
    # print("<----- Null Header updated  ----->")
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
def telescope_det(TELESCOPE,telescope_df):
    
    _, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, HORIZON, ZENITH = telescope_df.loc[TELESCOPE].values
    # instrument_df = instrument_df.loc[TELESCOPE]

    dict_append = {'TELESCOP': TELESCOPE, 'GEOLONG': OBS_LONG, 'GEOLAT': OBS_LAT, 'GEOELEV': OBS_ALT,
                       'HORIZON': HORIZON, 'ZENITH': ZENITH, 'TIMEZONE': OBS_TIMEZONE, 'INSTRUME' : "TANSPEC"}
    for keyword, value in dict_append.items():
        header[keyword] = value
    # print("<----- telescope Header updated  ----->")
def calculate_airmass(TELESCOPE, filename):
    """
    Calculates AIRMASS for the FITS file and appends respective details in the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        filename      : FITS file whose header has to be edited
    Returns:
        None
    """
    # hdulist = fits.open(filename, mode='update')
    # header = hdulist[0].header

    date_obs = header['DATE-OBS']
    object_ra = header['RA']
    object_dec = header['DEC']
    
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
    airmass =round(airmass,3)
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

    # hdulist.flush()
    # hdulist.close()
    # print("\n"+"<----- .  Airmass calculation updated  ----->")
def get_month_abbreviation(month_value):
    # Dictionary mapping numerical month values to their abbreviations
    month_abbr = {
        "01": "Jan", "02": "Feb", "03": "Mar", "04": "Apr",
        "05": "May", "06": "Jun", "07": "Jul", "08": "Aug",
        "09": "Sep", "10": "Oct", "11": "Nov", "12": "Dec"
    }

    # Fetch the abbreviation from the dictionary
    return month_abbr.get(month_value, "Invalid month")

copy_folder_from_server(server_ip, username, password)
log_data = merge_files_by_date(folder_path,yesterday)
# print(log_data)
# if log_data is not None:
ldate =lyear+lmonth+lday
timestamps =log_data["Timestamp"] = pd.to_datetime(log_data["Timestamp"])
# print('\n',timestamps,'\n')
for dirpath, dirnames, filenames in os.walk(os.path.join(local_path,ldate)):
    for file in filenames:
        file_path = os.path.join(dirpath, file)
        # print(file_path)
        if not "test" in file_path.lower():
            # print(file)

            with fits.open(file_path, mode='update') as hdul:
                header = hdul[0].header
                # append_nullheader(keywords_df,header)
                propid =PROPNO(file)
                # print(propid)
                if not propid == None:
                    header['PROPNO'] = propid
                else:
                    header['PROPNO'] =""
                keywords = ('lamp', 'ar_', 'arg_', 'ne_', 'neo_', 'cont', 'flat', 'dark', 'bias')
                if any(kw in file_path.lower() for kw in keywords):
                    header['CATEGORY'] = 'Calib'
                else:
                    header['CATEGORY'] = 'Science'
                header['DATE-OBS'] = header['DATE_OBS']+'T'+header['TIME_OBS']
                header['RA']= str(float(header['A_TRGTRA'])*15)
                header['DEC'] = header['A_TRGTDE']
                telescope_det(TELESCOPE,telescope_df)
                calculate_airmass(TELESCOPE, file_path)
                # extract_environmental_data(file_path, log_data)
                header['ORIGFILE'] = file
                header['EXPTIME'] =header['ITIMEREQ']
                header['GRISM'] = (header['GRATING'], 'GRATING')
                header['FILTER1'] = header['FILTER']
                header['FILTER2'] = 'NA'
                if propid == None:
                    # print(file)
                    propid ='PXX'
            # print(header['DATE-OBS'])
            rename_files(file_path,ldate,final_path)



    target_dir = os.path.join(final_path, ldate)
    # List only .fits files that are actual files
    fits_files = [f for f in os.listdir(target_dir) if f.endswith('.fits') and os.path.isfile(os.path.join(target_dir, f))]
    # print(fits_files)

    random_file = random.choice(fits_files)
    print(random_file)
    random_split = random_file.split('-')
    
    year = random_split[2]
    mon =random_split[3]
    month =get_month_abbreviation(mon)
    day =random_split[4].split('T')[0]
    date_folder=os.path.basename(str(target_dir))
    print(date_folder)
    telescope =random_split[5] 
    instrument =random_split[6].split('.')[0] 
    # print(instrument)
    db_dest_path = database_finaldata+'/'+str(year)+"/"+str(month)+"/"+str(telescope)+"/"+str(instrument)+"/"+str(date_folder)
    # print(db_dest_path)

    if os.path.exists(db_dest_path):
        print('Destination exist with path : '+db_dest_path)
    else:
        shutil.copytree(target_dir,db_dest_path)
        print(target_dir+ "--------------> copied to db")



    file_path = (os.path.relpath(db_dest_path, rel_path)) # Relative file path to be inserted in the database
    # print(file_path)
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
        connection = psycopg2.connect(user="userid",
                                    password="password",
                                    host="host ip",
                                    port="5432",
                                    database="db")
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
            
            TIMEINSERTED = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%s")				
            # print(TIMEINSERTED)
                    
            datarow=(DATEOBS,OBJECT,RA,DEC,PROPNO,PROG,OBSERVER,TELESCOP,INSTRUMENT,CATEGORY,TYPE,MODE,DATAID,ORIGFILE,EXPOSURE,FILTER1,FILTER2,GRISM,SLIT,AIRMASS,FILESIZE,FILEPATH,TIMEINSERTED)
            # print(datarow)
                    
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
    TIMEINSERTED = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S.%s")
    masterlog.write('{}\t {}\t {}\t {}\t \n'.format(TIMEINSERTED,db_dest_path,nofp,totfits))
    if connection:
        cursor.close()
        connection.close


        
    infile.close()
else:
    print("\n Check the log files are syncing")
            