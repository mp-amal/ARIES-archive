import numpy as np
from astropy.time import Time, TimeDelta
from astropy.table import Table, Column, hstack, vstack 
from astropy.io import fits
import glob, yaml
import os,sys
import time
import astropy.units as u
import argparse
from astropy.coordinates import ICRS, SkyCoord, EarthLocation, Angle
from astroquery.simbad import Simbad

# sys.path.insert(0,'/opt/anaconda3/lib/python3.8/site-packages')

master_header = "./ARIES-archive/scripts/ADFOSC-ARIES-main_V1.2/Master_headers.yaml"
site_telinfo = "./ARIES-archive/scripts/ADFOSC-ARIES-main_V1.2/site_telescope_info.yaml" 

parameters = {
    "NB": {
        0: "clear",
        1: "Small aperture",
        2: "H-beta",
        3: "BG3",
        4: "Low red H-alpha",
        5: "High red H-alpha",
        6: " ",
        7: " ",
        8: " ",
        9: " "
    },
    "SW": {
        0: "clear",
        1: "1.0'",  #"Brass large",
        2: "1.5'",  #"SS",
        3: "2.0'",  #"SS",
        4: "3.2'",  #"SS narrow",
        5: "Slit holder",
        6: " ",
        7: " ",
        8: " ",
        9: " "
    },
    "BB": {
        0: "z",
        1: "WP",
        2: "u",
        3: "BG39",
        4: "g",
        5: "RG610",
        6: "r",
        7: "Clear",
        8: "i",
        9: "Clear"
    },
    "GW": {
        0: "Clear",
        1: "770R-300g/mm",
        2: "132R-600g/mm",
        3: "460R-830.7g/mm",
        4: "676R-420g/mm",
        5: "Prism",
        6: "DIMM prism",
        7: "Clear",
        8: "100g/mm",
        9: "Clear"
    }
}

name_map = {
'Broadband position': 'BB',
'Narrowband position': 'NB',
'Grism wheel position': 'GW', 
'Slit wheel position': 'SW'
}

# The codes written by Vivek Jha and Sunil Chandra 
def ConvSysTime2Isot(date_string=None, time_string=None):
    if date_string and time_string:
        date_parts = date_string.split('/')
        if len(date_parts) > 2:
            year, month, day = date_parts[2], date_parts[1], date_parts[0]
            datetime_string = f"{year}-{month}-{day}T{time_string}"
            t = Time(datetime_string, format='isot', scale='local')
            return t
    print("ERROR: Invalid date & Time string")
    return None
    
def read_adfosc(filename = None):
        
    file1 = open(filename, 'r')
    lines = file1.readlines()
    count = 0
    broad_band_fltr_data = ([])
    narrow_band_fltr_data = ([])
    grism_pos_data = ([])
    slit_pos_data = ([])
    for line in lines:
        count += 1
        if line.strip != 0:
            components = line.strip().split()
            if len(components) > 2: 
                date_str = components[0]
                time_str = components[1]
                time2use = ConvSysTime2Isot(date_string = date_str, time_string = time_str)
                if time2use != None:
                    timeisot = time2use.isot
                    timemjd = time2use.mjd
                else:
                    timeisot = time2use
                    timemjd = time2use
                position = " ".join(components[2:-3])
                
                mapped_position = name_map.get(position, position)
                #if mapped_position:
                number = components[-1]
                if (mapped_position == "BB") and (timeisot != None):
                    if len(components) == 7:
                        broad_band_fltr_data.append((timeisot, timemjd, 
                                                     parameters[mapped_position][int(number)], 
                                                     number))
                elif (mapped_position == "NB") and (timeisot != None):
                    if len(components) == 7:
                        narrow_band_fltr_data.append((timeisot, timemjd, 
                                                     parameters[mapped_position][int(number)], 
                                                     number))
                elif (mapped_position == "GW") and (timeisot != None):
                    if len(components) == 8:
                        grism_pos_data.append((timeisot, timemjd, 
                                               parameters[mapped_position][int(number)], 
                                               number))
                elif (mapped_position == "SW") and (timeisot != None):
                    if len(components) == 8:
                        slit_pos_data.append((timeisot, timemjd, 
                                                parameters[mapped_position][int(number)], 
                                                number))
                else:
                    continue    
    broad_band_fltr_data = np.array(broad_band_fltr_data)
    narrow_band_fltr_data = np.array(narrow_band_fltr_data)
    grism_pos_data = np.array(grism_pos_data)
    slit_pos_data = np.array(slit_pos_data)
    #meta_data2ret = np.column_stack((broad_band_fltr_data, narrow_band_fltr_data, grism_pos_data, slit_pos_data))
    return broad_band_fltr_data, narrow_band_fltr_data, grism_pos_data, slit_pos_data
    

def read_telnet(filename = None):

    file1 = open(filename, 'r')
    lines = file1.readlines()
    count = 0
    tcsdata = ([])
    tcsmoutdata = ([])
    tcsweatherdata = ([])
    for line in lines:
        count += 1
        #print("Line{}: {}".format(count, line.strip()))
        if line.strip()!= "":
            components = line.strip().split()
            date_str = components[0]
            time_str = components[1]
            time2use = ConvSysTime2Isot(date_string = date_str, time_string = time_str)
            if time2use != None :
                timeisot = time2use.isot
                timemjd = time2use.mjd
            else :
                timeisot = time2use
                timeisot = None     		## DNT KNOW WHY DEFINING timeisot AGAIN AS None HERE. (MIGHT BE timemjd)
            #Saving data into lists
            if 'tcsdata' in components:
                mjd1, targetra, targetdec, airmass = components[6:]
                tcsdata.append((timeisot, timemjd, targetra, targetdec, airmass))
            elif 'tcsmount' in components:
                mjd2, altitude, azimuth, rotator = components[6:]
                tcsmoutdata.append((timeisot, timemjd, altitude, azimuth, rotator))
            elif 'tcsweather' in components:
                mjd3, temperature, humidity = components[6:]
                tcsweatherdata.append((timeisot, timemjd, temperature, humidity))
            else:
                continue
                
    tcsdata = np.array(tcsdata)
    tcsmoutdata = np.array(tcsmoutdata)
    tcsweatherdata = np.array(tcsweatherdata)
    return tcsdata, tcsmoutdata, tcsweatherdata


def time_from_file(filename = None):
    print(f"name of the fits file :{filename}")
    data=fits.open(filename)
    header=data[0].header
    time=header['DATE-OBS']
    data.close()
    return time


def create_logs_trange(log_path ="./", time=None, ttol=1):
    print(f' the log path : {log_path}')   ##JIN
    def convertisot2fname(time):
        # This returns the name of the file pattern up to hrs level
        t = Time(time, format='isot', scale='local')
        t2ret = "_".join(t.isot.replace('-', '_').replace(":", "_").split("_")[0:-2])
        return t2ret.replace('T', "_")
    
    #print(f'{time} in create_logs_trange function | before passing to Time function')     ## JIN
    time = Time(time, scale='local', format='isot')
    print(f'{time} in create_logs_trange function | after passing to Time function')     ## JIN
    delta = TimeDelta(ttol * u.minute)
    
    start_time = time - delta
    stop_time = time + delta
    #print(f'start and end time : {start_time} and {stop_time} | delta time: {delta} in create_log_trange')  ##JIN
    search_str1 = convertisot2fname(start_time.isot)
    search_str2 = convertisot2fname(stop_time.isot)
    #print(f' search string for start and stop time for the log files are {search_str1} & {search_str2}')	##JIN
    adfosc_log1 = glob.glob(os.path.join(log_path, f"adfosc_ics_log*{search_str1}*.txt"))
    adfosc_log2 = glob.glob(os.path.join(log_path, f"adfosc_ics_log*{search_str2}*.txt"))
    tcstel_log1 = glob.glob(os.path.join(log_path, f"tcs_telnet_log*{search_str1}*.txt"))
    tcstel_log2 = glob.glob(os.path.join(log_path, f"tcs_telnet_log*{search_str2}*.txt"))
    adfosc_log = adfosc_log1 + adfosc_log2
    tcstel_log = tcstel_log1 + tcstel_log2
    #print(f' looking for the files in format {os.path.join(log_path, f"adfosc_ics_log*{search_str1}*.txt")}')
    #print(f'the Number of log files for create_log_trange are {len(adfosc_log)} and {len(tcstel_log)}')		##JIN
    return adfosc_log, tcstel_log

def create_logs_all(log_path = "./"):
    def convertname2isot(fname):
        flist = os.path.basename(fname).replace(".txt","").split('_')[-6:]
        datevec = flist[0:3]
        timevec = flist[3:]
        stro = f"{'-'.join(datevec)}T{':'.join(timevec)}"
        t = Time(stro, format='isot', scale = 'local')
        return stro, t.mjd 
    def extractTimefromfpath(fname): 

        def run_forNonStandardNameOfFile(filename):
            fc = open(filename, 'r')
            linedata = fc.readlines()
            for line in linedata:
                if line.strip() != "":
                    first_entry = line.split()
                    if len(first_entry) > 2:
                        first_entry = first_entry[0:2]
                        timeconv = ConvSysTime2Isot(date_string = first_entry[0], 
                                                        time_string = first_entry[1])
                        return timeconv
        outvec = ([])   
        for filename in fname:
            filebasename = os.path.basename(filename)
            if len(filebasename.split("_")) < 7:
                timeconv = run_forNonStandardNameOfFile(filename)       
                if timeconv != None:
                    timeisot = timeconv.isot
                    timemjd =  timeconv.mjd
                else:
                    timeisot = "1995-01-01T00:00:00"
                    timemjd =  49718.0
                timeinfo = timeisot, timemjd  
            else:
                timeinfo = convertname2isot(filename)
            outvec.append([filename, timeinfo[0], timeinfo[1]])
        return np.array(outvec)    

    #Extracting the logfile from data base path 
    adfosc_log = sorted(glob.glob(os.path.join(log_path, '*adfosc_ics*.txt')))
    telnet_log = sorted(glob.glob(os.path.join(log_path, '*tcs_telnet*.txt')))
    #print(f' log files for the information extraction - {adfosc_log} and {telnet_log} in create_logs_all function')  ##JIN
    dbadfosc = extractTimefromfpath(adfosc_log)
    dbtelnet = extractTimefromfpath(telnet_log)
    return dbadfosc, dbtelnet

def applyFileTime2getdb(intime = None, logdb = None, tolerancetime = [10., 10.]):
    
    adfosc_toltime, tcs_toltime = tolerancetime
    def filterdbsOntolerance(intime, logdb, toltime = tcs_toltime):
        #intime should be in ISOT format 
        # tolerence should be in minutes
        timeformat = Time(intime, format='isot', scale='local')
        delta = TimeDelta(toltime * u.minute)
        time_start_w_tol = timeformat - delta
        time_stop_w_tol = timeformat + delta
        mjd_start_w_tol = time_start_w_tol.mjd
        mjd_stop_w_tol = time_stop_w_tol.mjd
        logdb_fltrd = logdb[(np.array(logdb[:,-1], float) >= mjd_start_w_tol) & 
                                  (np.array(logdb[:,-1], float) <= mjd_stop_w_tol)]

        return logdb_fltrd
    
    adfoscdb, tcstelnetdb = logdb
    adfoscdb_fltrd = filterdbsOntolerance(intime, adfoscdb, toltime = adfosc_toltime)
    while len(adfoscdb_fltrd) < 2 :
        tolerancetime = tolerancetime * 0.25
        adfoscdb_fltrd = filterdbsOntolerance(intime, adfoscdb, toltime = adfosc_toltime)

    tcstelnetdb_fltrd = filterdbsOntolerance(intime, tcstelnetdb, toltime = tcs_toltime)    
    while len(tcstelnetdb_fltrd) < 2 :
        tolerancetime = tolerancetime * 0.25
        tcstelnetdb_fltrd = filterdbsOntolerance(intime, tcstelnetdb, toltime = tcs_toltime)

    return adfoscdb_fltrd, tcstelnetdb_fltrd    


#Optimized version of above function using the bing AI

def read_adfosc_v02(filename = None):
    with open(filename, 'r') as file1:
        lines = file1.readlines()
    broad_band_fltr_data = []
    narrow_band_fltr_data = []
    grism_pos_data = []
    slit_pos_data = []
    # Code to process lines and append data to the lists above
    for line in lines:
        if line.strip():
            components = line.strip().split()
            if len(components) > 2:
                date_str, time_str = components[0], components[1]
                time2use = ConvSysTime2Isot(date_string=date_str, time_string=time_str)
                if time2use:
                    timeisot, timemjd = time2use.isot, time2use.mjd
                else:
                    timeisot = timemjd = time2use
                position = " ".join(components[2:-3])
                mapped_position = name_map.get(position, position)
                number = components[-1]
                if mapped_position == "BB" and timeisot and len(components) == 7:
                    broad_band_fltr_data.append(
                        (timeisot, timemjd, parameters[mapped_position][int(number)], number))
                elif mapped_position == "NB" and timeisot and len(components) == 7:
                    narrow_band_fltr_data.append(
                        (timeisot, timemjd, parameters[mapped_position][int(number)], number))
                elif mapped_position == "SW" and timeisot and len(components) == 8:
                    slit_pos_data.append(
                        (timeisot, timemjd, parameters[mapped_position][int(number)], number))    
                elif mapped_position == "GW" and timeisot and len(components) == 8:
                    grism_pos_data.append(
                        (timeisot, timemjd, parameters[mapped_position][int(number)], number))
                    
    broad_band_fltr_data = np.array(broad_band_fltr_data)
    narrow_band_fltr_data = np.array(narrow_band_fltr_data)
    grism_pos_data = np.array(grism_pos_data)
    slit_pos_data = np.array(slit_pos_data)
    return broad_band_fltr_data, narrow_band_fltr_data, grism_pos_data, slit_pos_data  


def read_telnet_v02(filename=None):
    with open(filename, 'r') as file1:
        lines = file1.readlines()
    tcsdata = []
    tcsmoutdata = []
    tcsweatherdata = []
    for line in lines:
        if line.strip() != "":
            components = line.strip().split()
            date_str = components[0]
            time_str = components[1]
            time2use = ConvSysTime2Isot(date_string=date_str, time_string=time_str)
            if time2use is not None:
                timeisot = time2use.isot
                timemjd = time2use.mjd
            else:
                timeisot = None
                timemjd = None
            if 'tcsdata' in components:
                mjd1, targetra, targetdec, airmass = components[6:]
                tcsdata.append((timeisot, timemjd, targetra, targetdec, airmass))
            elif 'tcsmount' in components:
                mjd2, altitude, azimuth, rotator = components[6:]
                tcsmoutdata.append((timeisot, timemjd, altitude, azimuth, rotator))
            elif 'tcsweather' in components:
                mjd3, temperature, humidity = components[6:]
                tcsweatherdata.append((timeisot, timemjd, temperature, humidity))
    tcsdata = np.array(tcsdata)
    tcsmoutdata = np.array(tcsmoutdata)
    tcsweatherdata = np.array(tcsweatherdata)
    return tcsdata, tcsmoutdata, tcsweatherdata

'''   
def log_file_used(filename,logfile):
    filename = os.path.basename(filename)
    adfosc_ics = "adfosc_ics_"+str(filename).replace((".fit" or ".fits"), ".txt")
    with (adfosc_ics,"w") as ics_log:
        log.write(f'{filename}:{logfile}\n')

'''

def DataBasefromFile_v02(adfosclist, telnetlist):
    adfoscmetadata = []
    used_adfosc = []
    unused_adfosc = []

    # Process ADFOSC files
    for adfoscfile in adfosclist:
        #print(f'The adfosclogfile for HDR inform: {adfoscfile}')  #JIN
        try:
            if os.path.getsize(adfoscfile) != 0:
                adfoscdata = read_adfosc(adfoscfile)
                used_adfosc.append(adfoscfile)  # Mark as used

                if not adfoscmetadata:
                    adfoscmetadata = adfoscdata
                else:
                    #print(f'adfoscmetadata: {len(adfoscmetadata)},\n adfoscdata: {len(adfoscdata)}')  #JIN
                    adfoscmetadata = [np.concatenate((adfoscmetadata[i], adfoscdata[i])) 
                                      for i in range(len(adfoscdata))]
            else:
                unused_adfosc.append(adfoscfile)  # File is empty, mark as unused
        except Exception as e:
            print(f'Error in file {adfoscfile}: {e}')
            unused_adfosc.append(adfoscfile)  # Mark as unused in case of exception

    # Save used and unused ADFOSC files
    with open(f'./ARIES-archive/logs/used_adfosc.txt', 'w') as used_file:
        for file in used_adfosc:
            used_file.write(f'{file}\n')
    with open(f'.ARIES-archive/logs/unused_adfosc.txt', 'w') as unused_file:
        for file in unused_adfosc:
            unused_file.write(f'{file}\n')

    tcsmetadata = []
    used_telnet = []
    unused_telnet = []

    # Process Telnet files
    for telnetfile in telnetlist:
        try:
            if os.path.getsize(telnetfile) != 0:
                tcsdata = read_telnet(telnetfile)
                used_telnet.append(telnetfile)  # Mark as used

                if not tcsmetadata:
                    tcsmetadata = tcsdata
                else:
                    tcsmetadata = [np.concatenate((tcsmetadata[i], tcsdata[i])) for i in range(len(tcsdata))]
            else:
                unused_telnet.append(telnetfile)  # File is empty, mark as unused
        except Exception as e:
            print(f'Error in file {telnetfile}: {e}')
            unused_telnet.append(telnetfile)  # Mark as unused in case of exception

    # Save used and unused Telnet files
    # with open('used_telnet.txt', 'a') as used_file:
    #     for file in used_telnet:
    #         used_file.write(f'{file}\n')
    # with open('unused_telnet.txt', 'a') as unused_file:
    #     for file in unused_telnet:
    #         unused_file.write(f'{file}\n')

    return adfoscmetadata, tcsmetadata


'''
def DataBasefromFile_v02(adfosclist, telnetlist):
    adfoscmetadata = []
    for adfoscfile in adfosclist:
        print(f' the adfosclogfile for HDR inform: {adfoscfile}')  #JIN
        try:           
            if os.path.getsize(adfoscfile) != 0:
                
                adfoscdata = read_adfosc(adfoscfile)
                if not adfoscmetadata:
                    adfoscmetadata = adfoscdata
                else:
                    print(f' adfoscmetadata : {len(adfoscmetadata)},\n adfoscdata : {len(adfoscdata)}')  #JIN
                    #print(f' adfoscmetadata : {adfoscmetadata},\n adfoscdata : {adfoscdata}')  #JIN
                    adfoscmetadata = [np.concatenate((adfoscmetadata[i], adfoscdata[i])) 
                                      for i in range(len(adfoscdata))]
        except Exception as e:
            print(f' error in file {adfoscfile} {e}')
    tcsmetadata = []
    for telnetfile in telnetlist:
        if os.path.getsize(telnetfile) != 0:
            tcsdata = read_telnet(telnetfile)
            if not tcsmetadata:
                tcsmetadata = tcsdata
            else:
                tcsmetadata = [np.concatenate((tcsmetadata[i], tcsdata[i])) for i in range(len(tcsdata))]            
    return adfoscmetadata, tcsmetadata
'''

"""
##original
def DataBasefromFile_v02(adfosclist, telnetlist):
    adfoscmetadata = []
    for adfoscfile in adfosclist:
        
        
        if os.path.getsize(adfoscfile) != 0:
            adfoscdata = read_adfosc(adfoscfile)
            if not adfoscmetadata:
                adfoscmetadata = adfoscdata
            else:
                adfoscmetadata = [np.concatenate((adfoscmetadata[i], adfoscdata[i])) 
                                  for i in range(len(adfoscdata))]
    tcsmetadata = []
    for telnetfile in telnetlist:
        if os.path.getsize(telnetfile) != 0:
            tcsdata = read_telnet(telnetfile)
            if not tcsmetadata:
                tcsmetadata = tcsdata
            else:
                tcsmetadata = [np.concatenate((tcsmetadata[i], tcsdata[i])) for i in range(len(tcsdata))]            
    return adfoscmetadata, tcsmetadata

"""
def create_adfosc_logs_all_v02(log_path='./'):
    def convertname2isot(fname):
        flist = os.path.basename(fname).replace(".txt", "").split('_')[-6:]
        datevec = flist[0:3]
        timevec = flist[3:]
        stro = f"{'-'.join(datevec)}T{':'.join(timevec)}"
        t = Time(stro, format='isot', scale='local')
        return stro, t.mjd

    def run_forNonStandardNameOfFile(filename):
        with open(filename, 'r') as fc:
            linedata = fc.readlines()
        for line in linedata:
            if line.strip():
                first_entry = line.split()
                if len(first_entry) > 2:
                    first_entry = first_entry[0:2]
                    timeconv = ConvSysTime2Isot(date_string=first_entry[0], 
                                                time_string=first_entry[1])
                    return timeconv
    def extractTimefromfpath(fname):
        outvec = []
        for filename in fname:
            filebasename = os.path.basename(filename)
            if len(filebasename.split("_")) < 7:
                timeconv = run_forNonStandardNameOfFile(filename)
                if timeconv:
                    timeisot, timemjd = timeconv.isot, timeconv.mjd
                else:
                    timeisot, timemjd = "1995-01-01T00:00:00", 49718.0
                timeinfo = timeisot, timemjd
            else:
                timeinfo = convertname2isot(filename)
            outvec.append([filename, timeinfo[0], timeinfo[1]])
        return np.array(outvec)

    adfosc_log = sorted(glob.glob(os.path.join(log_path, '*adfosc_ics*.txt')))
    telnet_log = sorted(glob.glob(os.path.join(log_path, '*tcs_telnet*.txt')))
    #print(f' log files for the information extraction - {adfosc_log} and {telnet_log} in create_logs_all_V02 function')  ##JIN
    dbadfosc = extractTimefromfpath(adfosc_log)
    dbtelnet = extractTimefromfpath(telnet_log)
    return dbadfosc, dbtelnet

def applyFltrTime2getdb_v02(intime = '2023-01-02T23:59:59', logdb = None, tolerancetime = 50):
    #The funtion to filter the database based on input DATE-OBS input
    def filterdbsOntolerance(intime, tolerancetime):
        timeformat = Time(intime, format='isot', scale='local')
        delta = TimeDelta(tolerancetime * u.minute)
        time_start_w_tol = timeformat - delta
        time_stop_w_tol = timeformat + delta
        mjd_start_w_tol = time_start_w_tol.mjd
        mjd_stop_w_tol = time_stop_w_tol.mjd
        adfoscdb, tcstelnetdb = logdb
        adfoscdb_fltrd = adfoscdb[(np.array(adfoscdb[:, -1], float) >= mjd_start_w_tol) &
                                  (np.array(adfoscdb[:, -1], float) <= mjd_stop_w_tol)]
        tcstelnetdb_fltrd = tcstelnetdb[(np.array(tcstelnetdb[:, -1], float) >= mjd_start_w_tol) &
                                        (np.array(tcstelnetdb[:, -1], float) <= mjd_stop_w_tol)]
        return adfoscdb_fltrd, tcstelnetdb_fltrd
    #Running the funtions with input variables 
    adfoscdb_fltrd, tcstelnetdb_fltrd = filterdbsOntolerance(intime, tolerancetime)
    while len(adfoscdb_fltrd) < 2 and len(tcstelnetdb_fltrd) < 2:
        tolerancetime *= 0.25
        adfoscdb_fltrd, tcstelnetdb_fltrd = filterdbsOntolerance(intime, tolerancetime)
    return adfoscdb_fltrd[:,0], tcstelnetdb_fltrd[:,0]

def GetMetric4romDb(intime = '2023-01-02T23:59:59', ttol = 10, db = None):
    adfocsdb, tcsdb = db
    t = Time(intime, format='isot', scale='local')
    delta = TimeDelta(ttol * u.minute)
    time_start = t - delta
    time_stop = t + delta
    mjd_start = time_start.mjd
    mjd_stop = time_stop.mjd
    
    def filterdb(db, mjdrange = [57200.000, 57800.00]):
        dbo = []
        for i in range(len(db)):
            db2use = db[i]
            #print (f"MJD RANGE INFO: {mjdrange}")
            #print (f"DB-PRINT INFO: {np.min(abs(np.array(db2use[:,1], float) - mjdrange[0])*24*60)}")
            db2use_fltrd = db2use[(np.array(db2use[:,1], float) >= mjdrange[0]) & (np.array(db2use[:,1], float) <= mjdrange[1]),]
            dbo.append(db2use_fltrd)
        return dbo  
    adfocsdb_fltrd =  filterdb(adfocsdb, mjdrange = [mjd_start, mjd_stop])
    tcsdb_fltrd =  filterdb(tcsdb, mjdrange = [mjd_start, mjd_stop])

    return adfocsdb_fltrd, tcsdb_fltrd


def metric_sort_down(intime = '2023-01-02T23:59:59', metricdb = None):
    
    def GetMinDiff(metricdb = None):
        nearestmetric = []
        for i in range(len(metricdb)):
            metric2use = metricdb[i]
            mjddiffvec = abs(np.array(metric2use[:,1], float) - mjd_image)
            final_metric4head = metric2use[mjddiffvec == np.min(mjddiffvec)]
            nearestmetric.append(final_metric4head)           
        return nearestmetric  

    def MakeFinalTable(metricdb = None, dbtype = 'adfosc'):
    
        if dbtype.upper() in ['ADFOSC', 'A', 'AD', 'ADF', 'ADFC']:
            for i in range(len(metricdb)):
                if i < 1:
                    metricvec4o = metricdb[i]
                else:
                    metricvec4o = np.concatenate((metricvec4o, metricdb[i][2:3]))
            #writing these as header info        
            headerinfo = ['TIME', 'MJD', 'FILTER1', 'FILTER2', 'GRISM', 'SLIT']
            for j in range(len(metricvec4o)):        
                metricvec4t = Table()
                metricvec4t.add_column(Column([metricvec4o[j]], names = headerinfo[j]))

        elif dbtype.upper() in ['TCS', 'T', 'TELESCOPE', 'TEL', 'SCOPE']:
            for i in range(len(metricdb)):
                if i < 1:
                    metricvec4o = metricdb[i]
                else:
                    metricvec4o = np.concatenate((metricvec4o, metricdb[i][2:]))
            #writing these as header info        
            headerinfo = ['TIME', 'MJD', 'RA', 'DEC', 'AIRMASS', 'ALTITUDE', 'AZIMUTH', 
                          'ROTANGLE', 'TC-AMBI', 'RHUM'] 
            for j in range(len(metricvec4o)):       
                metricvec4t = Table()
                metricvec4t.add_column(Column([metricvec4o[j]], names = headerinfo[j]))
        else:
            print ("ERROR: Unknown metric type!")
            sys.exit()
        return metricvec4t  
      
    t = Time(intime, format='isot', scale='local')
    mjd_image = t.mjd

    adfoscmetric, tcsmetric = metricdb
    adfc_metric4head = GetMinDiff(metricdb = adfoscmetric)
    tcs_metric4head = GetMinDiff(metricdb = tcsmetric)
    
    adfc_metric4head_t =  MakeFinalTable(metricdb = adfc_metric4head, dbtype = 'adfosc')
    tcs_metric4head_t =  MakeFinalTable(metricdb = tcs_metric4head, dbtype = 'tcs')
    final_metric_tbl = hstack([adfc_metric4head_t, tcs_metric4head_t]) 
    return adfc_metric4head, tcs_metric4head, final_metric_tbl

def metric_sort_down_v02(intime='2023-01-02T23:59:59', metricdb=None):
    def GetMinDiff(metricdb=None):
        nearestmetric = []
        for l, metric2use in enumerate(metricdb):
            mjddiffvec = abs(np.array(metric2use[:, 1], float) - mjd_image)
            final_metric4head = metric2use[mjddiffvec == np.min(mjddiffvec)]
            nearestmetric.append(final_metric4head[0])
        return nearestmetric

    def MakeFinalTable(metricdb4tbl=None, dbtype='adfosc', metricdicto = None):

        if dbtype.upper() in ['ADFOSC', 'A', 'AD', 'ADF', 'ADFC']:
            for i, metric in enumerate(metricdb4tbl):
                if i < 1:
                    metricvec4o = metric[:-1]
                else:
                    metricvec4o = np.concatenate((metricvec4o, metric[2:3])) 
                    
            headerinfo = ['TIME', 'MJD', 'FILTER1', 'FILTER2', 'GRISM', 'SLIT']
            #print(f' metricvec4o is {metricvec4o}')      ##JIN
            metricvec4t = Table()
            for j, header in enumerate(headerinfo):
                metricvec4t.add_column(Column([metricvec4o[j]], name=header))
                #print(f' updatation of header for {header} = {metricvec4o[j]}')    #JIN
                metricdicto[header] = metricvec4o[j]
        elif dbtype.upper() in ['TCS', 'T', 'TELESCOPE', 'TEL', 'SCOPE']:
            for i, metric in enumerate(metricdb4tbl):
                if i < 1:
                    metricvec4o = metric
                else:
                    metricvec4o = np.concatenate((metricvec4o, metric[2:]))         
            headerinfo = ['TIME', 'MJD', 'RA', 'DEC', 'AIRMASS', 'ALTITUDE', 'AZIMUTH',
                          'ROTANGLE', 'TC-AMBI', 'RHUM']
            #print(f' metricvec4o is {metricvec4o}')      ##JIN
            metricvec4t = Table()
            for j, header in enumerate(headerinfo):
            
                if header == 'RA':
                    ra_log = metricvec4o[j]
                    ra_degree = float(ra_log)*15
                    metricvec4t.add_column(Column([str(ra_degree)], name='RA')) 
                else:             
                    metricvec4t.add_column(Column([metricvec4o[j]], name=header))
                #print(f' updatation of header for {header} = {metricvec4o[j]}')    #JIN
                metricdicto[header] = metricvec4o[j]
        else:
            print("ERROR: Unknown metric type!")
            sys.exit()
            
        #print(f" final header information {metricvec4t}")
        return metricvec4t, metricdicto

    t = Time(intime, format='isot', scale='local')
    mjd_image = t.mjd
    adfc_metricdb4src, tcs_metricdb4src = metricdb
    adfoscmetric = GetMinDiff(metricdb=adfc_metricdb4src)
    tcsmetric = GetMinDiff(metricdb=tcs_metricdb4src)
    metricdicto = {}
    adfc_metric4head, metricdicto = MakeFinalTable(metricdb4tbl=adfoscmetric, 
                                                   dbtype='adfosc', metricdicto = metricdicto)
    tcs_metric4head, metricdicto = MakeFinalTable(metricdb4tbl=tcsmetric, 
                                                  dbtype='tcs', metricdicto = metricdicto)
    #print (adfc_metric4head)
    #print (tcs_metric4head)
    final_metric_tbl = hstack((adfc_metric4head, tcs_metric4head))
    return adfc_metric4head, tcs_metric4head, final_metric_tbl, metricdicto

def run_proc(infile = "infile.fits", log_path="./", make_list_mode = 0, listfltr_ttol = 50, 
             obstime_ttol = 10):
    #extarct the DATE-OBS from the input file
    dateobs = time_from_file(infile)
    
    if make_list_mode == 0:
        #Create the log lists using date based formatting - min tolrance is one hour
        adfosc_loglist, tcstel_loglist = create_logs_trange(log_path = log_path, time = dateobs, 
                                                            ttol = listfltr_ttol)
        adfosclist, telnetlist = adfosc_loglist, tcstel_loglist
    else:
        #Create the log lists using date based formatting - min tolrance is one hour
        adfosc_loglist, tcstel_loglist = create_adfosc_logs_all_v02(log_path = log_path)
        
        #Applying the filter on loglist based on the DATE-OBS from the input fits file
        logdb = adfosc_loglist, tcstel_loglist
        adfosclist, telnetlist = applyFltrTime2getdb_v02(intime = dateobs, logdb = logdb, 
                                                         tolerancetime = listfltr_ttol)
       
    adfoscdb_fltrd, tcstelnetdb_fltrd = DataBasefromFile_v02(adfosclist, telnetlist)
    #Get Metrics from the database to match with the DATE-OBS from the input fits file and ttol
    db = adfoscdb_fltrd, tcstelnetdb_fltrd
    
    adfocsmetric, tcsmetric = GetMetric4romDb(intime = dateobs, ttol = obstime_ttol, db = db)
    #print (adfocsmetric, tcsmetric)
    metricdb = adfocsmetric, tcsmetric
    adfc_metric4head, tcs_metric4head, final_metric_tbl, metricdicto = metric_sort_down_v02(intime = dateobs, 
                                                                                            metricdb = metricdb)
                                                                                           
    #	print(f' the outout TCS information from the logfile is \n {tcs_metric4head} \n and the final updated information \n {final_metric_tbl}')		##JIN
    return adfc_metric4head, tcs_metric4head, final_metric_tbl, metricdicto

# ------------------------------------------------------------------------------------------------------------------- #
# Functions for Homogenizing and Updating missing Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #

def append_nullheader(headertemp = "HEADER_TEMPELATE_YAML_FILE", filename = "INP_FITS_FILE", 
                      headtbl="HEADER_TABELE_WITH_PARAMS"):
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        headertemp : Name of the YAML input with template headers
        filename   : FITS file whose header has to be appended
        headtbl: The astropy table with correct header values created from the log database
    Returns:
        None
    """
    from datetime import datetime
    now = datetime.now()
    date2prnt = f'{now.year}-{now.month}-{now.day}T{now.strftime("%H:%M:%S")}'
    hdf_temp = yaml.load(open(master_header), Loader = yaml.FullLoader)
    #headtbl = Table.read(headtbl, format="ascii.ipac")

    with fits.open(filename, mode='update', ignore_missing_end=True) as hdulist:
        header = hdulist[0].header
        header['history'] = f"The file is updated by {sys.argv[0]} in batch mode on IST {date2prnt} to add hormonized header information"
        #header.remove('OBJECT', ignore_missing=True, remove_all=True)
        #header.remove('RA')
        #header.remove('DEC')
        #'''
        
        for ii in ['RA','DEC','AIRMASS','SLIT','GRISM']:
            ii = str(ii)
            if ii in header.keys():
                if header[ii]=='':
                    header.remove(ii)
                    print(f' removed vacant keyword {ii}')
                    
                elif ii in ['RA','DEC','AIRMASS','SLIT','GRISM']:
                # removing already existing SLIT ang GRISM informaion to make it in consist format
                    header.remove(ii)
        	        
        	        
        
        #''' 	
        	 
        for keyword in hdf_temp.keys():
        
            if keyword not in header.keys():
                #print (f"DEBUG:: KEYWORD: {keyword}")
                if keyword in headtbl.colnames:
                    header.append(card=(keyword, headtbl[keyword][0], hdf_temp[keyword][1]))
                else:
                    header.append(card=(keyword, hdf_temp[keyword][0], hdf_temp[keyword][1]))
                    
        update_sourcename(filename,header,header['RA'],header['DEC'])
                    

        print("Null Headers updated")

def coorddeg2sr(ra = None, dec = None):
    ra, dec = float(ra), float(dec)
    coord = ICRS(ra*u.degree, dec*u.degree)
    rahmsstr = coord.ra.to_string(u.hour)
    decdmsstr = coord.dec.to_string(u.degree, alwayssign=True)
    rahmsstr = rahmsstr.replace("h",":").replace("m",":").replace("s","")
    decdmsstr = decdmsstr.replace("d",":").replace("m",":").replace("s","")
    return rahmsstr, decdmsstr

### running astroquery (SIMBAD) to find name of the source for given RA and DEC

def name_frm_radec(ra,dec):
    
    ra = float(ra)*u.degree
    dec = float(dec)*u.degree
    print(f' RA = {ra} and DEC = {dec} of the source ')
    coord = SkyCoord(ra ,dec, frame='icrs')
    
    # Query SIMBAD
    result_table = Simbad.query_region(coord, radius='0d0m50s')

    # Check if results are found
    if result_table:
        #for i in range(len(result_table)):
         #     print(f"Object Name: {result_table['MAIN_ID'][i]}")
        kk = result_table['MAIN_ID'][0]
        print(f' name of the source {kk}')
        return result_table['MAIN_ID'][0]
    else:
        print("No objects found at the given coordinates.")

'''
def name_frm_radec(ra,dec):
    
    if isinstance(ra,str):
        try:
            sky_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg)).galactic
            #ra_deg,dec_deg = sky_coord.galactic
            #print(f' the skycoord {sky_coord} and ra {ra_deg} and dec {dec_deg}')
            ra_deg,dec_deg = sky_coord.l.deg*u.degree,sky_coord.b.deg*u.degree
            coord = SkyCoord(ra_deg,dec_deg,frame='galactic')
        except ValueError as e:
            print(f"Error converting RA: {e} in require format")
    else:
        ra_deg = float(ra)*u.degree
        dec_deg = float(dec)*u.degree
        coord = SkyCoord(ra_deg,dec_deg,frame='icrs')
        
    print(f' RA = {ra} and DEC = {dec} of the source ')
    
    
    # Query SIMBAD
    result_table = Simbad.query_region(coord, radius='0d0m2s')

    # Check if results are found
    if result_table:
        #for i in range(len(result_table)):
         #     print(f"Object Name: {result_table['MAIN_ID'][i]}")
        return result_table['MAIN_ID'][0]
    else:
        print("No objects found at the given coordinates.")
'''
'''
def update_sourcename2(header, TYPE, ra=None, dec=None):
    if header['OBJECT NAME'] == '':
        # Print the type of the source in the header for debugging
        print(f' the type of the source in the header: {TYPE}')  ##JIN
        
        # List of predefined source types
        lst = ['BIAS', 'FLAT', 'LAMP', 'LMP', 'ARC']
        
        # Check if TYPE is one of the predefined types in lst
        if TYPE in lst:
            print(f'{TYPE} object')
            header['OBJECT NAME'] = TYPE
            header['RA'] =''
            header['DEC'] = ''
        else:
            # Handle cases where TYPE is not in the predefined list
            if ra is not None and dec is not None:
                name = name_frm_radec(ra, dec)  # Assuming this function is defined elsewhere
                header['OBJECT NAME'] = name
            else:
                print("RA and Dec must be provided for non-standard object types.")
    else:
        print(f"OBJECT NAME mentioned in header: {header['OBJECT NAME']}")
'''

def type_from_Name(filename):
    #print(f' filename in type_from_Name is {filename}')
    filename = os.path.basename(filename)
    filename = filename.upper()
       
    print(f' TYPE for {filename}')
    if "_ARC" in filename:
        TYPE = "ARC"
        CATEGORY = "Calibration"
    elif "FLAT" in filename:
        TYPE = "FLAT"
        CATEGORY = "Calibration"
    elif "BIAS" in filename:
        TYPE = "BIAS"
        CATEGORY = "Calibration"
    elif "LAMP" in filename or "LMP" in filename:
        TYPE = "LAMP"
        CATEGORY = "Calibration"
    elif "LAMP1" in filename or "LMP1" in filename or "HGAR" in filename:
        TYPE = "LAMP1"
        CATEGORY = "Calibration"
    elif "LAMP2" in filename or "LMP2" in filename or "NEON" in filename:
        TYPE = "LAMP2"
        CATEGORY = "Calibration" 
    else:
        TYPE = ''
        CATEGORY = ''
        # TYPE = 'OBJECT'
        # CATEGORY = 'SCIENCE'

    return TYPE, CATEGORY 
            
def update_sourcename(filename,header, ra=None, dec=None):
    #print(f' the type of the source in the header: {TYPE}')  ##JIN
    TYPE,CATGRY = type_from_Name(filename)
    print(f' the type from the filename is {TYPE} and CATEGORY {CATGRY} in update_sourcename function') 
    # List of predefined source types
    lst = ['BIAS', 'FLAT', 'LAMP', 'LMP', 'ARC','HGAR','LAMP1','LAMP2','LMP1','LMP2','HGAR','NEON']
    
    if TYPE in lst:
        print(f'{TYPE} object and the object name in header is {header["OBJECT"]}')
        if header['OBJECT'] == '':
            header['OBJECT'] = TYPE
        else:
            print(' check for the type or object name')
    elif TYPE not in lst:
        # Handle cases where TYPE is not in the predefined list
        if ra is not None and dec is not None:
             name = name_frm_radec(ra, dec)  # Assuming this function is defined elsewhere
             header['OBJECT'] = name
        else:
             print("RA and Dec must be provided for non-standard object types.")

     
        
def GetModeProd_updateOthers(filename, telsiteinfo):

    def GetFileInforFromName(filename):
        #print(f' filename in GetFileInforFromName is {filename}')   ##JIN
        filename = os.path.basename(filename)
        filename = filename.upper()
       
        #print(f' TYPE for {filename}')
        if "_ARC" in filename:
            TYPE = "ARC"
            CATEGORY = "Calibration"
        elif "FLAT" in filename:
            TYPE = "FLAT"
            CATEGORY = "Calibration"
        elif "BIAS" in filename:
            TYPE = "BIAS"
            CATEGORY = "Calibration"
        elif "LAMP" in filename or "LMP" in filename:
            TYPE = "LAMP"
            CATEGORY = "Calibration"
        else:
            TYPE = 'OBJECT'     # JJ edited on 20250710
            CATEGORY = 'SCIENCE'   # JJ edited on 20250710
        return TYPE, CATEGORY

    #Function to decide mode, type category and update others
    type_key, category_key = GetFileInforFromName(filename)
    with fits.open(filename, mode='update', ignore_missing_end=True) as hdulist:
        header = hdulist[0].header
        grism = header['GRISM'].upper()
        if isinstance(header['SLIT'],str):
            slit = header['SLIT'].upper()
        else:
            slit = header['SLIT']
        filter1 = header['FILTER1'].upper()
        filter2 = header['FILTER2'].upper()
        '''
        if (grism is not 'CLEAR' ) and (slit is not 'CLEAR') :
            mode = "Spectroscopy"
        else:
            if (filter1 is not 'CLEAR') and (filter2 is not 'CLEAR') :
                mode = "Imaging" 
            else: 
                print(f' the filter information : grism:{grism} | slit:{slit} | filter1:{filter1} | filter2:{filter2}')  ##JIN
                mode = "NU"   
        '''
        print(f' the filter information : grism:{grism} | slit:{slit} | filter1:{filter1} | filter2:{filter2}')  ##JIN
        
        if (grism != 'CLEAR') and (slit != 'CLEAR'):
        	 mode ='Spectroscopy'
        elif (filter1 != 'CLEAR'):
            mode = 'Imaging'
        elif (filter2 != 'CLEAR'):
        	 mode = 'Imaging'
        elif (grism == 'CLEAR') and (slit =='CLEAR') and (filter1 == 'CLEAR') and (filter2 =='CLEAR'):
            mode = 'NULL'
        else:
        	 mode ='NULL'   ## added this line for those frames (especially test and IVT frames) where either slit is placed or grism or both slit and any one filter   ##JIN
            
        #'''    
        header.set("MODE", mode, "operational mode as per LOG database")
        header.set("TYPE", type_key)
        header.set("CATEGORY", category_key)

        header = UpdateOthers(header = header, telsiteinfo = telsiteinfo)

def UpdateOthers(header = "HEADER_OBJ", telsiteinfo = "sitelinfo" ):

    telsiteinfofile = yaml.load(open(telsiteinfo), Loader=yaml.FullLoader)

    def get_LST(telsiteinfo = telsiteinfofile, dateobs = None):
        """
        Defines a Telescope object in Ephem.
        Args:
           telescopename : Name of the Telescope from which the data was observed
        Returns:
           telescope     : Ephem.Observer object containing Telescope site details
        """
        OBS_LONG = telsiteinfofile['TELESCOPES']['DOT']['LONGITUDE']
        OBS_LAT = telsiteinfofile['TELESCOPES']['DOT']['LATITUDE']
        OBS_ALT = telsiteinfofile['TELESCOPES']['DOT']['ALTITUDE']

        def strlonglat2degree(lat = None, lon = None):
            latv = np.array(lat.split(':'), float)
            latindeg = latv[0] + latv[1]/60 + latv[2]/3600
            lonv = np.array(lon.split(':'), float)
            lonindeg = lonv[0] + lonv[1]/60 + lonv[2]/3600
            return latindeg, lonindeg

        OBS_LAT_deg, OBS_LONG_deg = strlonglat2degree(lat = OBS_LAT, lon = OBS_LONG)
        obslocation = EarthLocation(lat = OBS_LAT_deg*u.degree, lon = OBS_LONG_deg*u.degree, height = OBS_ALT*u.m)
        timenow = Time(dateobs, scale="utc", location = obslocation)
        lst_angle = timenow.sidereal_time('mean')
        lst_str = Angle(lst_angle).to_string(unit = u.hourangle, sep = ":", precision = 0, pad = True)
        return lst_angle, lst_str
    
    #Populating the FRAME with DATE-OBS and transform DATE-OBS to UTC
    ## If the samecode run twice then DATE_OBS will be modified, hence care should be taken. Delete the file and re-run. 
    
    Tist = Time(header['DATE-OBS'], format = 'isot', scale='local') 
    delta = TimeDelta(5.5*u.hour)
    Tutc = Tist - delta 
    header.set('FRAME',Tutc.isot, 'date & time (in UT) of the observations')
    #header['FRAME'] = 'date & time (in UT) of the observations'
    #header.set('DATE-OBS', Tutc.isot)
    #Extracting LST from FRAME and telescope site information
    lst = get_LST(telsiteinfo = telsiteinfo, dateobs = Tutc.isot)    
    header.set("LST", lst[1])
    #add telescope information
    header.set("TELESCOP", "DOT", telsiteinfofile['TELESCOPES']['DOT']['NAME'])
    header.set("INSTRUME", "ADFOSC", "ARIES-Devasthal Faint Object Spectrograph")
    header.set("GEOLONG", telsiteinfofile['TELESCOPES']['DOT']['LONGITUDE'])
    header.set("GEOLAT", telsiteinfofile['TELESCOPES']['DOT']['LATITUDE'])
    header.set("GEOELEV", telsiteinfofile['TELESCOPES']['DOT']['ALTITUDE'])
    #header.set("RA", telsiteinfofile['TELESCOPES']['DOT']['RA'])
    #header.set("DEC", telsiteinfofile['TELESCOPES']['DOT']['DEC'])
    header.set("JD", Tutc.jd, 'observation time in JD')
    header.set("DETECTOR", "ADFOSC-SCI-CAMERA", 'ANDOR Camera as main detector')
    header.set("CCD-TEMP", "-120", 'proxy temperature of the CCD')
    return header

        
def append_templatefiledetails(filename):
    """
    Append Filter details onto the header of the file 'filename'.
    Args:
        filename : FITS file to which filter details have to be appended.
    Returns:
        None
    """
    if filename in header_df.index:
        with fits.open(filename, mode='update') as hdulist:
            header = hdulist[0].header
            columns = header_df.columns
            for column in columns:
                val = header_df.loc[header_df.index == filename, column].values[0]
                if column.upper() in header.keys():
                    header.set(column.upper(), str(val).upper(), '')
    else:
        display_text("ERROR: File {0} not logged in {1}".format(filename, file_template))
        pass


def fix_invalidheaders(filename):
    """
    Filter invalid Header keywords from the FITS header.
    Args:
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    fixed_hdul = fits.HDUList()

    with fits.open(filename) as hdulist:
        for hdu in hdulist:
            try:
                # First try to fix the fixable cards using 'fix' verification
                hdu.verify('fix')
            except:
                # Handle exception thrown due to non-fixable cards
                try:
                    # The error message of 'exception' verification will inform which cards have problem
                    hdu.verify('exception')
                except fits.VerifyError as err:
                    err_str = str(err)
                    # Find all cards from the error message where FITS standard error occured
                    match_iterator = re.finditer(r"Card \'(.+)\' is not FITS standard", err_str)

                    count = 0  # Count items in match_iterator
                    for match in match_iterator:
                        count += 1
                        # Check if this error is due to invalid value string
                        # i.e. possibly non-ASCII string since it couldn't be fixed
                        if match.group() + " (invalid value string:" in err_str:
                            # Fetch card's key & replace its value with empty string
                            hdr_key = match.group(1)
                            hdu.header[hdr_key] = ""
                            print(repr(hdu.header))  # WORK-AROUND: The header has to be printed to save changes
                            print("-" * 60)

                    if count == 0:  # If the match_iterator is empty, there's some unhandlable error
                        raise err

            fixed_hdul.append(hdu)

        # Save the changes to fix fixable cards while saving
        fixed_hdul.writeto("fixed_{0}".format(filename), output_verify='fix', overwrite=True)


def format_header(telescopename, filename):
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
    header['ORIGFILE'] = filename

    print('step111111111111111',telescopename)
    if telescopename == 'ST':
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

    elif telescopename == 'DFOT':
        if 'FRAME' in header.keys():
            header[DATE_keyword] = format_dateobs(header['FRAME'])
        elif 'DATE' in header.keys():
            header[DATE_keyword] = format_dateobs(header['DATE'])
        else:
            display_text("ERROR: DATE keyword not found in the Header")
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


def append_telescopedetails(telescopename, instrument_df, filename):
    """
    Append Observatory details to the header of the file 'filename'
    Args:
        telescopename : Name of the Telescope from which the data was observed
        instrument_df : Pandas DataFrame containing Master-list of Instruments
        filename      : FITS file whose header has to be appended
    Returns:
        None
    """
    _, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, HORIZON, ZENITH = telescope_df.loc[telescopename].values
    instrument_df = instrument_df.loc[telescopename]

    instrument = identify_instrument(instrument_df, filename)
    if instrument in instrument_df['Instrument'].values:
        INSTRUME, PIXSCALE, RNOISE, GAIN, _, _ = instrument_df.loc[instrument_df['Instrument'] == instrument].values[0]

        dict_append = {'TELESCOP': telescopename, 'GEOLONG': OBS_LONG, 'GEOLAT': OBS_LAT, 'GEOELEV': OBS_ALT,
                       'HORIZON': HORIZON, 'ZENITH': ZENITH, 'TIMEZONE': OBS_TIMEZONE, 'INSTRUME': INSTRUME,
                       'PIXSCALE': PIXSCALE, 'DET-READ': RNOISE, 'DET-GAIN': GAIN}

        with fits.open(filename, mode='update') as hdulist:
            header = hdulist[0].header
            for keyword, value in dict_append.items():
                header[keyword] = value
    else:
        pass


def append_astrometryheader(filename):
    """
    Append Astrometry Keywords to the header of the file 'filename'. Also, compute central RA and DEC
    using those keywords and append them to the header
    Args:
        filename : FITS file to which Astrometry header details have to be appended
    Returns:
        None
    """
    new_filename = filename.split('.')[0] + '.new'
    headernew = fits.getheader(new_filename, ext=0)
    w = wcs.WCS(new_filename, naxis=2)
    radec = w.wcs_pix2world(headernew['NAXIS1'] / 2, headernew['NAXIS2'] / 2, 1)

    astrometrykeys = ['WCSAXES', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2',
                      'CUNIT1', 'CUNIT2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header

        for keyword in astrometrykeys:
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header.append(card=(keyword, headernew[keyword]))

        for idx, keyword in enumerate([RA_keyword, DEC_keyword]):
            header.remove(keyword, ignore_missing=True, remove_all=True)
            header. append(card=(keyword, str(radec[idx])))




# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Functions for Performing Astrometry
# ------------------------------------------------------------------------------------------------------------------- #


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

        LST_deg = get_LST(telescopename, date_obs)
        print(LST_deg,date_obs,'LST_deg & date_obs')
        # if telescopename == 'DFOT':
        try:
            subprocess.call('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIX_l + ' --scale-high '
                            + PIX_u + ' --scale-units app --config '+DIR_DOC+'Astrometry.cfg --ra ' + str(LST_deg) + ' --dec ' +
                            str(OBS_LAT) + ' --radius 100 ' + filename, timeout=60, shell=True)
        except subprocess.TimeoutExpired:
            display_text("Astrometry Timed Out For '{0}'".format(filename))
            return False
        else:
            display_text("Astrometry Ran Sucessfully For '{0}'".format(filename))
            os.system('rm -rf *.axy *.corr *.xyls *.match *.rdls *.solved *.wcs')
            return True
        print('Astrometry successfull---------------------------111--') 
        # elif telescopename == 'ST':
        #    os.system('solve-field --continue --downsample 2 --no-plots --scale-low ' + PIXSCALE_l + ' --scale-high '
        #              + PIXSCALE_u + ' --scale-units app --config Astrometry.cfg  ' + filename)
        # else:
        #    sys.exit(1)
    except KeyError:
        display_text("ERROR: Header Keyword 'PIXSCALE' not found")
        return False
    print('Astrometry successfull-----------------------------')

def main():
    #Input fits file which should be used to manipulate 
    parser = argparse.ArgumentParser(description='Handle input fits file and + to update headers')
    parser.add_argument('-i', '--inpfl', dest='inpfile', type=str, nargs='+',
                        help='pass the name of fits file of list to update the headers')
    parser.add_argument('-l', '--log_path', dest='log_path',  default=None, type = str,
                    help='Complete link of the LOG files database')
    parser.add_argument('-m', '--lmode', dest='logmode',  default=None, type = int,
                    help='The interger to supply the mode of handling log files.' 
                         +'\t Options are \'0\' or \'1\' with default of \'0\''
                         + '\n 0: using frame time and a tolarance to list the logs for database'
                         + '\n 1: using all logs in the database and create matrix for further processing')
    parser.add_argument('-t', '--ttol', dest='ttol',  default=None, type = float,
                    help='The tolerance time used to pinpoint the PARAMS to couple with the DATE-OBS in the frame')
    parser.add_argument('--listtol', dest='listtolt',  default=None, type = float,
                    help='The tolerance time used to list all the log files to create database')
    parser.add_argument('-c', '--conf', dest='conf',  default="/home/archive/Documents/program/DOT/ADFOSC/ADFOSC-ARIES-main_V1.2/conf_setting.yaml", type = str,
                    help='The YAML file containing the information about the site and telescopes')
    parser.add_argument('--sitelinfo', dest='sitelinfo',  default=None, type = str,
                    help='The YAML file containing the information about the site and telescopes')
    parser.add_argument('--headtemp', dest='headtemp',  default=None, type = str,
                    help="The YAML file containing the header template as per ARIES's convension")
    
    args = parser.parse_args()
    inpf = args.inpfile
    log_path = args.log_path
    logmode = args.logmode
    listfltr_ttol = args.listtolt
    obstime_ttol = args.ttol
    sitelinfo = args.sitelinfo
    headtemp = args.headtemp

    #managing input parameters through a configuration file
    conffile = yaml.load(open(args.conf), Loader = yaml.FullLoader)
    if log_path == None: 
        log_path = conffile['log_path']
        
        
    if not log_path:
        log_path = input("Please provide the full path of the LOG FILE database")
    if not logmode:
        logmode = conffile['log_access_mode']
    if not listfltr_ttol:  
        listfltr_ttol = conffile['list_tol_time']
    if not obstime_ttol:
        obstime_ttol = conffile['dt_access_time'] 
    if not sitelinfo:
        sitelinfo = conffile['SiteTelInfoFile']
    if not headtemp:
        headtemp = conffile['Tempelate_Headers']          

    txtInp4prnt =  "-----------\n"
    txtInp4prnt += f"INPUT FILE : {inpf}\n"
    txtInp4prnt += f"Configulation File : {args.conf}\n"
    txtInp4prnt += f"LOG PATH : {log_path}\n"
    txtInp4prnt += f"YAML FILE WITH SITE & TEL INFO : {sitelinfo}\n"
    txtInp4prnt += f"YAML FILE WITH DEFAULT HEADERS : {headtemp}\n"
    txtInp4prnt += f"LOG MODE : {logmode}\n"
    txtInp4prnt += f"LIST TOL : {listfltr_ttol}\n"
    txtInp4prnt += f"OBS TIME TOL : {obstime_ttol}\n"
    txtInp4prnt += "-----------\n"
    print (txtInp4prnt)
    
    
    if type(inpf) == list:
        cnt=0
        for infile in inpf:
            print (f"FILE : {infile}")
            cnt += 1   
            adfocsmetric, tcsmetric, finaltbl, metricdicto = run_proc(infile = infile, log_path=log_path, 
                                               make_list_mode = logmode, listfltr_ttol = listfltr_ttol, 
                                               obstime_ttol = obstime_ttol) 
            print (f"Final Table: {finaltbl.colnames}")
            print(metricdicto)
            #finaltbl.write("TRAUSNKSKSSKSK.ipac", format="ascii.ipac", overwrite=True)
            #sys.exit()
            append_nullheader(headtemp, infile, headtbl = finaltbl)
            GetModeProd_updateOthers(infile, sitelinfo)

            #print (f"ADFOCSMETRIC :\n {adfocsmetric} \n LENGTH: {len(adfocsmetric[0])}")
            #print (f"TCSMETRIC :\n {tcsmetric} \n LENGTH: {len(tcsmetric[0])}")
            

    else:
        print (f"INP FILE TYPE : {type(inpf)}")
        adfocsmetric, tcsmetric, finaltbl = run_proc(infile = infile, log_path=log_path, 
                                           make_list_mode = logmode, listfltr_ttol = listfltr_ttol, 
                                           obstime_ttol = obstime_ttol)
        #print (f"ADFOCSMETRIC :\n {adfocsmetric} \n LENGTH: {len(adfocsmetric)}")
        #print (f"TCSMETRIC :\n {tcsmetric} \t {len(tcsmetric)}")


if __name__ == "__main__":
    import sys
    sys.settrace
    scriptname = sys.argv[0]
    print (f"Running script {scriptname}")
    sys.setrecursionlimit(500000)
    main()    
