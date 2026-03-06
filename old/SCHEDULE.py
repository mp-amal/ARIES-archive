import os
import time
import schedule
import json
from datetime import datetime,date
import json
import schedule
# import aplpy
import paramiko
from scp import SCPClient
from datetime import datetime, timedelta

yesterday = datetime.now() - timedelta(days=1)
yesterday =  yesterday.strftime("%Y%m%d")
# print(yesterday)

#yesterday = 'DOT36_Images_for_poster'
# path links
remote_path = os.path.join('./2025-C2',yesterday)
remote2_path = os.path.join('./2025-C2',yesterday)
# rawpath = os.path.join('/data/DOT/ADFOSC/rawdata',yesterday)
# Processing_path = os.path.join('/data/DOT/ADFOSC/ProcessedData/2024-C2',yesterday)

def check_remote_folder(ssh, path):
    """Check if a remote folder exists."""
    stdin, stdout, stderr = ssh.exec_command(f"if [ -d '{path}' ]; then echo 'Exists'; else echo 'NotFound'; fi")
    return stdout.read().decode().strip()

def check_folder_from_server(remote_path, remote2_path):
    try:
        server_ip = "userip"  # Replace with your server's IP
        username = "username"    # Replace with your username
        password = "password"    # Replace with your password

        # Establish SSH connection
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(hostname=server_ip, username=username, password=password)
        # Check first remote path
        folder_status = check_remote_folder(ssh, remote_path)
        if folder_status == "Exists":
            instrument = "ADFOSC"
            # print('ADFOSC instrument mounted')

            # print(f"Folder not found at {remote_path}, checking {remote2_path}...")
        folder_status = check_remote_folder(ssh, remote2_path)
        if folder_status == "Exists":
            remote_path = remote2_path  # Use the second path
            instrument = 'TANSPEC'
        ssh.close()
        print(instrument)
        return instrument
    except paramiko.SSHException as e:
        print(f"SSH error: {e}")
    except FileNotFoundError as e:
        print(f"File error: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

def schedule_fn():
    
    telescopes = ['DFOT', 'ST','DOT']
    for telescope in telescopes:
        start_time = time.time()
        print(telescope+'\n')

        if telescope == 'DFOT':
          print('././scripts/DFOT.py')
          os.system(f'python3 ./scripts/DFOT.py')
        elif telescope == 'ST':
          print('\n././scripts/ST.py')
          os.system(f'python3 ./scripts/ST.py')
        elif telescope == 'DOT':
            if instrument == 'TANSPEC':
                print('\n././scripts/DFOT.py')
                os.system(f'python3 ./TANSPEC.py')
            if instrument == 'ADFOSC':
                print('\n././scripts/DFOT.py')
                os.system(f'python3 ./scripts/ADFOSC.py')
            else:
                print('No Observation')
        # main()
        end_time =time.time()
        execution_time = end_time - start_time
        print(f"Execution time: {execution_time} seconds")
    os.system(f'python3 ./scripts/json_update.py')
    os.system(f'python3 ./scripts/MAIL.py')


    
instrument = check_folder_from_server(remote_path, remote2_path)
# print(instrument)


# # schedule_fn()
# json_path =  '/home/archive/Documents/ADA_PROGRAM/output/JSON/daily_json.json'

# today_date = datetime.today().strftime('%Y%m%d')
# print('\n',today_date,'\n')

schedule_fn()






new_data = {
    "DFOT": {
        "status": None,
        "instrument": None,
        "file count": None,
        "manual-check": None
    },
    "ST": {
        "status": None,
        "instrument": None,
        "file count": None,
        "polarimetry": None,
        "manual-check": None
    },
    "DOT": {
        "status": None,
        "file count": '-',
        "instrument": instrument,
        "manual-check": 'Yes'
    }
}
