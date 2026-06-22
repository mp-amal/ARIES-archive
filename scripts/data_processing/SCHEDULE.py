#!/usr/bin/env python3

import subprocess
import sys
import yaml
# import  pathlib
from pathlib import Path

import sys
from datetime import datetime

def check_cycle():
    date_to_check = datetime.now().date()

    # Valid period: October 1 through June 15, inclusive
    if date_to_check.month >= 10 or (
        date_to_check.month < 6
        or (date_to_check.month == 6 and date_to_check.day <= 15)
    ):
        pass
    else:
        print(f"{date_to_check}")
        sys.exit("Cycle is Over.")

check_cycle()


scripts = [
    "/home/archive/Documents/ARIES-archive/scripts/sync/dfot.py",      # DFOT data sync
    "/home/archive/Documents/ARIES-archive/scripts/sync/adfosc.py",    
    "/home/archive/Documents/ARIES-archive/scripts/sync/st.py",
    "/home/archive/Documents/ARIES-archive/scripts/sync/tanspec.py",

    "/home/archive/Documents/ARIES-archive/scripts/data_processing/DFOT.py",
    "/home/archive/Documents/ARIES-archive/scripts/data_processing/ST.py",
    "/home/archive/Documents/ARIES-archive/scripts/data_processing/ADFOSC.py",
    "/home/archive/Documents/ARIES-archive/scripts/data_processing/TANSPEC.py",
    "/home/archive/Documents/ARIES-archive/scripts/Notifications_syatem/JSON_UPDATE.py",
    "/home/archive/Documents/ARIES-archive/scripts/Notifications_syatem/MAIL.py"

]

for script in scripts:
    path = Path(script)

    if not path.exists():
        print(f"File not found: {path}")
        sys.exit(1)

    print(f"Running: {path} ")
    result = subprocess.run([sys.executable, str(path)])

    if result.returncode != 0:
        print(f"Stopped because {path} failed with exit code {result.returncode}")
        sys.exit(result.returncode)

print("All scripts finished successfully.")



with open("/home/archive/Documents/ARIES-archive/config/credentials/1_9_archive.yaml", "r") as cred:
    arch_cred = yaml.safe_load(cred)



import paramiko
from scp import SCPClient
server_ip = arch_cred['ip']
username = arch_cred['user']
password = arch_cred['password']

# remote_folder = "/home/archive/data/Data_Share/DOT/TANSPEC"
# local_folder = os.path.join(final_path,date)

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
    # ssh.exec_command(f'mkdir -p "{remote_folder}"')
    ssh.exec_command('python3 "/home/archive/Documents/fileshare.py"')
    print()
    # # copy folder
    # with SCPClient(ssh.get_transport()) as scp:
    #     scp.put(local_folder, remote_path=remote_folder, recursive=True)

    # print("Folder copied successfully.")

    ssh.close()

except Exception as e:
    print("Error:", e)
