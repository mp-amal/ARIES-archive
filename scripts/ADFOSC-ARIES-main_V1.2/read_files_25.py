import glob
import os
import subprocess
import sys

script_dir = os.getcwd()
########################################################################################################################################
#  This script is to read files from folder of interest to update the entire (if specified) fits file.
#  Read the files from /observation_Data/ADFOSC/{cycle}/{folder} and save the updated fits file in /home/adfosctest/HEADER_updatation/{cycle}/{folder}.
# calls home/adfosctest/HEADER_updatation/ADFOSC-ARIES-main_2024/2024_Head_Manage_ULogs_v01.py which updates the FITS files.
########################################################################################################################################

####################################################
def save_unique_filenames(input_file, output_file):
    try:
        # Read the log file
        with open(input_file, 'r') as f:
            file_lines = f.readlines()

        # Remove duplicates by converting the list to a set, then back to a list
        unique_files = list(set(file_lines))

        # Sort the unique filenames (optional)
        unique_files.sort()

        # Write the unique filenames to the new output file
        with open(output_file, 'w') as f:
            f.writelines(unique_files)

        print(f"Unique filenames saved to {output_file}")

    except Exception as e:
        print(f"Error processing the file {input_file}: {e}")
        
    os.remove(input_file)


########################################################


main_dir = '/home/adfosctest/HEADER_updatation/2024-C1/'
log_header_dir = os.path.join(main_dir,"LOG_header_update")
os.makedirs(log_header_dir,exist_ok=True)
f = glob.glob(os.path.join(main_dir+'2024_03_0*'))

for fdate in f:
    foldername = os.path.basename(fdate)
    log_file = os.path.join(log_header_dir,f'headerupdated_files_{foldername}+_03012025.txt')
    err_log_file = os.path.join(log_header_dir,f'error_headerupdated_files_{foldername}+_03012025.txt')
    output_log_file = os.path.join(log_header_dir, f'{foldername}_output_log.txt')
    error_log_file = os.path.join(log_header_dir, f'{foldername}_error_log.txt')
    main_folders = os.path.join(main_dir, foldername)
    files_list = []  # Reset file list for each folder
    
    print(f'Processing folder {main_folders}')
    
    # Collect all .fit and .fits files
    for root, dirs, files in os.walk(main_folders):
        for file in files:
            if file.endswith((".fit", ".fits")):
                file_path = os.path.join(root, file)
                #print(f' the file path :{file_path}')
                files_list.append(file_path)

    print(f'The length of the file list is {len(files_list)}')

    # Process each file
    for file_path in files_list:
        cmd = ['python3.8', '2024_Head_Manage_ULogs_v01.py', '-i', file_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        print(f'Processing file: {file_path}')
        print(result.stdout)
        
        if result.returncode == 0:
            # Log the successfully processed file
            with open(log_file, 'a') as log:
                log.write(f'{file_path}\n')
        else:
            # Log errors
            with open(err_log_file, 'a') as elog:
                elog.write(f'{file_path}\n')
            print(f"Error occurred with file: {file_path}")
            print(result.stderr)

        # Log stdout and stderr for the current file
        if result.stdout:
            with open(output_log_file, 'a') as log:
                log.write(f'Output from {file_path}:\n{result.stdout}\n')
        if result.stderr:
            with open(error_log_file, 'a') as elog:
                elog.write(f'Error from {file_path}:\n{result.stderr}\n')
       
    save_unique_filenames(os.path.join(script_dir,"used_telnet.txt"),os.path.join(main_dir,foldername,"used_TCStelnet.txt"))
    save_unique_filenames(os.path.join(script_dir,"unused_telnet.txt"),os.path.join(main_dir,foldername,"unused_TCStelnet.txt"))
    save_unique_filenames(os.path.join(script_dir,"used_adfosc.txt"),os.path.join(main_dir,foldername,"used_ADics.txt"))
    save_unique_filenames(os.path.join(script_dir,"unused_adfosc.txt"),os.path.join(main_dir,foldername,"unused_ADics.txt"))
    
    sys.stdout.flush()
