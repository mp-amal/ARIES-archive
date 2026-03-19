# ADFOSC-ARIES
A patch to update the headers of ADFOSC observations

Head_Manage_ULogs_v01.py   -- Main script to be used from command line --- it can be easily integrated with a wrapper loop or cron job
It runs from command line with the following options

Running script Head_Manage_ULogs_v01.py
usage: Head_Manage_ULogs_v01.py [-h] [-i INPFILE [INPFILE ...]] [-l LOG_PATH]
                                [-m LOGMODE] [-t TTOL] [--listtol LISTTOLT]
                                [-c CONF] [--sitelinfo SITELINFO]
                                [--headtemp HEADTEMP]

Handle input fits file and + to update headers

options:
  -h, --help            show this help message and exit
  -i INPFILE [INPFILE ...], --inpfl INPFILE [INPFILE ...]
                        pass the name of fits file of list to update the
                        headers
  -l LOG_PATH, --log_path LOG_PATH
                        Complete link of the LOG files database
  -m LOGMODE, --lmode LOGMODE
                        The interger to supply the mode of handling log files.
                        Options are '0' or '1' with default of '0' 0: using
                        frame time and a tolarance to list the logs for
                        database 1: using all logs in the database and create
                        matrix for further processing
  -t TTOL, --ttol TTOL  The tolerance time used to pinpoint the PARAMS to
                        couple with the DATE-OBS in the frame
  --listtol LISTTOLT    The tolerance time used to list all the log files to
                        create database
  -c CONF, --conf CONF  The YAML file containing the information about the
                        site and telescopes
  --sitelinfo SITELINFO
                        The YAML file containing the information about the
                        site and telescopes
  --headtemp HEADTEMP   The YAML file containing the header template as per
                        ARIES's convension
#### Note that the conf file provided shall have many enrties which can also be parsed through command line. The enrties through command line will have precedence over the ones parsed through conf file. It also requires two more YAML files "sitelinfo" and "headtemp" which communicates some of the important general information about site, telescope and instruments at ARIES.   
