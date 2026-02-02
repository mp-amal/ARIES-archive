# ARIES-archive
This repository includes automation scripts for data archival from multiple ARIES-operated telescopes: DFOT, ST, and DOT (ADFOSC and TANSPEC instruments). Each telescope and instrument has its own processing requirements, which are handled through separate, dedicated scripts. All telescope- and instrument-specific pipelines are coordinated using a central schedule script.

# 1.3m DFOT (Devasthal Fast Optical Telescope)

The DFOT script automatically fetches science FITS data from the observatory PC. It performs preprocessing using standard FITS header keywords, applies precise astrometric calibration, and stores the processed data in the database with unique file names.

DFOT kinetic mode observations require special handling due to their data format. The script includes telescope-specific logic to correctly process these observations. All DFOT data processing is fully automated and does not require human intervention. If any error occurs during processing, the details are recorded in a JSON file for reporting.

# 1.04m ST (Sampurnanand Telescope)

For ST observations, the script handles both imaging and polarimetric data. Since polarimetric observations differ from standard image frames, the pipeline processes them using instrument-specific routines to ensure correct data organization and archival.

# 3.6m DOT (Devasthal Optical Telescope)

DOT observations are instrument-dependent and more complex. Data may come from different instruments mounted on the axial port, primarily ADFOSC and TANSPEC. The script first identifies which instrument was active during the observational night and then selects the corresponding instrument-specific pipeline for processing and archival.

# Control and Automation

All telescope- and instrument-specific scripts are controlled by a schedule script, which manages the execution order and ensures that the correct pipeline is applied for each observational night.

# Notification System

The MAIL.py script manages the notification system. It sends daily status updates on telescope observations, including processing success or failure and archival completion. Error information collected in JSON files is included in these notifications to support monitoring and troubleshooting.


All scripts generate log files that record processing status, warnings, and errors. These logs help verify whether all expected observations have been archived successfully.
