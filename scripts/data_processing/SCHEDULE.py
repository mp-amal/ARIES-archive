#!/usr/bin/env python3

import subprocess
import sys
from pathlib import Path

scripts = [
    "/home/archive/Documents/ARIES-archive/scripts/sync/dfot.py",
    "/home/archive/Documents/ARIES-archive/scripts/sync/adfosc.py",
    "/home/archive/Documents/ARIES-archive/scripts/data_processing/DFOT.py",
    "/home/archive/Documents/ARIES-archive/scripts/data_processing/ADFOSC.py",
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

print("Both scripts finished successfully.")