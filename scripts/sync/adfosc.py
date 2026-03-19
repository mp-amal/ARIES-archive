#!/usr/bin/env python3

import os
import sys
import subprocess
import yaml
from datetime import datetime
CONFIG_FILE = "/home/archive/Documents/ARIES-archive/config/credentials/NAS_MANORA.yaml"

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
        return f"{year}A"
    elif month in [10, 11, 12]:
        return f"{year}B"
    elif month == 1:
        return f"{year - 1}B"
    else:
        return None
    
cycle =get_current_cycle_name()


def main():
    try:
        with open(CONFIG_FILE, "r") as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: {CONFIG_FILE} not found.")
        sys.exit(1)

    try:
        source = config["adfosc_source"]
        destination = os.path.join(config["adfosc_destination"])
        password = config["password"]
    except KeyError as e:
        print(f"Missing key in YAML file: {e}")
        sys.exit(1)

    cmd = [
        "rsync",
        "-av",
        source,
        destination,
    ]

    env = os.environ.copy()
    env["RSYNC_PASSWORD"] = password

    print("Running:", " ".join(cmd))

    try:
        result = subprocess.run(cmd, env=env, check=True)
        print("\nSync completed successfully.")
        sys.exit(result.returncode)
    except subprocess.CalledProcessError as e:
        print(f"\nrsync failed with return code {e.returncode}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        print("\nError: rsync is not installed or not found in PATH.")
        sys.exit(1)

if __name__ == "__main__":
    main()