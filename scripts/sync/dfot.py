#!/usr/bin/env python3

import os
import sys
import subprocess
import yaml

CONFIG_FILE = "/home/archive/Documents/ARIES-archive/config/credentials/NAS_MANORA.yaml"

def main():
    try:
        with open(CONFIG_FILE, "r") as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: {CONFIG_FILE} not found.")
        sys.exit(1)

    try:
        source = config["rsync"]["source"]
        destination = config["rsync"]["destination"]
        password = config["rsync"]["password"]
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