import os
import json

folder_path = "/data/archived_data/astro_data/final_data/"  # change this
subfolders = [
    name for name in os.listdir(folder_path)
    if os.path.isdir(os.path.join(folder_path, name))
]

telescopes = [ 'DFOT','DOT','ST']
json_path = "/home/archive/Documents/ARIES-archive/logs"

def has_6digit_component(path):
    parts = path.split(os.sep)
    for part in parts:
        if len(parts) ==10:
            return path, parts
        else:
            # pass
            return None, parts    

def update_counts(date_str, instrument, total, science, flat, bias, lamp):
    # Ensure date key exists
    if date_str not in data:
        data[date_str] = {}

    # Add / overwrite instrument block
    data[date_str][instrument] = {
        "total": total,
        "science": science,
        "flat": flat,
        "bias": bias,
        "lamp": lamp
    }
i=0
for telescope in telescopes:
    # print(telescope)
    i=0
    json_file = os.path.join(json_path,telescope+'_json.json')
    if not os.path.exists(json_file):
        with open(json_file, "w") as f:
            json.dump({}, f, indent=4)   # empty JSON object
    if os.path.exists(json_file):
        with open(json_file, "r") as f:
            data = json.load(f)
    for year in subfolders:
        for path,dirs,files in os.walk(os.path.join(folder_path,year)):
            if telescope in path:
                # print(path)
            
                last_part = os.path.basename(path)

                # Check if last part is exactly an 8-digit date (YYYYMMDD)
                if last_part.isdigit() and len(last_part) == 8:
                    # print(path)
                    parts = path.split(os.sep)
                    # print(parts)

                    fits_files = [
                        f for f in os.listdir(path)
                        if os.path.isfile(os.path.join(path, f))
                        and f.lower().endswith((".fits", ".fit"))
                    ]
                    s_count = sum(1 for f in fits_files if f.startswith("S"))
                    f_count = sum(1 for f in fits_files if f.startswith("F"))
                    b_count = sum(1 for f in fits_files if f.startswith("B"))
                    l_count = sum(1 for f in fits_files if f.startswith("L"))
                    # print(parts[7],parts[8],parts[9])
                    # print("S* files:", s_count)
                    # print("F* files:", f_count)
                    # print("B* files:", b_count)
                    # print("L* files:", l_count)

                    total = s_count + f_count + b_count
                    # print("Total S/F/B files:", total)
                    update_counts(parts[9],parts[8], total, s_count, f_count, b_count, l_count)
                    # if os.path.exists(json_file):
                    with open(json_file, "w") as f:
                        json.dump(data, f, indent=4)



 
# json_file = os.path.join(json_path,'DOT_json.json')
# if not os.path.exists(json_file):
#     with open(json_file, "w") as f:
#         json.dump({}, f, indent=4)   # empty JSON object
# if os.path.exists(json_file):
#     with open(json_file, "r") as f:
#         data = json.load(f)

# for path,dirs,files in os.walk('/data/DOT/ADFOSC/Processed_Data/2026-C1'):
#     # print(path)


#     last_part = os.path.basename(path)

#     # Check if last part is exactly an 8-digit date (YYYYMMDD)
#     if last_part.isdigit() and len(last_part) == 8:
#         # print(path)
#         parts = path.split(os.sep)
#         # print(parts)

#         fits_files = [
#             f for f in os.listdir(path)
#             if os.path.isfile(os.path.join(path, f))
#             and f.lower().endswith((".fits", ".fit"))
#         ]
#         s_count = sum(1 for f in fits_files if f.startswith("S"))
#         f_count = sum(1 for f in fits_files if f.startswith("F"))
#         b_count = sum(1 for f in fits_files if f.startswith("B"))
#         l_count = sum(1 for f in fits_files if f.startswith("L"))
#         print(parts[6],parts[5],parts[3])
#         print("S* files:", s_count)
#         print("F* files:", f_count)
#         print("B* files:", b_count)
#         print("L* files:", l_count)

#         total = s_count + f_count + b_count
#         # print("Total S/F/B files:", total)
#         update_counts(parts[6],parts[3], total, s_count, f_count, b_count, l_count)
#         # if os.path.exists(json_file):
#         with open(json_file, "w") as f:
#             json.dump(data, f, indent=4)


# # ST
# s_count = 0
# f_count = 0
# b_count = 0
# l_count = 0
# # print(parts[7],parts[8],parts[9])
# # print("S* files:", s_count)
# # print("F* files:", f_count)
# # print("B* files:", b_count)
# # print("L* files:", l_count)

# total = 0

# json_file = os.path.join(json_path,'ST_json.json')
# if not os.path.exists(json_file):
#     with open(json_file, "w") as f:
#         json.dump({}, f, indent=4)   # empty JSON object
# if os.path.exists(json_file):
#     with open(json_file, "r") as f:
#         data = json.load(f)
# for year in subfolders:
#     for path,dirs,files in os.walk(os.path.join(folder_path,year)):
        
#         if 'ST' in path:
#             print(path)
#             # print(path)
#             last_part = os.path.basename(path)
#             # Check if last part is exactly an 8-digit date (YYYYMMDD)
#             if last_part.isdigit() and len(last_part) == 8:
#                 # print(path)
#                 parts = path.split(os.sep)
#                 # print(parts)

#                 fits_files = [
#                     f for f in os.listdir(path)
#                     if os.path.isfile(os.path.join(path, f))
#                     and f.lower().endswith((".fits", ".fit"))
#                 ]
#                 s_count = sum(1 for f in fits_files if f.startswith("S"))
#                 f_count = sum(1 for f in fits_files if f.startswith("F"))
#                 b_count = sum(1 for f in fits_files if f.startswith("B"))
#                 l_count = sum(1 for f in fits_files if f.startswith("L"))
#                 # print(parts[7],parts[8],parts[9])
#                 # print("S* files:", s_count)
#                 # print("F* files:", f_count)
#                 # print("B* files:", b_count)
#                 # print("L* files:", l_count)

#                 total = s_count + f_count + b_count
#                 # print("Total S/F/B files:", total)
#                 update_counts(parts[9],parts[8], total, s_count, f_count, b_count, l_count)
#                 # if os.path.exists(json_file):
#                 with open(json_file, "w") as f:
#                     json.dump(data, f, indent=4)