# cd C:\Users\Venuxk\Projects\data-visualizer\utils
# conda activate visualizer
# python download_sciencebase_data.py

import os
import requests
from sciencebasepy import SbSession
import re

def download_children(item_id, dir_path):
    item_json = sb.get_item(item_id)
    item_title = item_json["title"]
    print(item_json["hasChildren"])
    # subdir_path = re.sub("[ -,]", "_", os.path.join(dir_path, item_title))
    subdir_path = os.path.join(dir_path, item_title)
    # # Create a directory to store the given item's own children items.
    # print("dir", subdir_path, os.path.isdir(subdir_path))
    # if not os.path.exists(subdir_path):
    #     os.makedirs(subdir_path, exist_ok = True)
    #     print("making dir", subdir_path)
    if item_json["hasChildren"]:
        print("Looking in {} for data files...".format(item_title))
        # Also look through each child's own child items for raw data files.
        children_ids = sb.get_child_ids(item_id)
        for child_id in children_ids:
            download_children(child_id, subdir_path)
    else:
        print("Downloading data from {}...".format(item_title))
        # Create a directory to store the given item's own children items.
        print("dir", subdir_path, os.path.isdir(subdir_path))
        if not os.path.exists(subdir_path):
            os.makedirs(subdir_path, exist_ok = True)
            print("making dir", subdir_path)
        # # Get the path to the subdirectory, which will contain all the data files.
        # [data_subdir_path, _] = os.path.split(dir_path)
        # print(os.path.isdir(data_subdir_path), os.path.exists(dir_path))
        # Download raw data files that are attached to the current item.
        sb.get_item_files(item_json, destination = subdir_path, progress_bar = True)

sb = SbSession()
root_data_dir = os.path.abspath(".")
# 1. Ask the user for the item ID of the data that they want to download.
item_id = input("Enter the data's item ID: ")
# 2. Check if the item ID exists in ScienceBase.
print("Checking if {} is a valid item ID...".format(item_id))
item_exists = requests.Session().get("https://www.sciencebase.gov/catalog/item/" + item_id)
# 3. Show the download progress.
if item_exists:
    # Download raw data files for the inputted item ID.
    print("Inputted item ID exists in ScienceBase: Starting download...")
    download_children(item_id, root_data_dir)
    print("Download complete! All data files are saved at {}.".format(root_data_dir))
else:
    # Display an error message if the item ID doesn't exist in ScienceBase.
    print("Error: Inputted item ID doesn't exist in ScienceBase. Please run this script again with an existing item ID.")