# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# python ./utils/download_sciencebase_data.py

# Standard library imports
import os
import json

# External dependencies imports
import requests
from sciencebasepy import SbSession

# -------------------------------------------------- Constants --------------------------------------------------
sb = SbSession()
outputted_json_name = "sciencebase_id_to_title.json"

# -------------------------------------------------- Global Variables --------------------------------------------------
id_to_title = {}

# -------------------------------------------------- Helper Methods --------------------------------------------------
def download_children(item_id: str, parent_dir_path: str) -> None:
    """
    Recursively downloads all children of the item with the given ID.

    Args:
        item_id (str): ID of the item that potentially contains child items to download
        parent_dir_path (str): Location of the directory where child items are downloaded to
    """
    global sb, id_to_title
    item_json = sb.get_item(item_id)
    item_title = item_json["title"]
    # Create a directory to store the given item's own children/attached items.
    id_to_title[item_id] = item_title
    item_dir_path = os.path.join(parent_dir_path, item_id)
    if not os.path.exists(item_dir_path): os.makedirs(item_dir_path)
    if item_json["hasChildren"]:
        # print("Looking in {} for data files...".format(item_title))
        # Look through each child's own child items for raw data files.
        children_ids = sb.get_child_ids(item_id)
        for child_id in children_ids:
            download_children(item_id = child_id, parent_dir_path = item_dir_path)
    else:
        # print("Downloading data from {}...".format(item_title))
        # Download raw data files that are attached to the current item.
        sb.get_item_files(item_json, destination = item_dir_path, progress_bar = True)

# -------------------------------------------------- Main Program --------------------------------------------------
if __name__ == "__main__":
    # 1. Ask the user for the item ID of the data that they want to download.
    root_item_id = input("Please enter the data's item ID: ")
    # 2. Check if the item ID exists in ScienceBase.
    print("Checking if {} is a valid item ID...".format(root_item_id))
    item_exists = requests.Session().get("https://www.sciencebase.gov/catalog/item/" + root_item_id)
    # 3. Show the download progress.
    parent_data_dir_path = os.path.abspath("./utils")
    root_item_dir_path = os.path.join(parent_data_dir_path, root_item_id)
    if item_exists:
        # Download raw data files for the inputted item ID.
        print("Inputted item ID exists in ScienceBase: Starting download...")
        download_children(item_id = root_item_id, parent_dir_path = parent_data_dir_path)
        # Save dictionary which maps each item's ID to their title as a JSON.
        # ^ File paths with longer than 256 characters causes FileNotFoundErrors in Windows (https://github.com/python/cpython/issues/89935).
        with open(os.path.join(root_item_dir_path, outputted_json_name), "w") as json_file:
            json.dump(id_to_title, json_file, indent = 4)
        print("Download complete! All data files are saved in {}.".format(root_item_dir_path))
    else:
        # Display an error message if the item ID doesn't exist in ScienceBase.
        print("Error: Inputted item ID does not exist in ScienceBase. Please run this script again with an existing item ID.")