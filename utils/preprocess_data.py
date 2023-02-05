# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# python ./utils/preprocess_data.py

# Standard library imports
import os
import json
import xml.etree.ElementTree as ET

# External dependencies imports
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import rioxarray as rxr
from download_sciencebase_data import outputted_json_name as sb_download_output_json_name

# -------------------------------------------------- Constants (should match the constants used in DataMap.py) --------------------------------------------------
outputted_dataset_json_name = "dataset_info.json"
dataset_crs_property = "crs"
data_from_download_script_property = "from_download_sciencebase_data_script"
outputted_buffer_json_name = "buffer_config.json"

transects_subdir_name = "Transects"
transect_geojson_id_property = "Transect ID"
transect_geojson_start_point_property = "Start Point"
transect_geojson_end_point_property = "End Point"

# -------------------------------------------------- Global Variables --------------------------------------------------
dataset_crs = None
transects_dir_exists = False
buffer_config = {}

# -------------------------------------------------- Helper Methods --------------------------------------------------
def get_crs_from_xml_file(file_path: str) -> ccrs:
    """
    Read a data file's XML (extensible markup language) to get its CRS (coordinate reference system).

    Args:
        file_path (str): Path to a data file, which is in the same directory as its XML file
            ^ XML file contains information about the data's CRS
    """
    global dataset_crs
    if dataset_crs is None:
        # Get the path to the XML file.
        file_dir_path, _ = os.path.split(file_path)
        xml_file, *_ = [file for file in os.listdir(file_dir_path) if file.endswith(".xml")]
        xml_path = os.path.join(file_dir_path, xml_file)
        # Get proj4 parameters that correspond to metadata provided in the XML file.
        proj4_params = {"no_defs": True, "type": "crs"}
        tree = ET.parse(xml_path)
        root = tree.getroot()
        for projection_element in root.iter("mapprojn"):
            proj_name = projection_element.text
            if proj_name in ["Lambert Conformal Conic"]: proj4_params["proj"] = "lcc"
            else: print("Error getting the CRS from {}: Script did not account for the {} projection yet.".format(xml_path, proj_name))
        for grid_coordinate_system_element in root.iter("gridsysn"):
            proj_name = grid_coordinate_system_element.text
            if proj_name in ["Universal Transverse Mercator"]:
                proj4_params["proj"] = "utm"
                for utm_zone_element in root.iter("utmzone"): proj4_params["zone"] = utm_zone_element.text
            else:
                print("Error getting the CRS from {}: Script did not account for the {} projection yet.".format(xml_path, proj_name))
        for i, standard_parallel_element in enumerate(root.iter("stdparll")):
            if i < 2:
                standard_parallel_param = "lat_{}".format(i + 1)
                proj4_params[standard_parallel_param] = standard_parallel_element.text
            else:
                print("Error getting the CRS from {}: More than 2 standard parallels were provided in the XML file.".format(xml_path))
        for scaling_factor_at_central_meridian_element in root.iter("sfctrmer"): proj4_params["k_0"] = scaling_factor_at_central_meridian_element.text
        for longitude_of_central_meridian_element in root.iter("longcm"): proj4_params["lon_0"] = longitude_of_central_meridian_element.text
        for latitude_of_projection_origin_element in root.iter("latprjo"): proj4_params["lat_0"] = latitude_of_projection_origin_element.text
        for false_easting_element in root.iter("feast"): proj4_params["x_0"] = false_easting_element.text
        for false_northing_element in root.iter("fnorth"): proj4_params["y_0"] = false_northing_element.text
        for planar_distance_units_element in root.iter("plandu"):
            units = planar_distance_units_element.text.lower()
            if units in ["meter", "meters"]: proj4_params["units"] = "m"
            else: print("Error getting the CRS from {}: Script did not account for the {} unit yet.".format(xml_path, planar_distance_units_element.text))
        for horizontal_datum_element in root.iter("horizdn"):
            datum_name = horizontal_datum_element.text
            if datum_name in ["NAD83 (CORS96)", "D_North_American_1983"]: proj4_params["datum"] = "NAD83"
            elif datum_name in ["WGS1984"]: proj4_params["datum"], proj4_params["proj"] = "WGS84", "longlat"
            else: print("Error getting the CRS from {}: Script did not account for the {} datum yet.".format(xml_path, datum_name))
        for ellipsoid_element in root.iter("ellips"):
            ellips_name = ellipsoid_element.text
            if ellips_name in ["GRS_1980"]: proj4_params["ellps"] = "GRS80"
            elif ellips_name in ["WGS1984"]: proj4_params["ellps"] = "WGS84"
            else: print("Error getting the CRS from {}: Script did not account for the {} ellipsoid yet.".format(xml_path, ellips_name))
        # Create a CRS with the proj4 parameters (https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.CRS.html#cartopy.crs.CRS.__init__).
        xml_crs = ccrs.CRS(proj4_params)
        xml_epsg_code = xml_crs.to_epsg()
        # Return the created CRS.
        if xml_epsg_code is not None:
            dataset_crs = ccrs.epsg(xml_epsg_code)
            return dataset_crs
        else:
            return xml_crs
    else:
        # Datasets should have the same CRS for all their data files, so return the found CRS if it was already previously computed.
        return dataset_crs

def convert_csv_txt_data_into_geojson(file_path: str, geojson_path: str) -> None:
    """
    Creates and saves a GeoJSON file containing Points for each data point in the given dataframe.
    
    Args:
        file_path (str): Path to the file containing data points
        geojson_path (str): Path to the newly created GeoJSON file, which is a FeatureCollection of Points
    """
    if not os.path.exists(geojson_path):
        # Get the file's CRS in case there's only point data in the dataset (only converting ASCII grid files requires the returned CRS).
        global dataset_crs
        if dataset_crs is None: _ = get_crs_from_xml_file(file_path)
        # Read the data file as a DataFrame, drop any rows with all NaN values, and replace any NaN values with "N/A".
        dataframe = pd.read_csv(file_path).dropna(axis = 0, how = "all").fillna("N/A")
        # Ignore any unnamed columns.
        dataframe = dataframe.loc[:, ~dataframe.columns.str.match("Unnamed")]
        # Get the latitude and longitude column names.
        latitude_col, *_ = [col for col in dataframe.columns if "lat" in col.lower()]
        longitude_col, *_ = [col for col in dataframe.columns if "lon" in col.lower()]
        # Convert the DataFrame into a GeoDataFrame.
        geodataframe = gpd.GeoDataFrame(
            data = dataframe,
            geometry = gpd.points_from_xy(
                x = dataframe[longitude_col],
                y = dataframe[latitude_col],
                crs = ccrs.PlateCarree()
            )
        )
        # Save the GeoDataFrame into a GeoJSON file to skip converting the data file again.
        geodataframe.to_file(geojson_path, driver = "GeoJSON")

def convert_ascii_grid_data_into_geotiff(file_path: str, geotiff_path: str) -> None:
    """
    Converts an ASCII grid file into a GeoTIFF file.

    Args:
        file_path (str): Path to the ASCII grid file
        geotiff_path (str): Path to the newly created GeoTIFF file
    """
    if not os.path.exists(geotiff_path):
        # Get the CRS of the data file using its XML file.
        file_crs = get_crs_from_xml_file(file_path)
        # Add found CRS to the data file.
        dataset = rxr.open_rasterio(file_path)
        dataset.rio.write_crs(file_crs, inplace = True)
        # Save the data as a GeoTIFF.
        dataset.rio.to_raster(
            raster_path = geotiff_path,
            driver = "GTiff"
        )

def convert_transect_data_into_geojson(file_path: str, geojson_path: str) -> None:
    """
    Creates and saves a GeoJSON file containing LineStrings for each transect in the given transect file.
    Note: This method was specifically written to read transect files that are in the same format as ew_lines.txt.

    Args:
        file_path (str): Path to the file containing transect data
        geojson_path (str): Path to the newly created GeoJSON file
    """
    if not os.path.exists(geojson_path):
        # Create the GeoData folder if it doesn't exist yet.
        geodata_dir_path, _ = os.path.split(geojson_path)
        if not os.path.isdir(geodata_dir_path): os.makedirs(geodata_dir_path)
        # Create a FeatureCollection of LineStrings based on the data file.
        features_list = []
        with open(file_path, "r") as file:
            transect_feature = None
            for line in file:
                [point_id, x, y, _] = line.split(",")
                point = [float(x), float(y)]
                if transect_feature is None:
                    # Initialize a new transect feature.
                    id = int("".join([char for char in point_id if char.isdigit()]))
                    transect_feature = {
                        "type": "Feature",
                        "properties": {transect_geojson_id_property: id},
                        "geometry": {
                            "type": "LineString",
                            "coordinates": []
                        }
                    }
                    # Add the transect's start point.
                    transect_feature["geometry"]["coordinates"].append(point)
                    transect_feature["properties"][transect_geojson_start_point_property] = "({}, {})".format(x, y)
                else:
                    # Add the transect's end point.
                    transect_feature["geometry"]["coordinates"].append(point)
                    transect_feature["properties"][transect_geojson_end_point_property] = "({}, {})".format(x, y)
                    # Save the transect to the FeatureCollection.
                    features_list.append(transect_feature)
                    # Reset the feature for the next transect.
                    transect_feature = None
        # Convert the FeatureCollection into a GeoJSON.
        geodataframe = gpd.GeoDataFrame.from_features(
            {"type": "FeatureCollection", "features": features_list},
            crs = dataset_crs if dataset_crs is not None else ccrs.PlateCarree()
        )
        # Save the GeoJSON file to skip converting the data file again.
        geodataframe.to_file(geojson_path, driver = "GeoJSON")

def preprocess_data(src_dir_path: str, dest_dir_path: str, dir_level: int = 1) -> None:
    """
    Recursively preprocess all the data in the given data directory.

    Args:
        src_dir_path (str): Location of the source directory containing data files to convert into formats that are compatible with DataMap
        dest_dir_path (str): Location of the destination directory where newly created subdirectories or processed data files are saved
        dir_level (int): Tree level of the given source data directory
            ^ used to prevent creating nested subdirectories (not compatible with DataMap) within the root destination directory
            ^ root destination directory is the only directory that contains subdirectories, each holding processed data files
    """
    # Create a new directory in the destination directory if it's still compatible with DataMap after it gets added.
    src_dir_name = os.path.basename(src_dir_path)
    new_dest_dir_path = os.path.join(dest_dir_path, src_dir_name)
    if dir_level < 3:
        if not os.path.exists(new_dest_dir_path): os.makedirs(new_dest_dir_path)
    else:
        new_dest_dir_path = dest_dir_path
    # Look through source directory for raw data files.
    subdir_level = dir_level + 1
    print("Searching for data files to convert from {}...".format(src_dir_path))
    for file in os.listdir(src_dir_path):
        file_path = os.path.join(src_dir_path, file)
        # Look through subdirectories for raw data files as well.
        if os.path.isdir(file_path):
            preprocess_data(src_dir_path = file_path, dest_dir_path = new_dest_dir_path, dir_level = subdir_level)
        elif file != sb_download_output_json_name:
            name, extension = os.path.splitext(file)
            file_format = extension.lower()
            # Convert data file into a format that is more compatible for DataMap.
            if file_format in [".csv", ".txt"]:
                geojson_file_path = os.path.join(new_dest_dir_path, name + ".geojson")
                print("\t{} -> {}".format(file, geojson_file_path))
                if transects_dir_exists and (src_dir_name == transects_subdir_name):
                    convert_transect_data_into_geojson(file_path, geojson_file_path)
                else:
                    convert_csv_txt_data_into_geojson(file_path, geojson_file_path)
                    buffer_config[geojson_file_path] = 3
            elif file_format == ".asc":
                geotiff_file_path = os.path.join(new_dest_dir_path, name + ".tif")
                print("\t{} -> {}".format(file, geotiff_file_path))
                convert_ascii_grid_data_into_geotiff(file_path, geotiff_file_path)
                buffer_config[geotiff_file_path] = 0
            elif file_format in [".geojson", ".tif", ".tiff"]:
                print("TODO: copy file")
            elif file_format not in [".xml", ".png"]:
                print("Error converting {}: Data files with the {} file format are not supported yet.".format(file_path, file_format))

# -------------------------------------------------- Main Program --------------------------------------------------
if __name__ == "__main__":
    parent_data_dir_path = os.path.abspath("./utils")
    unprocessed_data_dirs = [file for file in os.listdir(parent_data_dir_path) if os.path.isdir(os.path.join(parent_data_dir_path, file)) and (file != "__pycache__")]
    num_unprocessed_data_dirs = len(unprocessed_data_dirs)
    if num_unprocessed_data_dirs > 0:
        # 1. Get the path to the data directory that the user wants to preprocess.
        print("Data directories to preprocess:")
        for i, data_dir in enumerate(unprocessed_data_dirs): print("\t[{}] {}".format(i + 1, data_dir))
        dir_index = input("Please enter your numeric choice: ")
        # 2. Check if the user inputted a valid choice.
        if dir_index.isnumeric() and (0 < int(dir_index) <= num_unprocessed_data_dirs):
            selected_data_dir = unprocessed_data_dirs[int(dir_index) - 1]
            transects_input = input("Does {} contain a `Transects` subdirectory?\n\t[y] Yes\n\t[n] No\nPlease enter your alphabetic choice: ".format(selected_data_dir))
            if transects_input in ["y", "n"]:
                if transects_input == "y": transects_dir_exists = True
                elif transects_input == "n": transects_dir_exists = False
                # 3. Iterate through data directories and convert data files into formats that are compatible with DataMap.
                data_dir_path = os.path.join(parent_data_dir_path, selected_data_dir)
                print("All data from {} will be preprocessed momentarily...".format(data_dir_path))
                root_output_dir_path = os.path.abspath("./data")
                preprocess_data(src_dir_path = data_dir_path, dest_dir_path = root_output_dir_path)
                # # 4. Rename directories if they're outputted from download_sciencebase_data.py.
                # sb_download_output_json_file_path = os.path.join(data_dir_path, sb_download_output_json_name)
                # if os.path.exists(sb_download_output_json_file_path):
                #     # Open the JSON file that maps ScienceBase item IDs to their titles.
                #     json_file = open(sb_download_output_json_file_path)
                #     item_id_to_title = json.load(json_file)
                #     # Replace each item ID directory with the item's title instead.

                #     # Get the title of the root item that contained all the downloaded items.
                #     selected_data_dir = item_id_to_title[selected_data_dir]
                # 4. Save data's CRS in an outputted data_info.json file.
                # TODO: account for when there's no CRS found (b/c data files have been converted already)
                data_info = {data_from_download_script_property: False}
                if dataset_crs is None: data_info[dataset_crs_property] = None
                else: data_info[dataset_crs_property] = dataset_crs.to_epsg()
                # Also save contents from sciencebase_id_to_title.json if the data was downloaded with download_sciencebase_data.py.
                sb_download_output_json_file_path = os.path.join(data_dir_path, sb_download_output_json_name)
                if os.path.exists(sb_download_output_json_file_path):
                    data_info[data_from_download_script_property] = True
                    # Open the JSON file that maps ScienceBase item IDs to their titles.
                    json_file = open(sb_download_output_json_file_path)
                    item_id_to_title = json.load(json_file)
                    data_info.update(item_id_to_title)
                preprocessed_data_path = os.path.join(root_output_dir_path, selected_data_dir)
                with open(os.path.join(preprocessed_data_path, outputted_dataset_json_name), "w") as dataset_json_file:
                    json.dump(data_info, dataset_json_file, indent = 4)
                # 5. Save buffer configurations for each data file, which is later used to extract data along or near a transect.
                with open(os.path.join(preprocessed_data_path, outputted_buffer_json_name), "w") as buffer_json_file:
                    json.dump(buffer_config, buffer_json_file, indent = 4)
                print("Converting data complete! All preprocessed data files are saved in {}.".format(preprocessed_data_path))
            else:
                print("Invalid choice: Your choice {} did not match any of the ones provided above. Please run this script again with a valid alphabetic choice.".format(transects_input))
        else:
            print("Invalid choice: Your choice {} did not match any of the ones provided above. Please run this script again with a valid numeric choice.".format(dir_index))
    else:
        print("Data not found: There are no data directories to preprocess. Make sure your data is placed in {}.".format(parent_data_dir_path))