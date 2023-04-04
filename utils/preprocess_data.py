# cd C:\Users\Venuxk\Projects\data-visualizer
# conda activate visualizer
# python ./utils/preprocess_data.py

# Standard library imports
import os
import json
import xml.etree.ElementTree as ET
import shutil
import math
from collections import defaultdict

# External dependencies imports
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import rioxarray as rxr
import dask_geopandas
from download_sciencebase_data import outputted_json_name as sb_download_output_json_name

# -------------------------------------------------- Constants (should match the constants used in DataMap.py) --------------------------------------------------
outputted_collection_json_name = "collection_info.json"
collection_epsg_property = "epsg"
collection_data_categories_property = "categories"
outputted_buffer_json_name = "buffer_config.json"

transects_subdir_name = "Transects"
transect_geojson_id_property = "Transect ID"
transect_geojson_start_point_property = "Start Point"
transect_geojson_end_point_property = "End Point"

elwha_river_delta_item_id = "5a01f6d0e4b0531197b72cfe"
elwha_epsg = 32148

# -------------------------------------------------- Global Variables --------------------------------------------------
collection_dir_name = None
collection_crs = None
transects_dir_exists = False
collection_info = {
    collection_epsg_property: 4326,
    collection_data_categories_property: defaultdict(list)
}
buffer_config = {}

# -------------------------------------------------- Helper Methods --------------------------------------------------
def get_crs_from_xml_file(file_path: str) -> ccrs:
    """
    Read a data file's XML (extensible markup language) to get its CRS (coordinate reference system).

    Args:
        file_path (str): Path to a data file, which is in the same directory as its XML file
            ^ XML file contains information about the data's CRS
    """
    global collection_crs
    if collection_crs is None:
        if os.path.exists(file_path):
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
                collection_crs = ccrs.epsg(xml_epsg_code)
                return collection_crs
            else:
                return xml_crs
        else:
            return ccrs.epsg(4326)
    else:
        # Collections should have the same CRS for all their data files, so return the found CRS if it was already previously computed.
        return collection_crs

def convert_csv_txt_data_into_parquet(file_path: str, parquet_path: str) -> None:
    """
    Creates and saves a directory with Parquet files containing Points for each data point in the given dataframe.
    
    Args:
        file_path (str): Path to the file containing data points
        parquet_path (str): Path to the newly created Parquet directory, which is a FeatureCollection of Points
    """
    if not os.path.exists(parquet_path):
        # Get the file's CRS in case there's only point data in the collection (only converting ASCII grid files requires the returned CRS).
        global collection_crs
        if collection_crs is None: _ = get_crs_from_xml_file(file_path)
        # Read the data file as a pandas DataFrame, drop any rows with all NaN values.
        pd_dataframe = pd.read_csv(file_path).dropna(axis = 0, how = "all")
        # Ignore any unnamed columns.
        pd_dataframe = pd_dataframe.loc[:, ~pd_dataframe.columns.str.match("Unnamed")]
        # Get the latitude and longitude column names.
        latitude_col, *_ = [col for col in pd_dataframe.columns if "lat" in col.lower()]
        longitude_col, *_ = [col for col in pd_dataframe.columns if "lon" in col.lower()]
        # Convert the pandas DataFrame into a geopandas GeoDataFrame.
        gpd_geodataframe = gpd.GeoDataFrame(
            data = pd_dataframe,
            geometry = gpd.points_from_xy(
                x = pd_dataframe[longitude_col],
                y = pd_dataframe[latitude_col],
                crs = ccrs.PlateCarree()
            )
        )
        num_parquet_partitions = math.ceil(gpd_geodataframe.memory_usage(deep = True).sum() / 1e9)
        # Convert the geopandas GeoDataFrame into a Dask GeoDataFrame.
        dask_geodataframe = dask_geopandas.from_geopandas(gpd_geodataframe, npartitions = num_parquet_partitions)
        # Spatially optimize partitions of the DaskGeoDataFrame by using a Hilbert R-tree packing method, which groups neighboring data into the same partition.
        dask_geodataframe = dask_geodataframe.spatial_shuffle(by = "hilbert", npartitions = num_parquet_partitions)
        # Save the spatially optimized Dask GeoDataFrame.
        dask_geodataframe.to_parquet(path = parquet_path)

def convert_ascii_grid_data_into_geotiff(file_path: str, geotiff_path: str) -> None:
    """
    Converts an ASCII grid file into a cloud optimized GeoTIFF file.

    Args:
        file_path (str): Path to the ASCII grid file
        geotiff_path (str): Path to the newly created cloud optimized GeoTIFF file
    """
    if not os.path.exists(geotiff_path):
        # Get the CRS of the data file using its XML file.
        file_crs = get_crs_from_xml_file(file_path)
        # Add found CRS to the data file.
        dataset = rxr.open_rasterio(file_path)
        dataset.rio.write_crs(file_crs, inplace = True)
        # Save the data as a cloud optimized GeoTIFF.
        dataset.rio.to_raster(
            raster_path = geotiff_path,
            tiled = True,
            driver = "COG"
        )

def convert_transect_data_into_parquet(file_path: str, parquet_path: str) -> None:
    """
    Creates and saves a directory with Parquet files containing LineStrings for each transect in the given transect file.
    Note: This method was specifically written to read transect files that are in the same format as ew_lines.txt.

    Args:
        file_path (str): Path to the file containing transect data
        parquet_path (str): Path to the newly created Parquet directory
    """
    if not os.path.exists(parquet_path):
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
        # Convert the FeatureCollection into a geopandas GeoDataFrame.
        gpd_geodataframe = gpd.GeoDataFrame.from_features({"type": "FeatureCollection", "features": features_list})
        if collection_dir_name and (collection_dir_name == elwha_river_delta_item_id or collection_dir_name == "points_test"): transect_epsg = elwha_epsg
        else: transect_epsg = 4326
        gpd_geodataframe = gpd_geodataframe.set_crs(transect_epsg)
        num_parquet_partitions = math.ceil(gpd_geodataframe.memory_usage(deep = True).sum() / 1e9)
        # Convert the geopandas GeoDataFrame into a Dask GeoDataFrame.
        dask_geodataframe = dask_geopandas.from_geopandas(gpd_geodataframe, npartitions = num_parquet_partitions)
        # Save the spatially optimized Dask GeoDataFrame.
        dask_geodataframe.to_parquet(path = parquet_path)

def set_readable_file_name(file_path: str) -> None:
    """
    Sets the human-readable name for the given data file, which will be saved to the JSON file that stores all information about the collection.

    Args:
        file_path (str): Path to the data file that needs to be assigned a human-readable name
    """
    global collection_info
    if collection_dir_name and collection_dir_name == elwha_river_delta_item_id:
        # Assign each Elwha River delta data file's name by month/year collected and the data type.
        month = {
            "jan": "January",
            "feb": "February",
            "mar": "March",
            "apr": "April",
            "may": "May",
            "june": "June",
            "july": "July",
            "aug": "August",
            "sept": "September",
            "oct": "October",
            "nov": "November",
            "dec": "December"
        }
        type = {
            "1m.tif": "Digital Elevation Model (1-meter resolution DEM)",
            "5m.tif": "Digital Elevation Model (5-meter resolution DEM)",
            "grainsize.parq": "Surface-Sediment Grain-Size Distribution",
            "kayak.parq": "Bathymetry (Kayak)",
            "pwc.parq": "Bathymetry (Personal Watercraft)",
            "bathy.parq": "Bathymetry (Personal Watercraft)",
            "topo.parq": "Topography"
        }
        file_name = os.path.basename(file_path)
        mon, yr, typ = "", "", ""
        for subname in file_name.split("_"):
            if "ew" in subname:
                yr = "20" + subname[2:]
            elif subname in month:
                mon = month[subname]
            elif subname in type:
                typ = type[subname]
                collection_info[collection_data_categories_property][typ].append(file_path)
        collection_info[file_path] = "{} {}".format(mon if mon else "July", yr)

def sort_data_files_by_collection_date(categories_to_files: dict) -> None:
    """
    Sorts given data categories and their files by the collection date.

    Args:
        categories_to_files (dict): Dictionary mapping each data category (key) to a list of paths (value) that leads to data files belonging to the category
    """
    global collection_info
    if collection_dir_name == elwha_river_delta_item_id:
        months = {
            "January": 1, "February": 2, "March": 3,
            "April": 4, "May": 5, "June": 6,
            "July": 7, "August": 8, "September": 9,
            "October": 10, "November": 11, "December": 12
        }
        for category, file_paths in categories_to_files.items():
            # Create a collection date ID for each data file in the category.
            file_ids, ids_to_paths = [], {}
            for path in file_paths:
                if path in collection_info:
                    collection_date = collection_info[path].split(" - ")[0]
                    month_name, year = collection_date.split()
                    file_id = date_id = str(months[month_name] + (int(year) * 100))
                    file_ids.append(date_id)
                else:
                    file_id = file_name = os.path.basename(path)
                    file_ids.append(file_name)
                ids_to_paths[file_id] = path
            # Sort data file paths by collection month and year.
            sorted_file_paths = [ids_to_paths[id] for id in sorted(file_ids)]
            # Save the sorted list of data file paths.
            collection_info[collection_data_categories_property][category] = sorted_file_paths

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
    global collection_dir_name
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
            # Create a subdirectory within the outputted directory with the same name if the inputted item doesn't have any children (only has attached files).
            # ^ ensures all preprocessed data has exactly one outputted directory containing subdirectories with data files (no nested subdirectories)
            if dir_level == 1:
                subdir_path = os.path.join(new_dest_dir_path, src_dir_name)
                if not os.path.exists(subdir_path): os.makedirs(subdir_path)
                new_dest_dir_path = subdir_path
            # Convert data file into a format that is more compatible for DataMap.
            if file_format in [".csv", ".txt"]:
                parquet_files_path = os.path.join(new_dest_dir_path, name + ".parq")
                print("\t{} -> {}".format(file, parquet_files_path))
                if transects_dir_exists and (src_dir_name == transects_subdir_name):
                    convert_transect_data_into_parquet(file_path, parquet_files_path)
                else:
                    convert_csv_txt_data_into_parquet(file_path, parquet_files_path)
                    buffer_config[parquet_files_path] = 3
                    set_readable_file_name(parquet_files_path)
            elif file_format == ".asc":
                geotiff_file_path = os.path.join(new_dest_dir_path, name + ".tif")
                print("\t{} -> {}".format(file, geotiff_file_path))
                convert_ascii_grid_data_into_geotiff(file_path, geotiff_file_path)
                buffer_config[geotiff_file_path] = 0
                set_readable_file_name(geotiff_file_path)
            elif file_format in [".parq", ".tif", ".tiff"]:
                geodata_file_path = os.path.join(new_dest_dir_path, file)
                if not os.path.exists(geodata_file_path): shutil.copy2(file_path, new_dest_dir_path)
                set_readable_file_name(geodata_file_path)
            elif file_format not in [".xml", ".png"]:
                print("Error converting {}: Data files with the {} file format are not supported yet.".format(file_path, file_format))

# -------------------------------------------------- Main Program --------------------------------------------------
if __name__ == "__main__":
    parent_data_dir_path = os.path.relpath("./utils")
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
            collection_dir_name = selected_data_dir
            transects_input = input("Does {} contain a `Transects` subdirectory?\n\t[y] Yes\n\t[n] No\nPlease enter your alphabetic choice: ".format(selected_data_dir))
            if transects_input in ["y", "n"]:
                if transects_input == "y": transects_dir_exists = True
                elif transects_input == "n": transects_dir_exists = False
                # 3. Iterate through data directories and convert data files into formats that are compatible with DataMap.
                data_dir_path = os.path.join(parent_data_dir_path, selected_data_dir)
                print("All data from {} will be preprocessed momentarily...".format(data_dir_path))
                root_output_dir_path = os.path.relpath("./data")
                preprocess_data(src_dir_path = data_dir_path, dest_dir_path = root_output_dir_path)
                # 4. Save data's CRS in an outputted collection_info.json file.
                if collection_crs is not None: collection_info[collection_epsg_property] = collection_crs.to_epsg()
                # Also save contents from sciencebase_id_to_title.json if the data was downloaded with download_sciencebase_data.py.
                sb_download_output_json_file_path = os.path.join(data_dir_path, sb_download_output_json_name)
                if os.path.exists(sb_download_output_json_file_path):
                    # Open the JSON file that maps ScienceBase item IDs to their titles.
                    json_file = open(sb_download_output_json_file_path)
                    item_id_to_title = json.load(json_file)
                    collection_info.update(item_id_to_title)
                # Sort data files by their collection date, if possible.
                sort_data_files_by_collection_date(collection_info[collection_data_categories_property])
                preprocessed_data_path = os.path.join(root_output_dir_path, selected_data_dir)
                with open(os.path.join(preprocessed_data_path, outputted_collection_json_name), "w") as collection_json_file:
                    json.dump(collection_info, collection_json_file, indent = 4)
                # 5. Save buffer configurations for each data file, which is later used to extract data along or near a transect.
                with open(os.path.join(preprocessed_data_path, outputted_buffer_json_name), "w") as buffer_json_file:
                    json.dump(buffer_config, buffer_json_file, indent = 4)
                print("Converting data complete! All preprocessed data files are saved as a collection in {}.".format(preprocessed_data_path))
            else:
                print("Invalid choice: Your choice {} did not match any of the ones provided above. Please run this script again with a valid alphabetic choice.".format(transects_input))
        else:
            print("Invalid choice: Your choice {} did not match any of the ones provided above. Please run this script again with a valid numeric choice.".format(dir_index))
    else:
        print("Data not found: There are no data directories to preprocess. Make sure your data is placed in {}.".format(parent_data_dir_path))