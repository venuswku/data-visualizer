# Standard library imports
import os
import xml.dom.minidom

# External dependencies imports
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import rioxarray as rxr

# -------------------------------------------------- Global Variables --------------------------------------------------

# -------------------------------------------------- Helper Methods --------------------------------------------------
def convert_csv_txt_data_into_geojson(file_path: str, geojson_path: str) -> None:
    """
    Creates and saves a GeoJSON file containing Points for each data point in the given dataframe.
    
    Args:
        file_path (str): Path to the file containing data points
        geojson_path (str): Path to the newly created GeoJSON file, which is a FeatureCollection of Points
    """
    if not os.path.exists(geojson_path):
        # Read the data file as a DataFrame and replace any NaN values with "N/A".
        dataframe = pd.read_csv(file_path).fillna("N/A")
        # Ignore any unnamed columns.
        dataframe = dataframe.loc[:, ~dataframe.columns.str.match("Unnamed")]
        # Get the latitude and longitude column names.
        latitude_col, *other_lat_cols = [col for col in dataframe.columns if "lat" in col.lower()]
        longitude_col, *other_long_cols = [col for col in dataframe.columns if "lon" in col.lower()]
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
        # Read the data file's XML (extensible markup language) to get its CRS (coordinate reference system).
        xml_file = xml.dom.minidom.parse("")
        column_element = xml_file.getElementsByTagName("planar")
        column_names = [element.firstChild.data for element in column_elements]
        # https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.CRS.html#cartopy.crs.CRS.__init__
        geotiff_crs = ccrs.CRS({
            "proj": "lcc",
            "lat_1": 47.5,
            "lat_2": 48.73333333333333,
            "lon_0": -120.8333333333333,
            "lat_0": 47.0,
            "x_0": 500000.0,
            "y_0": 0.0,
            "units": "m",
            "datum": "NAD83",
            "ellps": "GRS80",
            "no_defs": True,
            "type": "crs"
        })
        # Add projection to the data file based on the XML's metadata.
        dataset = rxr.open_rasterio(file_path)
        dataset.rio.write_crs(self._epsg, inplace = True)
        # Save the data as a GeoTIFF.
        dataset.rio.to_raster(
            raster_path = geotiff_path,
            driver = "GTiff"
        )

def convert_transect_data_into_geojson(self, file_path: str, geojson_path: str) -> None:
    """
    Creates and saves a GeoJSON file containing LineStrings for each transect in the given transect file.
    Note this method was specifically written to read transect files that are in the same format as ew_lines.txt.

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
                        "properties": {self._transects_id_col_name: id},
                        "geometry": {
                            "type": "LineString",
                            "coordinates": []
                        }
                    }
                    # Add the transect's start point.
                    transect_feature["geometry"]["coordinates"].append(point)
                    transect_feature["properties"][self._transect_start_point_prop_name] = "({}, {})".format(x, y)
                else:
                    # Add the transect's end point.
                    transect_feature["geometry"]["coordinates"].append(point)
                    transect_feature["properties"][self._transect_end_point_prop_name] = "({}, {})".format(x, y)
                    # Save the transect to the FeatureCollection.
                    features_list.append(transect_feature)
                    # Reset the feature for the next transect.
                    transect_feature = None
        # Convert the FeatureCollection into a GeoJSON.
        geodataframe = gpd.GeoDataFrame.from_features(
            {"type": "FeatureCollection", "features": features_list},
            crs = self._epsg
        )#.to_crs(crs = self._default_crs)
        # Save the GeoJSON file to skip converting the data file again.
        geodataframe.to_file(geojson_path, driver = "GeoJSON")

# -------------------------------------------------- Main Program --------------------------------------------------
# 1. Get the path to the data directory that the user wants to preprocess.

# 2. Iterate through data directories and convert data files into formats that are compatible with DataMap.

# 3. 
# Rename directories if they're outputted from download_sciencebase_data.py.
if os.path.exists(path):