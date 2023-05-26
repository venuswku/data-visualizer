# Standard library imports
import os
import json
from collections import Counter
import asyncio
import datetime as dt
import time

# External dependencies imports
import param
import panel as pn
import holoviews as hv
import rioxarray as rxr
import geopandas as gpd
import pandas as pd
import dask_geopandas
from shapely.geometry import Point, LineString
import cartopy.crs as ccrs
import numpy as np
from bokeh.models.formatters import PrintfTickFormatter
from bokeh.palettes import Turbo256
from .DataMap import DataMap

### PopupModal is used to display a time-series plot or any other data/message in the app's modal. ###
class PopupModal(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    clicked_transects_info = param.Dict(default = {}, label = "Information About the Recently Clicked Transect(s)")
    start_data_collection_date = param.Date(default = dt.date(2010, 9, 1), label = "Start Date of Collected Data")
    end_data_collection_date = param.Date(default = dt.date(2018, 8, 1), label = "End Date of Collected Data")
    data_category = param.Selector(label = "Category of Selected Time-Series Data")
    displayed_data_file = param.Selector(label = "Path to a Displayed Time-Series Data File")
    
    update_collection_dir_path = param.Event(label = "Action that Triggers the Updating of the Collection Directory and Its Related Objects")
    update_buffer_config = param.Event(label = "Action that Triggers Updating the Buffer Config File")
    update_accordion_section = param.Event(label = "Indicator for Updating the PopupModal's Accordion Sections")
    download_time_series = param.Event(label = "Action that Triggers Downloading the Computed Time-Series for a Selected Transect")

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, data_map: DataMap, elwha_collection_name: str, **params) -> None:
        """
        Creates a new instance of the PopupModal class with its instance variables.

        Args:
            data_map (DataMap): Instance containing methods for converting data files to allow quicker loading onto a map
            elwha_collection_name (str): ScienceBase item id for the root item containing all the Elwha River Delta data
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._data_map = data_map
        self._elwha_collection_name = elwha_collection_name
        
        # _outputs_dir_path = path to directory containing all downloaded time-series outputs
        self._outputs_dir_path = os.path.relpath("./outputs")
        # _dist_col_name = name of the column that stores the x-axis values (distance from shore) for the time-series plot
        self._dist_col_name = "Across-Shore Distance (m)"
        # _default_y_axis_data_col_name = default name of the column that stores the y-axis values for the time-series plot (default is often used for data in ASCII grid files)
        self._default_y_axis_data_col_name = "Elevation (m)"
        # _point_type_col_name = name of the column that stores the type of transect point (either start or end)
        self._point_type_col_name = "Point Type"
        # The following list of constant variables are keys that appear in the dictionary that DataMap sends into PopupModal's _clicked_transects_pipe stream.
        # ^ When the _clicked_transects_pipe stream gets sent a new dictionary, the dictionary is passed into the _create_time_series_plot() callback as the `data` keyword argument.
        [self._clicked_transects_file, self._num_clicked_transects, self._clicked_transects_crs, self._clicked_transects_longitude_col,
        self._clicked_transects_latitude_col, self._clicked_transects_table_cols, self._clicked_transects_id_col] = data_map.clicked_transects_info_keys
        # _preprocessed_data_buffer_output = name of the buffer config file outputted by the script from utils/preprocess_data.py (should be the same as `outputted_buffer_json_name`)
        self._preprocessed_data_buffer_output = "buffer_config.json"
        # _placeholder_data_category = default option for not selecting a data category
        self._placeholder_data_category = "None"
        # _placeholder_displayed_data = default data file option when no data category has been selected
        # ^ else the last selected data file (most likely the most recent data) would be the default
        self._placeholder_displayed_data = "None"

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------        
        # _collection_dir_path = path to the directory containing all the data files used for the time-series
        self._collection_dir_path = data_map.selected_collection_dir_path
        # _y_axis_data_col_name = name of the column that stores the y-axis values for the time-series plot
        self._y_axis_data_col_name = self._default_y_axis_data_col_name
        # _time_series_plot = time-series plot for data collected along the most recently clicked transect/path
        self._time_series_plot = pn.pane.HoloViews(object = None)
        # _category_multiselect_widget = dictionary mapping each data category (key) to its MultiSelect widget (value)
        self._category_multiselect_widget = {}
        # _user_selected_data_files = list of paths to data files that are used for the time-series
        self._user_selected_data_files = []
        # _data_file_colors = dictionary mapping each selected data file's path (key) to their color (value) in the time-series
        self._data_file_colors = {}
        # _buffers = dictionary mapping each data file's path (key) to the selected transect's buffer/search radius (value) when extracting data around the transect
        self._buffers = {}
        # _buffer_widget_file_path = dictionary mapping the name of each float input widget (key) to the path (value) of the data file that uses this buffer when extracting data around the transect
        self._buffer_widget_file_path = {}

        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # _data_files_widgets = column layout containing widgets for the "Time-Series Data" accordion section
        self._data_files_widgets = pn.Column(objects = [])
        # _time_series_data_wiki_info_button = button that opens a tab to the GitHub Wiki page that describes how to use the time-series data controls
        self._time_series_data_wiki_info_button = pn.widgets.Button(name = "\u2139", button_type = "light", width = 30)
        self._time_series_data_wiki_info_button.js_on_click(
            code = "window.open('https://github.com/venuswku/data-visualizer/wiki/Sidebar-Controls#time-series-data')"
        )
        # _start_data_collection_date_picker = date picker for choosing the start date of the collected data, which will be the earliest date in the time-series
        self._start_data_collection_date_picker = pn.widgets.DatePicker.from_param(
            parameter = self.param.start_data_collection_date,
            name = "Time-Series Start Date"
        )
        # _end_data_collection_date_picker = date picker for choosing the end date of the collected data, which will be the latest date in the time-series
        self._end_data_collection_date_picker = pn.widgets.DatePicker.from_param(
            parameter = self.param.end_data_collection_date,
            name = "Time-Series End Date"
        )
        # _data_category_select = widget for selecting a category of data used for the time-series
        self._data_category_select = pn.widgets.Select.from_param(
            parameter = self.param.data_category,
            name = "Time-Series Data Category",
            options = [self._placeholder_data_category],
            value = self._placeholder_data_category
        )
        # _displayed_data_select = widget for selecting a path to a data file that can be displayed on a map
        self._displayed_data_select = pn.widgets.Select.from_param(
            parameter = self.param.displayed_data_file,
            name = "",
            options = {self._placeholder_displayed_data: None}
        )
        # _time_series_data_constant_widgets = list of widgets that always appear at the top of the "Time-Series Data" accordion section
        self._time_series_data_constant_widgets = [
            pn.Row(
                pn.widgets.StaticText(value = "The data files highlighted below are used when creating a time-series of how data changes over time along a chosen transect.", width = 250),
                self._time_series_data_wiki_info_button
            ),
            pn.pane.Markdown(object = "**Displayed Data**", sizing_mode = "stretch_width", margin = (10, 10, -15, 10)),
            pn.widgets.StaticText(value = "In order to prevent overcrowding data, select only one of the highlighted data files to display on the map."),
            self._displayed_data_select,
            pn.pane.Markdown(object = "**All Available Data**", sizing_mode = "stretch_width", margin = (10, 10, -15, 10))
        ]
        
        # _transect_search_radius_widgets = column layout containing widgets for the "Transect Search Radius" accordion section
        self._transect_search_radius_widgets = pn.Column(objects = [])
        # _transect_search_radius_wiki_info_button = button that opens a tab to the GitHub Wiki page that describes how to use the transect search radius controls
        self._transect_search_radius_wiki_info_button = pn.widgets.Button(name = "\u2139", button_type = "light", width = 30)
        self._transect_search_radius_wiki_info_button.js_on_click(
            code = "window.open('https://github.com/venuswku/data-visualizer/wiki/Sidebar-Controls#transect-search-radius')"
        )
        # _update_buffer_config_file_button = button for updating the collection's buffer config file with values from the buffer float widgets for each data file
        self._update_buffer_config_file_button = pn.widgets.Button.from_param(
            parameter = self.param.update_buffer_config,
            name = "Save Current Values to {}".format(self._preprocessed_data_buffer_output),
            button_type = "primary", disabled = True
        )
        # _transect_search_radius_constant_widgets = list of widgets that always appear at the top of the "Transect Search Radius" accordion section
        self._transect_search_radius_constant_widgets = [
            pn.Row(
                pn.widgets.StaticText(value = "Adjust the search radius for extracting time-series data around a selected transect.", width = 250),
                self._transect_search_radius_wiki_info_button
            ),
            self._update_buffer_config_file_button
        ]

        # _modal_heading = list containing markdown objects for the modal heading (title for first markdown, details for second markdown)
        self._modal_heading = pn.Column(objects = [pn.pane.Markdown(object = ""), pn.pane.Markdown(object = "", margin = (-20, 5, 0, 5))])
        # _clicked_transects_table = table containing information about the clicked transect(s)'s start and end points
        self._clicked_transects_table = pn.widgets.DataFrame(
            value = pd.DataFrame(),
            name = "Selected Transect(s) Data",
            show_index = True, auto_edit = False, text_align = "center",
            sizing_mode = "stretch_width", margin = (-20, 5, 10, 5)
        )
        # _time_series_dataframes = list of pandas DataFrames containing time-series data for each collection date
        self._time_series_dataframes = []
        # _time_series_download_button = button for downloading the time-series plot
        self._time_series_download_button = pn.widgets.Button.from_param(
            parameter = self.param.download_time_series,
            name = "Save Time-Series",
            button_type = "primary", visible = False
        )
        # Initialize widgets that depend on the selected collection from DataMap.
        self._update_collection_objects()

    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    def _categorize_data_files(self) -> None:
        """
        Assigns the collection's data files to the multiselect widget that corresponds to their data category.
        """
        widgets = []
        collection_name = os.path.basename(self._collection_dir_path)
        # Group data files by the type of data for the Elwha River delta collection.
        if collection_name == self._elwha_collection_name:
            # Create a widget for each data category in the collection.
            for category_name, file_paths in self._data_map.selected_collection_json_info["categories"].items():
                options_dict = {}
                for path in file_paths:
                    option_name = self._data_map.selected_collection_json_info[path]
                    options_dict[option_name] = path
                category_multiselect = pn.widgets.MultiSelect(name = category_name, options = options_dict, value = [], disabled = True)
                widgets.append(category_multiselect)
                self._category_multiselect_widget[category_name] = category_multiselect
        else:
            single_category_multiselect = pn.widgets.MultiSelect(options = self._data_map.data_file_options, value = [], disabled = True)
            self._category_multiselect_widget["Other"] = single_category_multiselect
            widgets.append(single_category_multiselect)
        # Assign new widgets for allowing the user to choose time-series data files.
        self._data_files_widgets.objects = self._time_series_data_constant_widgets + widgets
    
    def _within_selected_time_period(self, file_option: str) -> bool:
        """
        Checks if the data file of the given option contains data collected within the selected start and end date.

        Args:
            file_option (str): Option name of a data file, which is a human-readable name for the file
        """
        collection_name = os.path.basename(self._collection_dir_path)
        if collection_name == self._elwha_collection_name:
            month_str, year_str = file_option.split()
            months = {
                "January": 1, "February": 2, "March": 3,
                "April": 4, "May": 5, "June": 6,
                "July": 7, "August": 8, "September": 9,
                "October": 10, "November": 11, "December": 12
            }
            date = dt.date(int(year_str), months[month_str], 1)
            start_date = dt.date(self.start_data_collection_date.year, self.start_data_collection_date.month, 1)
            end_date = dt.date(self.end_data_collection_date.year, self.end_data_collection_date.month, 1)
            return start_date <= date <= end_date
        return True

    def _save_changed_buffer_val(self, event: param.parameterized.Event) -> None:
        """
        Updates the buffers dictionary whenever any of the float input widgets (for each data file) change value.

        Args:
            event (param.parameterized.Event): Information about a buffer value change to a float input widget
        """
        self._update_buffer_config_file_button.disabled = False
        # Get the path to the data file that corresponds to the widget with the new buffer value.
        data_file_path = self._buffer_widget_file_path[event.obj.name]
        # Save the new buffer value.
        self._buffers[data_file_path] = event.new

    def _get_transect_search_radius_float_inputs(self) -> pn.Column:
        """
        Returns a list of float input widgets for each data file's buffer value.
        """
        widgets = []
        visited_subdirs = []
        for data_file_path, buffer in self._buffers.items():
            data_base_path, data_file = os.path.split(data_file_path)
            subdir_name = os.path.basename(data_base_path)
            # Create a heading for each subdirectory within the collection.
            if subdir_name not in visited_subdirs:
                visited_subdirs.append(subdir_name)
                # Get the actual title if a ScienceBase ID was used for the subdirectory name.
                if subdir_name in self._data_map.selected_collection_json_info:
                    subdir_name = self._data_map.selected_collection_json_info[subdir_name]
                subdir_heading = pn.pane.Markdown(
                    object = "**{}**".format(subdir_name),
                    sizing_mode = "stretch_width", margin = (10, 10, -15, 10)
                )
                widgets.append(subdir_heading)
            # Create float input widget.
            file_buffer_float_input = pn.widgets.FloatInput(
                name = data_file, value = buffer, start = 0, step = 1e-2,
                format = PrintfTickFormatter(format = "%.2f")
            )
            # Map the name of the float input widget to its data file path.
            self._buffer_widget_file_path[data_file] = data_file_path
            # Save new buffer values when a float input widget's value changes.
            file_buffer_float_input.param.watch(self._save_changed_buffer_val, "value")
            widgets.append(file_buffer_float_input)
        return widgets
    
    @param.depends("update_buffer_config", watch = True)
    def _update_buffer_config_file(self) -> None:
        """
        Saves the new buffer configurations to the JSON file storing buffer values used for the transect of each data file.
        """
        with pn.param.set_values(self._update_buffer_config_file_button, loading = True):
            with open(os.path.join(self._collection_dir_path, self._preprocessed_data_buffer_output), "w") as buffers_json_file:
                json.dump(self._buffers, buffers_json_file, indent = 4)
        self._update_buffer_config_file_button.disabled = True
    
    def _get_data_col_name(self, possible_data_cols: list[str]) -> str:
        """
        Gets the column name that exists in _data_map.all_data_cols, which is a list of column names provided by the user
        (any column in _data_map.all_data_cols could be used for the time-series plot's y-axis values).

        Args:
            possible_data_cols (list[str]): List of all column names in a time-series data file
        """
        for col in possible_data_cols:
            if col in self._data_map.all_data_cols: return col
        return self._default_y_axis_data_col_name
    
    def _get_coordinates_in_meters(self, x_coords: list[float], y_coords: list[float], crs: any) -> list[list[float]]:
        """
        Gets the given coordinates in meters, and transforms coordinates without a CRS into the time-series data's CRS.

        Args:
            x_coords (list[float]): List of x, longitude, or easting values (in degrees) to convert into meters
            y_coords (list[float]): List of y, latitude, or northing values (in degrees) to convert into meters
            crs (cartopy.crs or None): Source coordinate reference system of the given coordinates
        """
        if crs is None:
            crs = self._data_map.selected_collection_crs
            transformed_points = crs.transform_points(
                src_crs = self._data_map.map_default_crs,
                x = x_coords, y = y_coords
            )
            easting_vals = [point[0] for point in transformed_points]
            northing_vals = [point[1] for point in transformed_points]
            return [easting_vals, northing_vals]
        else:
            return [x_coords, y_coords]

    def _data_within_crs_bounds(self, x_data: list[float], y_data: list[float], crs: ccrs) -> bool:
        """
        Checks if the given data is within the bounds of the given CRS.

        Args:
            x_data (list[float]): List of x, longitude, or easting data values
            y_data (list[float]): List of y, latitude, or northing data values
            crs (cartopy.crs): Coordinate reference system used for checking if the data values lie within its bounds
        """
        if x_data and y_data:
            x_start, y_start, x_end, y_end = crs.boundary.bounds
            data_west_boundary, data_east_boundary = min(x_data), max(x_data)
            data_south_boundary, data_north_boundary = min(y_data), max(y_data)
            is_within_longitude_bounds = (x_start <= data_west_boundary <= data_east_boundary <= x_end)
            is_within_latitude_bounds = (y_start <= data_south_boundary <= data_north_boundary <= y_end)
            return (is_within_longitude_bounds and is_within_latitude_bounds)
        else:
            return False
    
    def _get_data_along_transect(self, data_file_path: str, transect_points: list[list[float]], long_col_name: str, lat_col_name: str, transect_crs: ccrs) -> pd.DataFrame:
        """
        Gets all data that was collected along the given transect and returns that data as a dataframe.
        Returns None if no data could be extracted with the given transect.

        Args:
            data_file_path (str): Path to the file containing data to extract for the time-series plot
            transect_points (list[list[float]]): List of coordinates for the transect's start (first item/list) and end (second item/list) points
                ^ [
                    [start point's longitude/easting, start point's latitude/northing],
                    [end point's longitude/easting, end point's latitude/northing]
                ]
            long_col_name (str): Name of the column containing the longitude/easting of each data point
            lat_col_name (str): Name of the column containing the latitude/northing of each data point
            transect_crs (cartopy.crs): Coordinate reference system of the given transect
        """
        _, data_file = os.path.split(data_file_path)
        _, extension = os.path.splitext(data_file)
        extension = extension.lower()
        if extension in [".tif", ".tiff"]:
            dataset = rxr.open_rasterio(data_file_path)
            # Clip data collected along the clicked transect from the given data file.
            try:
                transect_buffer = self._buffers.get(data_file_path, 0)
                if transect_buffer > 0:
                    padded_transect_polygon = LineString(transect_points).buffer(transect_buffer, cap_style = 2)
                    clipped_dataset = dataset.rio.clip(
                        geometries = [padded_transect_polygon],
                        from_disk = True
                    )
                else:
                    clipped_dataset = dataset.rio.clip(
                        geometries = [{
                            "type": "LineString",
                            "coordinates": transect_points
                        }],
                        from_disk = True
                    )
            except ValueError:
                # Given transect doesn't overlap data file, so return None early since the clipped dataset would be empty.
                return None
            # Convert clipped data into a GeoDataFrame to easily get each data point's distance from the transect's start point.
            clipped_dataset = clipped_dataset.squeeze().drop("spatial_ref").drop("band")
            # Set name of the column with time-series' y-axis values to the default value because ASCII grid files don't have data columns.
            self._y_axis_data_col_name = self._default_y_axis_data_col_name
            clipped_dataset.name = self._y_axis_data_col_name
            clipped_dataframe = clipped_dataset.to_dataframe().reset_index()
            no_data_val = clipped_dataset.attrs["_FillValue"]
            clipped_dataframe = clipped_dataframe[clipped_dataframe[self._y_axis_data_col_name] != no_data_val]
            clipped_geodataframe = gpd.GeoDataFrame(
                data = clipped_dataframe,
                geometry = gpd.points_from_xy(
                    x = clipped_dataframe["x"],
                    y = clipped_dataframe["y"],
                    crs = transect_crs
                )
            )
            # Calculate each point's distance from the transect's start point.
            transect_start_point = Point(transect_points[0])
            clipped_geodataframe[self._dist_col_name] = [point.distance(transect_start_point) for point in clipped_geodataframe.geometry]
            # Convert clipped data into a DataFrame for easier plotting.
            clipped_data_dataframe = clipped_geodataframe.drop(columns = "geometry").rename(
                columns = {
                    "x": long_col_name,
                    "y": lat_col_name
                }
            ).sort_values(by = self._dist_col_name).reset_index(drop = True)
            return clipped_data_dataframe
        elif extension == ".geojson":
            data_geodataframe = gpd.read_file(filename = data_file_path)
            # Reproject the data file to match the transect's projection, if necessary.
            if data_geodataframe.crs is None: data_geodataframe = data_geodataframe.set_crs(crs = self._data_map.map_default_crs)
            data_crs = ccrs.CRS(data_geodataframe.crs)
            if not data_crs.is_exact_same(transect_crs): data_geodataframe = data_geodataframe.to_crs(crs = transect_crs)
            # Add buffer/padding to the clicked transect, which is created with the given transect's start and end point coordinates.
            # ^ Buffer allows data points within a certain distance from the clicked transect to be included in the time-series (since it's rare for data points to lie exactly on a transect).
            padded_transect = LineString(transect_points).buffer(self._buffers.get(data_file_path, 3), cap_style = 2)
            # Create GeoDataFrame for the padded transect.
            clicked_transect_geodataframe = gpd.GeoDataFrame(
                data = {"geometry": [padded_transect]},
                geometry = "geometry",
                crs = transect_crs
            )
            # Clip data collected along the clicked transect from the given data file.
            clipped_geodataframe = data_geodataframe.clip(mask = clicked_transect_geodataframe)
            # Given transect doesn't overlap data file, so return None early since the clipped geodataframe would be empty.
            if clipped_geodataframe.empty: return None
            # Calculate each point's distance from the transect's start point.
            transect_start_point = Point(transect_points[0])
            clipped_geodataframe.insert(
                loc = len(clipped_geodataframe.columns),
                column = self._dist_col_name,
                value = [point.distance(transect_start_point) for point in clipped_geodataframe.geometry]
            )
            # Convert clipped data into a DataFrame for easier plotting.
            clipped_data_dataframe = clipped_geodataframe.drop(columns = "geometry").sort_values(by = self._dist_col_name).reset_index(drop = True)
            # Get name of the column with time-series' y-axis values.
            self._y_axis_data_col_name = self._get_data_col_name(list(clipped_data_dataframe.columns))
            return clipped_data_dataframe
        elif extension in [".parq", ".parquet"]:
            data_dask_geodataframe = dask_geopandas.read_parquet(data_file_path).to_crs(self._data_map.selected_collection_json_info.get("epsg", 4326))
            # Add buffer/padding to the clicked transect, which is created with the given transect's start and end point coordinates.
            # ^ Buffer allows data points within a certain distance from the clicked transect to be included in the time-series (since it's rare for data points to lie exactly on a transect).
            padded_transect = LineString(transect_points).buffer(self._buffers.get(data_file_path, 3), cap_style = 2)
            # Create Dask GeoDataFrame for the padded transect.
            clicked_transect_dask_geodataframe = dask_geopandas.from_geopandas(
                gpd.GeoDataFrame(
                    data = {"geometry": [padded_transect]},
                    geometry = "geometry",
                    crs = transect_crs
                ),
                npartitions = 1
            )
            # Clip data collected along the clicked transect from the given data file.
            clipped_geodataframe = data_dask_geodataframe.sjoin(
                df = clicked_transect_dask_geodataframe,
                how = "inner", predicate = "intersects"
            )
            # Calculate each point's distance from the transect's start point.
            transect_start_point = Point(transect_points[0])
            clipped_data_dataframe = clipped_geodataframe.assign(
                new_dist_col_name = clipped_geodataframe["geometry"].distance(transect_start_point)
            ).rename(columns = {"new_dist_col_name": self._dist_col_name}).sort_values(by = self._dist_col_name).drop(columns = ["index_right", "geometry"]).compute()
            # When clipped dataframe is empty, then return None early.
            if clipped_data_dataframe.empty:
                return None
            else:
                # Get name of the column with time-series' y-axis values.
                self._y_axis_data_col_name = self._get_data_col_name(list(clipped_data_dataframe.columns))
                return clipped_data_dataframe
        # Return None if there's currently no implementation to extract data from the data file yet.
        print("Error extracting data along a transect from", data_file, ":", "Files with the", extension, "file format are not supported yet.")
        return None

    async def _clip_data(self, file_path: str, easting_data: list[float], northing_data: list[float], long_col_name: str, lat_col_name: str, transect_crs: ccrs) -> hv.Overlay:
        """
        Clips data from the given file path with the selected transect.

        Args:
            file_path (str): Path to the data file, which is used to extract data for the time-series
            easting_data (list[float]): List of longitude/easting values (in meters) for each transect point
            northing_data (list[float]): List of latitude/northing values (in meters) for each transect point
            long_col_name (str): Name of the column containing the longitude/easting of each data point
            lat_col_name (str): Name of the column containing the latitude/northing of each data point
            transect_crs (cartopy.crs): Coordinate reference system of the given transect
        """
        start_time = time.time()
        subdir_path, filename = os.path.split(file_path)
        subdir = os.path.basename(subdir_path)
        file_option = " - ".join([subdir, filename])
        if file_path in self._data_map.selected_collection_json_info:
            file_option = self._data_map.selected_collection_json_info[file_path]
        # Clip data along the selected transect for each data file.
        clipped_dataframe = self._get_data_along_transect(
            data_file_path = file_path,
            transect_points = list(zip(easting_data, northing_data, strict = True)),
            long_col_name = long_col_name,
            lat_col_name = lat_col_name,
            transect_crs = transect_crs
        )
        if clipped_dataframe is not None:
            # Assign the time-series plot's options.
            x_axis_col = self._dist_col_name
            y_axis_col = self._y_axis_data_col_name
            other_val_cols = [col for col in clipped_dataframe.columns if col not in [x_axis_col, y_axis_col]]
            # Plot clipped data.
            clipped_data_curve_plot = hv.Curve(
                data = clipped_dataframe,
                kdims = x_axis_col,
                vdims = y_axis_col,
                label = file_option
            ).opts(color = self._data_file_colors[file_path])
            clipped_data_point_plot = hv.Points(
                data = clipped_dataframe,
                kdims = [x_axis_col, y_axis_col],
                vdims = other_val_cols,
                label = file_option
            ).opts(
                color = self._data_file_colors[file_path],
                tools = ["hover"],
                size = 5
            )
            # Save the time-series data for the given file as a pandas DataFrame.
            all_dims = clipped_data_curve_plot.dimensions(selection = "all")
            data_col_name = "{}: {}".format(file_option, y_axis_col)
            dataframe = clipped_data_curve_plot.dframe(dimensions = all_dims).rename(columns = {y_axis_col: data_col_name})
            self._time_series_dataframes.append(dataframe)
            # Return the clipped data file's plot.
            clipped_data_plot = clipped_data_curve_plot * clipped_data_point_plot
            end_time = time.time()
            print("Creating time-series plot for {} took {} seconds.".format(file_path, end_time - start_time))
            return clipped_data_plot
        else:
            return None

    async def _create_time_series_plot(self) -> pn.pane.HoloViews:
        """
        Creates and returns a time-series plot for data collected along a clicked transect on the map.
        """
        # Get informational key-value pairs that aren't part of the time-series plot.
        transect_file = self.clicked_transects_info.get(self._clicked_transects_file, None)
        num_transects = self.clicked_transects_info.get(self._num_clicked_transects, 0)
        transect_crs = self.clicked_transects_info.get(self._clicked_transects_crs, self._data_map.selected_collection_crs)
        long_col_name = self.clicked_transects_info.get(self._clicked_transects_longitude_col, "Longitude")
        lat_col_name = self.clicked_transects_info.get(self._clicked_transects_latitude_col, "Latitude")
        transect_id_col_name = self.clicked_transects_info.get(self._clicked_transects_id_col, "Transect ID")
        if num_transects == 1:
            # Get the ID of the selected transect without the set's brackets.
            transect_id = str(set(self.clicked_transects_info[transect_id_col_name]))[1:-1]
            # Transform any user-drawn transect's coordinates into a CRS with meters as a unit.
            easting_data, northing_data = [], []
            if (long_col_name in self.clicked_transects_info) and (lat_col_name in self.clicked_transects_info):
                easting_data, northing_data = self._get_coordinates_in_meters(
                    x_coords = self.clicked_transects_info[long_col_name],
                    y_coords = self.clicked_transects_info[lat_col_name],
                    crs = transect_crs if self._clicked_transects_crs in self.clicked_transects_info else None
                )
            # For each data file, plot its data collected along the clicked transect.
            plot = None
            if self._data_within_crs_bounds(x_data = easting_data, y_data = northing_data, crs = transect_crs):
                # Create a list of tasks (clip all selected data files with the selected transect) to run asynchronously.
                self._time_series_dataframes = []
                tasks = [asyncio.create_task(self._clip_data(file_path, easting_data, northing_data, long_col_name, lat_col_name, transect_crs)) for file_path in self._user_selected_data_files]
                # Gather the returned results of each task.
                start_time = time.time()
                results = await asyncio.gather(*tasks)
                end_time = time.time()
                print("Creating all time-series plots took {} seconds.".format(end_time - start_time))
                # Overlay all clipped data files' plots.
                start_time = time.time()
                for clipped_data_plot in results:
                    if clipped_data_plot is not None:
                        if plot is None:
                            self._update_heading_text(
                                title = "Time-Series of Data Collected Along Transect {} from {}".format(
                                    transect_id, transect_file
                                ),
                                details = "Scroll on the axes or data area to zoom in and out of the plot."
                            )
                            plot = clipped_data_plot
                        else:
                            plot = plot * clipped_data_plot
                end_time = time.time()
                print("Overlaying all time-series plots took {} seconds.".format(end_time - start_time))
            if plot is not None:
                # Return the overlay plot containing data collected along the transect for all data files.
                self._time_series_plot = pn.pane.HoloViews(
                    object = plot.opts(
                        title = "Time-Series of {} Data".format(self.data_category),
                        xlabel = self._dist_col_name,
                        ylabel = self._y_axis_data_col_name,
                        active_tools = ["pan", "wheel_zoom"],
                        legend_position = "bottom", legend_cols = 3,# or any integer can eventually be used in Bokeh 3.1.0 and Holoviews 1.16.0 so (TODO) make sure to update the Bokeh package when it becomes compatible with HoloViews!
                        legend_limit = len(self._user_selected_data_files) * 2,
                        show_legend = True, toolbar = None,
                        height = 500, responsive = True, padding = 0.1
                    ),
                    visible = True
                )
            else:
                self._update_heading_text(
                    title = "No Time-Series Available",
                    details = "Unfortunately, no data has been collected along your selected transect (Transect {} from {}). Please select another transect or create your own transect.".format(
                        transect_id, transect_file
                    )
                )
                self._time_series_plot = pn.pane.HoloViews(object = None, visible = False)
        elif num_transects > 1:
            # Make the selected transects' IDs readable for the error message.
            ids = sorted([str(id) for id in set(self.clicked_transects_info[transect_id_col_name])])
            transect_ids = ""
            if len(ids) == 2: transect_ids = "{} and {}".format(*ids)
            else: transect_ids = ", ".join(ids[:-1]) + ", and " + str(ids[-1])
            self._update_heading_text(
                title = "More than One Selected Transect",
                details = "Transects with IDs {} from {} have been selected. Please select only one transect because the time-series of data can (currently) only be generated from one transect.".format(
                    transect_ids, transect_file
                )
            )
            self._time_series_plot = pn.pane.HoloViews(object = None, visible = False)
        # Set the visibility of the time-series download button based on whether a time-series is computed.
        self._time_series_download_button.visible = self._time_series_plot.visible
        # Return the newly computed time-series plot.
        return self._time_series_plot

    @param.depends("download_time_series", watch = True)
    def _download_time_series(self) -> None:
        """
        Saves different versions of the time-series plot.
        """
        current_time = dt.datetime.now()
        # Create directory to hold all downloaded files.
        downloads_dir_name = current_time.strftime("%b-%d-%Y_%Hhr%Mmin%Ssec")
        downloads_dir_path = os.path.join(self._outputs_dir_path, downloads_dir_name)
        if not os.path.exists(downloads_dir_path): os.makedirs(downloads_dir_path)
        # Create downloaded files's name.
        category_name = self.data_category.replace(" ", "_")
        if "(" in category_name: category_name = category_name.split("(")[1].replace(")", "")
        transect_name = self._modal_heading.objects[0].object.split(" Transect ")[1].replace(" ", "_").replace(".", "_")
        filename = "{}_Time_Series_Along_Transect_{}".format(category_name, transect_name)
        # Save HTML version.
        html_name = filename + ".html"
        html_path = os.path.join(downloads_dir_path, html_name)
        html_content = pn.Column(self._time_series_plot, "Selected Transect(s) Data", self._clicked_transects_table)
        html_content.save(filename = html_path)
        # Save CSV version (combine all time-series dataframes together, group the rows with the same distance together, and round the precision of the distances to at most 2 decimal places).
        csv_name = filename + ".csv"
        csv_path = os.path.join(downloads_dir_path, csv_name)
        all_time_series_data = pd.concat(objs = self._time_series_dataframes, axis = 0, ignore_index = True).sort_values(by = self._dist_col_name).reset_index(drop = True)
        all_time_series_data = all_time_series_data.groupby(by = self._dist_col_name, as_index = False).aggregate("first")
        all_time_series_data[self._dist_col_name] = all_time_series_data[self._dist_col_name].apply(lambda dist_val: round(number = dist_val, ndigits = 2))
        with open(csv_path, "a+") as csv_file:
            # Use \" (escape character followed by the double-quote character) at the start and end of each comment to prevent commas from breaking sentences into several columns in Excel.
            # ^ Make sure to read the downloaded CSV file with pandas.read_csv(filepath_or_buffer = ..., comment = "\"").
            csv_file.write("\"# ============================================================\"\n")
            csv_file.write("\"# The data contained in this comma separated value file is a subset of data that originally come from the U.S. Geological Survey (USGS), published as: Stevens, A.W., Gelfenbaum, G., Warrick, J.A., Miller, I.M., and Weiner, H.M., 2017, Bathymetry, topography, and sediment grain-size data from the Elwha River delta, Washington: U.S. Geological Survey data release, https://doi.org/10.5066/F72N51GC.\"\n")
            csv_file.write("\"# Unless otherwise stated, all data, metadata and related materials are considered to satisfy the quality standards relative to the purpose for which the data were collected. Although these data and associated metadata have been reviewed for accuracy and completeness and approved for release by the U.S. Geological Survey (USGS), no warranty expressed or implied is made regarding the display or utility of the data for other purposes, nor on all computer systems, nor shall the act of distribution constitute any such warranty.\"\n")
            csv_file.write("\"# The data contained in this comma separated value file may have been modified minimally by the software, Elwha Topo-Bathy Data Viewer, that generated this comma separated value file. This version of the data is provided to meet the need for timely best science.\"\n")
            csv_file.write("\"# This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.\"\n")
            csv_file.write("\"# ============================================================\"\n")
            all_time_series_data.to_csv(path_or_buf = csv_file, sep = ",", na_rep = np.nan, lineterminator = "\n", index = False)
        # Save the buffer configurations used for creating the time-series.
        with open(os.path.join(downloads_dir_path, "buffer_config.json"), "w") as buffer_json_file:
            json.dump(self._buffers, buffer_json_file, indent = 4)

    def _update_clicked_transects_table(self) -> pn.widgets.DataFrame:
        """
        Updates and returns the Panel DataFrame widget with new information about the newly clicked transect(s).
        """
        # Create a new dictionary containing the new dataframe's columns.
        dataframe_cols = self.clicked_transects_info.get(self._clicked_transects_table_cols, [])
        new_transects_data = {col: vals for col, vals in self.clicked_transects_info.items() if col in dataframe_cols}
        # Add a new "Point Type" column to differentiate each transect's points.
        point_type_col_vals = []
        if self.clicked_transects_info:
            transect_id_col_name = self.clicked_transects_info.get(self._clicked_transects_id_col, "Transect ID")
            transect_id_col_vals, seen_ids = self.clicked_transects_info[transect_id_col_name], set()
            transect_id_to_num_points_dict = dict(Counter(transect_id_col_vals))
            for id in transect_id_col_vals:
                if id not in seen_ids:
                    seen_ids.add(id)
                    transect_point_type_vals = ["middle"] * transect_id_to_num_points_dict[id]
                    transect_point_type_vals[0] = "start"
                    transect_point_type_vals[-1] = "end"
                    point_type_col_vals = point_type_col_vals + transect_point_type_vals
        new_transects_data[self._point_type_col_name] = point_type_col_vals
        # Create a new dataframe with information about the newly clicked transect(s).
        new_transects_dataframe = pd.DataFrame.from_dict(new_transects_data).set_index(self._point_type_col_name)
        # Update the dataframe for the transects table.
        self._clicked_transects_table.value = new_transects_dataframe
        return self._clicked_transects_table

    def _update_heading_text(self, title: str = "Loading...", details: str = "") -> None:
        """
        Updates the heading text at the top of the popup modal by setting new values for the rendered Panel markdown objects.

        Args:
            title (str): Heading of modal
            details (str): Subheading of modal
        """
        title_markdown = "#### {}".format(title)
        details_markdown = "##### {}".format(details)
        self._modal_heading.objects[0].object = title_markdown
        self._modal_heading.objects[1].object = details_markdown
    
    @param.depends("data_category", "start_data_collection_date", "end_data_collection_date", watch = True)
    def _update_selected_data_files(self) -> None:
        """
        Selects data files for the time-series based on the selected data category and time period.
        """
        new_selected_data_files_paths = []
        for category in self._data_category_select.options:
            if (self.data_category is not None) and (category == self.data_category):
                if category != self._placeholder_data_category:
                    # Select data files that belong to the selected data category and lie within the selected time period.
                    valid_category_data_files_paths = []
                    for file_path in self._data_map.selected_collection_json_info["categories"][category]:
                        file_option = self._data_map.selected_collection_json_info[file_path]
                        if self._within_selected_time_period(file_option): valid_category_data_files_paths.append(file_path)
                    # Automatically select all valid data files in the widget that corresponds to their data category.
                    if category in self._category_multiselect_widget: self._category_multiselect_widget[category].value = valid_category_data_files_paths
                    new_selected_data_files_paths.extend(valid_category_data_files_paths)
            else:
                # Unselect data files that do not belong to the selected data category.
                if category in self._category_multiselect_widget: self._category_multiselect_widget[category].value = []
        # Assign a color for each selected data file that might appear in the time-series plot.
        self._data_file_colors = {}
        indices = np.round(np.linspace(0, len(Turbo256) - 1, len(new_selected_data_files_paths))).astype(int)
        for path_idx, data_file_path in enumerate(new_selected_data_files_paths):
            color_idx = indices[path_idx]
            self._data_file_colors[data_file_path] = Turbo256[color_idx]
        # Set new options for the widget that lets the user specify which selected data file to display.
        new_selected_data_options = {self._placeholder_displayed_data: None}
        for data_file_path in new_selected_data_files_paths:
            option_name = self._data_map.selected_collection_json_info[data_file_path]
            new_selected_data_options[option_name] = data_file_path
        self._displayed_data_select.options = new_selected_data_options
        # Set the default data file to display.
        if new_selected_data_files_paths and (self.data_category is not None) and (self.data_category != self._placeholder_data_category):
            self._displayed_data_select.value = new_selected_data_files_paths[-1]
        else:
            self._displayed_data_select.value = None
        # Save the newly selected data files.
        self._user_selected_data_files = new_selected_data_files_paths
    
    @param.depends("update_collection_dir_path", watch = True)
    def _update_collection_objects(self) -> None:
        """
        Assign styles for each time-series data file in the newly selected collection directory.
        """
        self._collection_dir_path = self._data_map.selected_collection_dir_path
        # Update widgets in the "Time-Series Data" section.
        data_categories = list(self._data_map.selected_collection_json_info["categories"].keys())     # name of the key should be same as `collection_data_categories_property` in utils/preprocess_data.py
        self._data_category_select.options = [self._placeholder_data_category] + data_categories
        self._categorize_data_files()
        self._update_selected_data_files()
        # Load buffer configuration file's values.
        json_file = open(os.path.join(self._collection_dir_path, self._preprocessed_data_buffer_output))
        self._buffers = json.load(json_file)
        # Update widgets in the "Transect Search Radius" section.
        self._transect_search_radius_widgets.objects = self._transect_search_radius_constant_widgets + self._get_transect_search_radius_float_inputs()
    
    # -------------------------------------------------- Public Class Properties & Methods --------------------------------------------------
    @param.depends("clicked_transects_info")
    def content(self) -> pn.Column:
        """
        Returns a Panel column with components to display in the popup modal whenever a new transect is selected.
        """
        # Update the modal heading to tell the user that its contents are being computed.
        self._update_heading_text()
        # Return the new modal contents.
        return pn.Column(
            objects = [
                *(self._modal_heading),
                pn.panel(self._create_time_series_plot, loading_indicator = True),
                pn.Column(
                    "Selected Transect(s) Data",
                    pn.panel(self._update_clicked_transects_table, loading_indicator = True),
                    sizing_mode = "stretch_width"
                ),
                self._time_series_download_button
            ],
            sizing_mode = "stretch_width"
        )
    
    def get_sidebar_widgets(self) -> list:
        """
        Returns a list of widgets that appear in the sidebar of the app.
        """
        collection_name = os.path.basename(self._collection_dir_path)
        if collection_name == self._elwha_collection_name:
            return [
                self._start_data_collection_date_picker,
                self._end_data_collection_date_picker,
                self._data_category_select
            ]
        else:
            return []

    def get_accordion_sections(self) -> list[tuple]:
        """
        Returns a list of tuples, each containing the name of the accordion section and its content.
        """
        collection_name = os.path.basename(self._collection_dir_path)
        if collection_name == self._elwha_collection_name:
            return [
                ("Time-Series Data", self._data_files_widgets),
                ("Transect Search Radius", self._transect_search_radius_widgets)
            ]
        else:
            return []