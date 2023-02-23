# Standard library imports
import os
import json
from collections import defaultdict, Counter
import asyncio

# External dependencies imports
import param
import panel as pn
import holoviews as hv
import rioxarray as rxr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString
import cartopy.crs as ccrs
from bokeh.models.formatters import PrintfTickFormatter
from .DataMap import DataMap

### PopupModal is used to display a time-series plot or any other data/message in the app's modal. ###
class PopupModal(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    start_data_collection_date = param.Date(label = "Start Datetime of Time-Series Data")
    end_data_collection_date = param.Date(label = "End Datetime of Time-Series Data")
    update_collection_dir_path = param.Event(label = "Action that Triggers the Updating of the Collection Directory and Its Related Objects")
    user_selected_data_files = param.ListSelector(default = [], label = "Time-Series Data Files")
    update_buffer_config = param.Event(label = "Action that Triggers Updating the Buffer Config File")
    update_modal = param.Event(label = "Action that Triggers the Updating of Modal Contents")

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, data_map: DataMap, template: pn.template, time_series_data_col_names: list[str] = [], **params) -> None:
        """
        Creates a new instance of the PopupModal class with its instance variables.

        Args:
            data_map (DataMap): Instance containing methods for converting data files to allow quicker loading onto a map
            template (panel.template): Data visualizer app's template
            time_series_data_col_names (list[str]): Optional list of column names for columns containing data for the time-series' y-axis
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._data_map = data_map
        self._app_template = template
        self._all_data_cols = time_series_data_col_names
        
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

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _modal_heading_pipe = pipe stream that contains text for the popup modal's heading
        # ^ first string in list is the title of the modal
        # ^ second string in list is more detail about the modal's contents
        self._modal_heading_pipe = hv.streams.Pipe(data = [])
        # Update the _modal_heading property whenever data in the modal heading pipe changes.
        self._modal_heading_pipe.add_subscriber(self._update_heading_text)
        
        # _collection_dir_path = path to the directory containing all the data files used for the time-series
        self._collection_dir_path = data_map.selected_collection_dir_path
        # _clicked_transects_pipe = pipe stream that contains info about the most recently clicked transect(s)
        self._clicked_transects_pipe = hv.streams.Pipe(data = {})
        # _y_axis_data_col_name = name of the column that stores the y-axis values for the time-series plot
        self._y_axis_data_col_name = self._default_y_axis_data_col_name
        # _time_series_plot = time-series plot for data collected along the most recently clicked transect/path
        self._time_series_plot = hv.DynamicMap(self._create_time_series_plot, streams = [self._clicked_transects_pipe]).opts(
            title = "Time-Series",
            xlabel = self._dist_col_name,
            hooks = [self._update_modal_content],
            active_tools = ["pan", "wheel_zoom"],
            show_legend = True, toolbar = None,
            height = 500, responsive = True, padding = 0.1
        )
        # _selected_file_groups = set containing unique group names that each selected data file belongs to or is similar to
        # ^ data files that use the same/similar measurement for the time-series' y-axis values are in the same group
        self._selected_file_groups = set()
        # _checkbox_group_widget = dictionary mapping each data file group name (key) to its checkbox widget (value)
        self._checkbox_group_widget = {}
        # _buffers = dictionary mapping each data file's path (key) to the selected transect's buffer/search radius (value) when extracting data around the transect
        self._buffers = {}
        # _buffer_widget_file_path = dictionary mapping the name of each float input widget (key) to the path (value) of the data file that uses this buffer when extracting data around the transect
        self._buffer_widget_file_path = {}

        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # _data_files_widgets = column layout containing widgets for the "Time-Series Data" accordion section
        self._data_files_widgets = pn.Column(objects = [])
        # _update_buffer_config_file_button = button for updating the collection's buffer config file with values from the buffer float widgets for each data file
        self._update_buffer_config_file_button = pn.widgets.Button.from_param(
            parameter = self.param.update_buffer_config,
            name = "Save Current Values to {}".format(self._preprocessed_data_buffer_output),
            button_type = "primary", disabled = True
        )
        # _wiki_info_button = button that opens a tab to the GitHub Wiki page that describes how to use the time-series controls
        self._wiki_info_button = pn.widgets.Button(name = "\u2139", button_type = "light", width = 30)
        self._wiki_info_button.js_on_click(
            args = {"wiki_url": "https://github.com/venuswku/data-visualizer/wiki/Features#visualize-data"},
            code = "window.open(wiki_url)"
        )
        # _time_series_data_constant_widgets = list of widgets that always appear at the top of the "Time-Series Data" accordion section
        self._time_series_data_constant_widgets = [
            pn.Row(
                pn.widgets.StaticText(value = "Select data files to use when creating a time-series of how data changes over time along a chosen transect."),
                # self._wiki_info_button
            )
        ]
        # _transect_search_radius_constant_widgets = list of widgets that always appear at the top of the "Transect Search Radius" accordion section
        self._transect_search_radius_constant_widgets = [
            pn.Row(
                pn.widgets.StaticText(value = "Adjust the search radius for extracting time-series data around a selected transect."),
                # self._wiki_info_button
            ),
            self._update_buffer_config_file_button
        ]
        # _transect_search_radius_widgets = column layout containing widgets for the "Transect Search Radius" accordion section
        self._transect_search_radius_widgets = pn.Column(objects = [])
        # _time_series_controls_accordion = accordion layout widget allowing the user to change settings for the time-series
        self._time_series_controls_accordion = pn.Accordion(
            objects = [
                ("Time-Series Data", self._data_files_widgets),
                ("Transect Search Radius", self._transect_search_radius_widgets)
            ],
            active = [], toggle = True, sizing_mode = "stretch_width"
        )

        # _modal_heading = list containing markdown objects for the modal heading (title for first markdown, details for second markdown)
        self._modal_heading = [pn.pane.Markdown(object = ""), pn.pane.Markdown(object = "", margin = (-20, 5, 0, 5))]
        # _clicked_transects_table = table containing information about the clicked transect(s)'s start and end points
        self._clicked_transects_table = pn.widgets.DataFrame(
            value = pd.DataFrame(),
            name = "Selected Transect(s) Data",
            show_index = True, auto_edit = False, text_align = "center",
            sizing_mode = "stretch_width", margin = (-20, 5, 10, 5)
        )

        # Initialize widgets that depend on the selected collection from DataMap.
        self._update_collection_objects()

    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    def _save_selected_data_files(self, event: param.parameterized.Event) -> None:
        """
        Check if the most recently selected data file has a measurement (for the time-series' y-axis) that is compatible with other selected data files.
        If the data file has a matching measurement, then save the file's path to the user_selected_data_files parameter; else display an error message.

        Args:
            event (param.parameterized.Event): Event caused by a value change to one of the checkbox group widgets
        """
        collection_name = os.path.basename(self._collection_dir_path)
        if collection_name == "5a01f6d0e4b0531197b72cfe":
            elevation_groups = ["Digital Elevation Models (DEMs)", "Bathymetry (Kayak)", "Bathymetry (Personal Watercraft)", "Topography"]
            f_w_mean_group = "Surface-Sediment Grain-Size Distributions"
            if (event.old is None) or ((event.old is not None) and (len(event.new) > len(event.old))):
                # Check if the newly selected data file is in the same or similar group as the selected ones.
                newly_selected_file_group = event.obj.name
                newly_selected_file_path = event.new[-1]
                newly_selected_file_name = os.path.basename(newly_selected_file_path)
                if self._selected_file_groups:
                    # When a grain-size data file is recently selected but data files from any of the elevation groups were already selected...
                    if (newly_selected_file_group == f_w_mean_group) and (f_w_mean_group not in self._selected_file_groups):
                        self._data_map.error_messages.append(" ".join([
                            "You can only select data files with matching measurements for the time-series.",
                            "{}'s `F-W Mean` measurements are not compatible with other selected data's `Elevation` measurements.".format(newly_selected_file_name),
                            "Please either unselect all the currently selected data file(s) in order to select {} or continue selecting data files under any of the following sections: {}.".format(newly_selected_file_name, ", ".join(elevation_groups))
                        ]))
                        self._checkbox_group_widget[newly_selected_file_group].value = []
                        self._data_map.data_file_path = None
                    # When a data file from one of the elevation groups is recently selected but one or more grain-size data files was already selected...
                    elif (newly_selected_file_group in elevation_groups) and (f_w_mean_group in self._selected_file_groups):
                        self._data_map.error_messages.append(" ".join([
                            "You can only select data files with matching measurements for the time-series.",
                            "{}'s `Elevation` measurements are not compatible with other selected data's `F-W Mean` measurements.".format(newly_selected_file_name),
                            "Please either unselect all the currently selected data file(s) in order to select {} or continue selecting data files under {}.".format(newly_selected_file_name, f_w_mean_group)
                        ]))
                        self._checkbox_group_widget[newly_selected_file_group].value = []
                        self._data_map.data_file_path = None
                    # When the recently selected data file belongs to a group that is compatible with the already selected data files' groups...
                    else:
                        if newly_selected_file_group in elevation_groups: self._selected_file_groups.update(elevation_groups)
                        elif newly_selected_file_group == f_w_mean_group: self._selected_file_groups.add(f_w_mean_group)
                        self.user_selected_data_files = self.user_selected_data_files + [newly_selected_file_path]
                else:
                    # Add the data file if there are no selected files.
                    if newly_selected_file_group in elevation_groups: self._selected_file_groups.update(elevation_groups)
                    elif newly_selected_file_group == f_w_mean_group: self._selected_file_groups.add(f_w_mean_group)
                    self.user_selected_data_files = self.user_selected_data_files + [newly_selected_file_path]
            else:
                # A data file was unselected, so remove it from the user_selected_data_files parameter.
                newly_unselected_file_group = event.obj.name
                if newly_unselected_file_group in elevation_groups: self._selected_file_groups.difference_update(set(elevation_groups))
                elif newly_unselected_file_group == f_w_mean_group: self._selected_file_groups.discard(f_w_mean_group)
                newly_unselected_file_path = event.old[-1]
                self.user_selected_data_files = [path for path in self.user_selected_data_files if path != newly_unselected_file_path]
        else:
            # TODO: finish
            print("either add or remove the new data file")

    def _group_data_files(self, all_data_files: dict) -> None:
        """
        Assigns the given data files into groups of checkboxes.

        Args:
            all_data_files (dict): Dictionary mapping each data file's name (key) to its path/location (value) on your local machine
        """
        widgets = self._time_series_data_constant_widgets
        collection_name = os.path.basename(self._collection_dir_path)
        # Group data files by the type of data for the Elwha collection.
        if collection_name == "5a01f6d0e4b0531197b72cfe":
            checkbox_group_options = defaultdict(lambda: {})
            # Group data files.
            for filename, file_path in all_data_files.items():
                if "_dem_" in filename: checkbox_group_options["Digital Elevation Models (DEMs)"][filename] = file_path
                elif "_kayak" in filename: checkbox_group_options["Bathymetry (Kayak)"][filename] = file_path
                elif ("_pwc" in filename) or ("_bathy" in filename): checkbox_group_options["Bathymetry (Personal Watercraft)"][filename] = file_path
                elif "_topo" in filename: checkbox_group_options["Topography"][filename] = file_path
                elif "_grainsize" in filename: checkbox_group_options["Surface-Sediment Grain-Size Distributions"][filename] = file_path
                else: checkbox_group_options["Other"][filename] = file_path
            # Create a checkbox group widget for each group.
            for group_name, options_dict in checkbox_group_options.items():
                if options_dict:
                    group_heading = pn.pane.Markdown(object = "**{}**".format(group_name), sizing_mode = "stretch_width", margin = (10, 10, -10, 10))
                    widgets.append(group_heading)
                    checkbox_group = pn.widgets.CheckBoxGroup(name = group_name, options = options_dict, value = [])
                    widgets.append(checkbox_group)
                    self._checkbox_group_widget[group_name] = checkbox_group
                    # Save newly selected data files when a checkbox group widget's value changes.
                    checkbox_group.param.watch(self._save_selected_data_files, "value")
        else:
            single_checkbox_group = pn.widgets.CheckBoxGroup(name = "Other", options = all_data_files, value = [])
            self._checkbox_group_widget["Other"] = single_checkbox_group
            # Save newly selected data files when the checkbox group widget's value changes.
            single_checkbox_group.param.watch(self._save_selected_data_files, "value")
            widgets.append(single_checkbox_group)
        # Assign new widgets for allowing the user to choose time-series data files.
        self._data_files_widgets.objects = widgets

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
                    sizing_mode = "stretch_width", margin = (10, 10, -10, 10)
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
        Gets the column name that exists in _all_data_cols, which is a list of column names provided by the user
        (any column in _all_data_cols could be used for the time-series plot's y-axis values).

        Args:
            possible_data_cols (list[str]): List of all column names in a time-series data file
        """
        for col in possible_data_cols:
            if col in self._all_data_cols: return col
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
                    padded_transect_polygon = LineString(transect_points).buffer(transect_buffer)
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
            padded_transect = LineString(transect_points).buffer(self._buffers.get(data_file_path, 3))
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
        subdir_path, filename = os.path.split(file_path)
        subdir = os.path.basename(subdir_path)
        file_option = ": ".join([subdir, filename])
        if file_option in self._data_map.selected_collection_json_info:
            file_option = self._data_map.selected_collection_json_info[file_option]
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
            ).opts(
                color = self._data_map.data_file_color[file_path],
                line_dash = self._data_map.data_file_line_style[file_path]
            )
            clipped_data_point_plot = hv.Points(
                data = clipped_dataframe,
                kdims = [x_axis_col, y_axis_col],
                vdims = other_val_cols,
                label = file_option
            ).opts(
                color = self._data_map.data_file_color[file_path],
                marker = self._data_map.data_file_marker[file_path],
                tools = ["hover"],
                size = 10
            )
            # Return the clipped data file's plot.
            clipped_data_plot = clipped_data_curve_plot * clipped_data_point_plot
            return clipped_data_plot
    
    # async def _get_clipped_data_plots(self, easting_data: list[float], northing_data: list[float], long_col_name: str, lat_col_name: str, transect_crs: ccrs) -> hv.Overlay:
    #     """
    #     Asynchronously clips all selected data files with the selected transect.

    #     Args:
    #         easting_data (list[float]): List of longitude/easting values (in meters) for each transect point
    #         northing_data (list[float]): List of latitude/northing values (in meters) for each transect point
    #         long_col_name (str): Name of the column containing the longitude/easting of each data point
    #         lat_col_name (str): Name of the column containing the latitude/northing of each data point
    #         transect_crs (cartopy.crs): Coordinate reference system of the given transect
    #     """
    #     # Create a list of tasks to run asynchronously.
    #     tasks = [asyncio.create_task(self._clip_data(file_path, easting_data, northing_data, long_col_name, lat_col_name, transect_crs)) for file_path in self.user_selected_data_files]
    #     # Gather the returned results of each task.
    #     results = await asyncio.gather(*tasks)
    #     print(results)
    #     # Overlay all clipped data files' plots.
    #     plot = None
    #     for clipped_data_plot in results:
    #         if plot is None: plot = clipped_data_plot
    #         else: plot = plot * clipped_data_plot
    #     return plot
    
    async def _create_time_series_plot(self, data: dict = {}) -> hv.Overlay:
        """
        Creates a time-series plot for data collected along a clicked transect on the map.

        Args:
            data (dict): Dictionary containing information about the clicked transect(s)
                and maps each data column (keys) to a list of values for that column (values)
        """
        # Get informational key-value pairs that aren't part of the time-series plot.
        transect_file = data.get(self._clicked_transects_file, None)
        num_transects = data.get(self._num_clicked_transects, 0)
        transect_crs = data.get(self._clicked_transects_crs, self._data_map.selected_collection_crs)
        long_col_name = data.get(self._clicked_transects_longitude_col, "Longitude")
        lat_col_name = data.get(self._clicked_transects_latitude_col, "Latitude")
        transect_id_col_name = data.get(self._clicked_transects_id_col, "Transect ID")
        self._update_clicked_transects_table(info = data)
        if num_transects == 1:
            # Get the ID of the selected transect without the set's brackets.
            transect_id = str(set(data[transect_id_col_name]))[1:-1]
            # Transform any user-drawn transect's coordinates into a CRS with meters as a unit.
            easting_data, northing_data = [], []
            if (long_col_name in data) and (lat_col_name in data):
                easting_data, northing_data = self._get_coordinates_in_meters(
                    x_coords = data[long_col_name],
                    y_coords = data[lat_col_name],
                    crs = transect_crs if self._clicked_transects_crs in data else None
                )
            # For each data file, plot its data collected along the clicked transect.
            plot = None
            if self._data_within_crs_bounds(x_data = easting_data, y_data = northing_data, crs = transect_crs):
                # plot = await self._get_clipped_data_plots(easting_data, northing_data, long_col_name, lat_col_name, transect_crs)
                # Create a list of tasks (clip all selected data files with the selected transect) to run asynchronously.
                tasks = [asyncio.create_task(self._clip_data(file_path, easting_data, northing_data, long_col_name, lat_col_name, transect_crs)) for file_path in self.user_selected_data_files]
                # Gather the returned results of each task.
                results = asyncio.gather(*tasks)
                print(await results)
                # Overlay all clipped data files' plots.
                for clipped_data_plot in results:
                    if plot is None: plot = clipped_data_plot
                    else: plot = plot * clipped_data_plot
            print(plot)
            if plot is not None:
                self._modal_heading_pipe.event(data = [
                    "Time-Series of Data Collected Along Transect {} from {}".format(
                        transect_id, transect_file
                    ),
                    "Scroll on the axes or data area to zoom in and out of the plot."
                ])
                # Return the overlay plot containing data collected along the transect for all data files.
                return plot.opts(ylabel = self._y_axis_data_col_name)
            else:
                self._modal_heading_pipe.event(data = [
                    "No Time-Series Available",
                    "Unfortunately, no data has been collected along your selected transect (Transect {} from {}). Please select another transect or create your own transect.".format(
                        transect_id, transect_file
                    )
                ])
        elif num_transects > 1:
            # Make the selected transects' IDs readable for the error message.
            ids = sorted([str(id) for id in set(data[transect_id_col_name])])
            transect_ids = ""
            if len(ids) == 2: transect_ids = "{} and {}".format(*ids)
            else: transect_ids = ", ".join(ids[:-1]) + ", and " + str(ids[-1])
            self._modal_heading_pipe.event(data = [
                "More than One Selected Transect",
                "Transects with IDs {} from {} have been selected. Please select only one transect because the time-series of data can (currently) only be generated from one transect.".format(
                    transect_ids, transect_file
                )
            ])
        # Return an overlay plot with placeholder plots for each data file if exactly 1 transect has not been selected yet or there's no data overlapping the clicked transect.
        # ^ since DynamicMap requires callback to always return the same viewable element (in this case, Overlay)
        # ^ DynamicMap currently doesn't update when new plots are added to the initially returned plots, so placeholder/empty plots are created for each data file
        return hv.Curve(data = []) * hv.Points(data = [])

    def _update_clicked_transects_table(self, info: dict = {}) -> None:
        """
        Updates the Panel DataFrame widget with new information about the newly clicked transect(s).

        Args:
            info (dict): Dictionary containing information about the clicked transect(s)'s start and end points
        """
        # Create a new dictionary containing the new dataframe's columns.
        dataframe_cols = info.get(self._clicked_transects_table_cols, [])
        new_transects_data = {col: vals for col, vals in info.items() if col in dataframe_cols}
        # Add a new "Point Type" column to differentiate each transect's points.
        point_type_col_vals = []
        if info:
            transect_id_col_name = info.get(self._clicked_transects_id_col, "Transect ID")
            transect_id_col_vals, seen_ids = info[transect_id_col_name], set()
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

    def _update_heading_text(self, data: list[str] = ["", ""]) -> None:
        """
        Updates the heading text at the top of the popup modal by setting new values for the rendered Panel markdown objects.

        Args:
            data (list[str]): List of strings to display in the modal's heading
        """
        if data:
            title = "#### {}".format(data[0])
            details = "##### {}".format(data[1])
            self._modal_heading[0].object = title
            self._modal_heading[1].object = details

    def _update_modal_content(self, plot: any, element: any) -> None:
        """
        Triggers an event for the update_modal parameter, which in turn invokes the content() method and updates the modal content whenever an event is triggered.

        Args:
            plot (any): HoloViews object rendering the plot; this hook/method is applied after the plot is rendered
            element (any): Element rendered in the plot
        """
        self.update_modal = True
    
    @param.depends("update_collection_dir_path", watch = True)
    def _update_collection_objects(self) -> None:
        """
        Assign styles for each time-series data file in the newly selected collection directory.
        """
        self._collection_dir_path = self._data_map.selected_collection_dir_path
        # Update widgets in the "Time-Series Data" section.
        self._selected_file_groups = set()
        self._group_data_files(all_data_files = self._data_map.data_file_options)
        # Load buffer configuration file's values.
        json_file = open(os.path.join(self._collection_dir_path, self._preprocessed_data_buffer_output))
        self._buffers = json.load(json_file)
        # Update widgets in the "Transect Search Radius" section.
        self._transect_search_radius_widgets.objects = self._transect_search_radius_constant_widgets + self._get_transect_search_radius_float_inputs()
        # If any section of the time-series accordion is open, close the accordion and reopen it so that the layout resizes to fit all the new widgets.
        current_active_cards = self._time_series_controls_accordion.active
        if current_active_cards:
            self._time_series_controls_accordion.active = []
            self._time_series_controls_accordion.active = current_active_cards
    
    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @param.depends("update_modal")
    def content(self) -> pn.Column:
        """
        Returns a Panel column with components to display in the popup modal.
        """
        # Open the app's modal to display info/error message about the selected transect(s).
        self._app_template.open_modal()
        # Return the new modal content.
        display_time_series_plot = "Time-Series of Data Collected Along Transect" in self._modal_heading[0].object
        return pn.Column(
            objects = [
                *(self._modal_heading),
                pn.panel(self._time_series_plot, visible = display_time_series_plot),
                pn.Column("Selected Transect(s) Data", self._clicked_transects_table)
            ],
            sizing_mode = "stretch_width"
        )

    @property
    def clicked_transects_pipe(self) -> hv.streams.Pipe:
        """
        Returns a pipe stream that contains information about the most recently clicked transect(s).
        Other classes can pass the data into this stream by triggering an event.
        """
        return self._clicked_transects_pipe
    
    @property
    def time_series_controls(self) -> pn.Accordion:
        """
        Returns the accordion layout containing controls/widgets that let the user change settings for the time-series.
        """
        return self._time_series_controls_accordion