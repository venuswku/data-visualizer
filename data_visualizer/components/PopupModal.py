# Standard library imports
import os
import json
from collections import Counter

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import holoviews as hv
import rioxarray as rxr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString
import cartopy.crs as ccrs
import pyproj
from shapely.ops import transform
from bokeh.palettes import Set2
from bokeh.models.formatters import PrintfTickFormatter
from .DataMap import DataMap

### PopupModal is used to display a time-series plot or any other data/message in the app's modal. ###
class PopupModal(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    update_collection_dir_path = param.Event(label = "Action that Triggers the Updating of the Collection Directory and Its Related Objects")
    user_selected_data_files = param.ListSelector(label = "Time-Series Data Files")
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
        # _file_color = directory mapping each data file option (key) in _all_data_files to a color (value) in the time-series plot
        self._file_color = {}
        # _file_line = directory mapping each data file option (key) in _all_data_files to a line style (value) in the time-series plot
        self._file_line = {}
        # _file_marker = directory mapping each data file option (key) in _all_data_files to a marker (value) in the time-series plot
        self._file_marker = {}

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
        # # _clicked_transects_plot = plot containing the user's clicked transect(s) with its buffer if extracting point data for the time-series
        # self._clicked_transects_plot = hv.DynamicMap(self._create_clicked_transects_plot, streams = [self._clicked_transects_pipe])
        # _buffers = dictionary mapping each data file's path (key) to the selected transect's buffer/search radius (value) when extracting data along the transect
        self._buffers = {}
        
        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        self._point_colors = list(Set2[8])
        self._curve_styles = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
        self._point_markers = ["o", "^", "s", "d", "*", "+"]
        self._total_colors, self._total_styles, self._total_markers = len(self._point_colors), len(self._curve_styles), len(self._point_markers)
        
        # _data_files_checkbox_group = custom widget that stores the user's selected data files for the time-series
        self._data_files_checkbox_group = pn.widgets.CheckBoxGroup.from_param(parameter = self.param.user_selected_data_files)
        # _update_buffer_config_file_button = button for updating the collection's buffer config file with values from the buffer float widgets for each data file
        self._update_buffer_config_file_button = pn.widgets.Button.from_param(
            parameter = self.param.update_buffer_config,
            name = "Save Current Values to {}".format(self._preprocessed_data_buffer_output),
            button_type = "primary", disabled = True
        )
        # _transect_search_radius_constant_widgets = list of widgets that always appear at the top of the "Transect Search Radius" accordion section
        self._transect_search_radius_constant_widgets = [
            pn.widgets.StaticText(value = "Adjust the search radius for extracting time-series data around a selected transect."),
            self._update_buffer_config_file_button
        ]
        # _transect_search_radius_widgets = column layout containing widgets for the "Transect Search Radius" accordion section
        self._transect_search_radius_widgets = pn.Column(objects = [])
        # _time_series_controls_accordion = accordion layout widget allowing the user to change settings for the time-series
        self._time_series_controls_accordion = pn.Accordion(
            objects = [
                ("Time-Series Data", self._data_files_checkbox_group),
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
    def _save_changed_buffer_val(self, event: param.parameterized.Event) -> None:
        """
        Updates the buffers dictionary whenever any of the float input widgets (for each data file) change value.
        """
        self._update_buffer_config_file_button.disabled = False
        # Get the path to the data file that corresponds to the float input widget with the changed value.
        subdir, data_file = (event.obj.name).split(": ")
        data_file_path = os.path.join(self._collection_dir_path, subdir, data_file)
        # Save the new buffer value.
        self._buffers[data_file_path] = event.new

    def _get_transect_search_radius_float_inputs(self) -> pn.Column:
        """
        Returns a list of float input widgets for each data file's buffer value.
        """
        widgets = []
        for data_file_path, buffer in self._buffers.items():
            # Create float input widget.
            data_base_path, data_file = os.path.split(data_file_path)
            subdir_name = os.path.basename(data_base_path)
            file_buffer_float_input = pn.widgets.FloatInput(name = "{}: {}".format(subdir_name, data_file), value = buffer, start = 0)
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
            crs = self._data_map.collection_crs
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
    
    def _get_data_along_transect(self, data_file_option: str, transect_points: list[list[float]], long_col_name: str, lat_col_name: str, transect_crs: ccrs) -> pd.DataFrame:
        """
        Gets all data that was collected along the given transect and returns that data as a dataframe.
        Returns None if no data could be extracted with the given transect.

        Args:
            data_file_option (str): Name of the option for the file containing data to extract for the time-series plot
            transect_points (list[list[float]]): List of coordinates for the transect's start (first item/list) and end (second item/list) points
                ^ [
                    [start point's longitude/easting, start point's latitude/northing],
                    [end point's longitude/easting, end point's latitude/northing]
                ]
            long_col_name (str): Name of the column containing the longitude/easting of each data point
            lat_col_name (str): Name of the column containing the latitude/northing of each data point
            transect_crs (cartopy.crs): Coordinate reference system of the given transect
        """
        subdir, data_file = data_file_option.split(": ")
        data_file_subpath = os.path.join(subdir, data_file)
        _, extension = os.path.splitext(data_file)
        extension = extension.lower()
        if extension in [".tif", ".tiff"]:
            data_file_path = os.path.join(self._collection_dir_path, data_file_subpath)
            dataset = rxr.open_rasterio(data_file_path)
            # Clip data collected along the clicked transect from the given data file.
            try:
                transect_buffer = self._buffers.get(data_file_path, 0)
                if transect_buffer > 0:
                    padded_transect_points = LineString(transect_points).buffer(transect_buffer)
                    clipped_dataset = dataset.rio.clip(
                        geometries = [{
                            "type": "Polygon",
                            "coordinates": padded_transect_points
                        }],
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
            data_file_path = os.path.join(self._collection_dir_path, data_file_subpath)
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

    def _create_time_series_plot(self, data: dict = {}) -> hv.Overlay:
        """
        Creates a time-series plot for data collected along a clicked transect on the map.

        Args:
            data (dict): Dictionary containing information about the clicked transect(s)
                and maps each data column (keys) to a list of values for that column (values)
        """
        # Get informational key-value pairs that aren't part of the time-series plot.
        transect_file = data.get(self._clicked_transects_file, None)
        num_transects = data.get(self._num_clicked_transects, 0)
        transect_crs = data.get(self._clicked_transects_crs, self._data_map.collection_crs)
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
                for file_option in self._data_files_checkbox_group.value:
                    # Clip data along the selected transect for each data file.
                    clipped_dataframe = self._get_data_along_transect(
                        data_file_option = file_option,
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
                            color = self._file_color[file_option],
                            line_dash = self._file_line[file_option]
                        )
                        clipped_data_point_plot = hv.Points(
                            data = clipped_dataframe,
                            kdims = [x_axis_col, y_axis_col],
                            vdims = other_val_cols,
                            label = file_option
                        ).opts(
                            color = self._file_color[file_option],
                            marker = self._file_marker[file_option],
                            tools = ["hover"],
                            size = 10
                        )
                        # Add the data file's plot to the overlay plot.
                        clipped_data_plot = clipped_data_curve_plot * clipped_data_point_plot
                        if plot is None: plot = clipped_data_plot
                        else: plot = plot * clipped_data_plot
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
    
    # def _create_clicked_transects_plot(self, data: dict = {}) -> hv.Overlay:
    #     """
    #     Creates a plot containing the user's clicked transect(s) with its buffer if extracting point data for the time-series.

    #     Args:
    #         data (dict): Dictionary containing information about the clicked transect(s)
    #             and maps each data column (keys) to a list of values for that column (values)
    #     """
    #     all_transect_plots = None
    #     num_transects = data.get(self._num_clicked_transects, 0)
    #     transect_crs = data.get(self._clicked_transects_crs, self._data_map.map_default_crs)
    #     long_col_name = data.get(self._clicked_transects_longitude_col, "Longitude")
    #     lat_col_name = data.get(self._clicked_transects_latitude_col, "Latitude")
    #     transect_id_col_name = data.get(self._clicked_transects_id_col, "Transect ID")
    #     geometry_col_name, buffer_col_name = "geometry", "Search Radius"
    #     start_pt_col_name, end_pt_col_name = "Start Point", "End Point"
    #     # Plot each clicked transect.
    #     start_pt_index, next_start_pt_index = 0, 1
    #     for _ in range(num_transects):
    #         transect_id_col_vals = data[transect_id_col_name]
    #         while (next_start_pt_index < len(transect_id_col_vals)) and (transect_id_col_vals[start_pt_index] == transect_id_col_vals[next_start_pt_index]):
    #             next_start_pt_index += 1
    #         transect_pts = list(zip(
    #             data[long_col_name][start_pt_index: next_start_pt_index],
    #             data[lat_col_name][start_pt_index: next_start_pt_index],
    #             strict = True
    #         ))
    #         transect_id = transect_id_col_vals[start_pt_index]
    #         # Create a GeoDataFrame to plot the clicked transect as a path.
    #         clicked_transect = LineString(transect_pts)
    #         clicked_transect_geodataframe = gpd.GeoDataFrame(
    #             data = {
    #                 geometry_col_name: [clicked_transect],
    #                 transect_id_col_name: [transect_id],
    #                 start_pt_col_name: ["({}, {})".format(*(transect_pts[0]))],
    #                 end_pt_col_name: ["({}, {})".format(*(transect_pts[-1]))]
    #             },
    #             geometry = geometry_col_name,
    #             crs = transect_crs
    #         )
    #         transect_plot = gv.Path(
    #             data = clicked_transect_geodataframe,
    #             vdims = [transect_id_col_name, start_pt_col_name, end_pt_col_name],
    #             crs = transect_crs
    #         ).opts(
    #             color = self._data_map.app_main_color,
    #             tools = ["hover"]
    #         )
    #         plot = transect_plot
    #         # # Plot the buffer/padding/search radius around the clicked transect only if there's just one selected transect.
    #         # if num_transects == 1:
    #         #     clicked_transect_buffer = 3
    #         #     if self._clicked_transects_crs in data:
    #         #         buffer = clicked_transect.buffer(clicked_transect_buffer)
    #         #     else:
    #         #         # Transform any user-drawn transect's coordinates into a CRS with meters as a unit (for plotting an accurate buffer).
    #         #         easting_data, northing_data = [], []
    #         #         if (long_col_name in data) and (lat_col_name in data):
    #         #             easting_data, northing_data = self._get_coordinates_in_meters(
    #         #                 x_coords = data[long_col_name],
    #         #                 y_coords = data[lat_col_name],
    #         #                 crs = None
    #         #             )
    #         #         transformed_transect_pts = list(zip(easting_data, northing_data, strict = True))
    #         #         transformed_buffer = LineString(transformed_transect_pts).buffer(clicked_transect_buffer)
    #         #         # Transform the buffer back into the data map's default CRS in case it lies outside of the data CRS's bounds.
    #         #         projection = pyproj.Transformer.from_crs(
    #         #             crs_from = self._data_map.collection_crs,
    #         #             crs_to = self._data_map.map_default_crs,
    #         #             always_xy = True
    #         #         ).transform
    #         #         buffer = transform(projection, transformed_buffer)
    #         #     # Create the buffer plot.
    #         #     buffer_plot = gv.Polygons(
    #         #         data = gpd.GeoDataFrame(
    #         #             data = {
    #         #                 geometry_col_name: [buffer],
    #         #                 buffer_col_name: [clicked_transect_buffer]
    #         #             },
    #         #             geometry = geometry_col_name,
    #         #             crs = transect_crs
    #         #         ),
    #         #         vdims = [buffer_col_name],
    #         #         crs = transect_crs
    #         #     ).opts(
    #         #         color = self._data_map.app_main_color,
    #         #         line_color = self._data_map.app_main_color,
    #         #         alpha = 0.1, tools = ["hover"]
    #         #     )
    #         #     plot = plot * buffer_plot
    #         # else:
    #         #     plot = plot * gv.Polygons(data = [], crs = transect_crs)
    #         # Save the plots for each transect.
    #         if all_transect_plots is None: all_transect_plots = plot
    #         else: all_transect_plots = all_transect_plots * plot
    #         # Assign the next transect's start point index.
    #         start_pt_index = next_start_pt_index
    #     # Return an overlay plot containing all clicked transects.
    #     if all_transect_plots is None: return gv.Path(data = [])# * gv.Polygons(data = [])
    #     else: return all_transect_plots

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
        # Load buffer configuration file's values.
        json_file = open(os.path.join(self._collection_dir_path, self._preprocessed_data_buffer_output))
        self._buffers = json.load(json_file)
        # Update widgets in "Transect Search Radius" section.
        self._transect_search_radius_widgets.objects = self._transect_search_radius_constant_widgets + self._get_transect_search_radius_float_inputs()
        # Get all data files' widget option names (i.e. "{subdirectory in _collection_dir_path}: {data file name}") from collection directory.
        i, data_file_options_dict = 0, {}
        collection_subdirs = [file for file in os.listdir(self._collection_dir_path) if os.path.isdir(os.path.join(self._collection_dir_path, file)) and (file != self._data_map.transects_dir_name)]
        for subdir in collection_subdirs:
            subdir_path = os.path.join(self._collection_dir_path, subdir)
            # data_file_options_dict[subdir] = subdir_path
            for file in [file for file in os.listdir(subdir_path) if os.path.isfile(os.path.join(subdir_path, file))]:
                file_option_name = ": ".join([subdir, file])
                data_file_options_dict[file] = file_option_name
                # Set styles for each data file's plot in the time-series.
                self._file_color[file_option_name] = self._point_colors[i % self._total_colors]
                self._file_line[file_option_name] = self._curve_styles[i % self._total_styles]
                self._file_marker[file_option_name] = self._point_markers[i % self._total_markers]
                i += 1
        # Set available options for the widget that lets the user choose what data to display in the time-series plot.
        self._data_files_checkbox_group.options = data_file_options_dict
        self._data_files_checkbox_group.value = []
        # If the time-series data accordion is open, close the accordion and reopen it so that the layout resizes to fit all the new data files.
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
                pn.Row(
                    pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
                    # pn.panel(
                    #     (self._data_map.selected_basemap * self._clicked_transects_plot).opts(
                    #         toolbar = None, xaxis = None, yaxis = None, title = "", responsive = True
                    #     ),
                    #     sizing_mode = "stretch_both"
                    # )
                )
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