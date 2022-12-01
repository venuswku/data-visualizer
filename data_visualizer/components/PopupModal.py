# Standard library imports
import os

# External dependencies imports
import param
import panel as pn
import holoviews as hv
from holoviews import opts
import rioxarray as rxr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from bokeh.palettes import Set2
from .DataMap import DataMap

### PopupModal is used to display a time-series plot, error message, or any other message in the app's modal. ###
class PopupModal(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, data_dir_path: str, data_converter: DataMap, **params) -> None:
        """
        Creates a new instance of the PopupModal class with its instance variables.

        Args:
            data_dir_path (str): Path to the directory containing all the data for the time-series plot
            data_converter (DataMap): Instance containing methods for converting data files to allow quicker loading onto a map
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._data_dir_path = data_dir_path
        self._data_converter = data_converter
        
        # _data_col_name = name of the column that stores the y-axis values for the time-series plot
        self._data_col_name = "Elevation"
        # _dist_col_name = name of the column that stores the x-axis values (across-shore distance) for the time-series plot
        self._dist_col_name = "Distance from Shore"
        # _point_type_col_name = name of the column that stores the type of transect point (either start or end)
        self._point_type_col_name = "Point Type"

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _clicked_transect_pipe = pipe stream that contains info about the most recently clicked transect
        self._clicked_transect_pipe = hv.streams.Pipe(data = {})
        # _time_series_plot = time-series plot for data collected along the most recently clicked transect (path)
        self._time_series_plot = hv.DynamicMap(self._create_time_series_plot, streams = [self._clicked_transect_pipe])
        # _clicked_transect_table = table containing information about the clicked transect's start and end points
        self._clicked_transect_table = hv.DynamicMap(self._create_clicked_transect_table, streams = [self._clicked_transect_pipe])
        
        # _modal_heading_pipe = pipe stream that contains text for the popup modal's heading
        # ^ first string in list is the title of the modal
        # ^ second string in list is more detail about the modal's contents
        self._modal_heading_pipe = hv.streams.Pipe(data = [])
        # Update the _modal_heading property whenever data in the modal heading pipe changes.
        self._modal_heading_pipe.add_subscriber(self._update_heading_text)
        # _modal_heading = list containing markdown objects for the modal heading (title for first markdown, details for second markdown)
        self._modal_heading = [pn.pane.Markdown(object = ""), pn.pane.Markdown(object = "")]
        
        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # Assign styles for each data file in the time-series plot.
        point_colors = list(Set2[8])
        curve_styles = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
        point_markers = ["o", "^", "s", "d", "x", ">", "*", "v", "+", "<"]
        total_colors, total_styles, total_markers = len(point_colors), len(curve_styles), len(point_markers)
        self._data_files, self._file_color, self._file_line, self._file_marker, i = [], {}, {}, {}, 0
        for file in os.listdir(self._data_dir_path):
            if os.path.isfile(os.path.join(self._data_dir_path, file)):
                self._data_files.append(file)
                self._file_color[file] = point_colors[i % total_colors]
                self._file_line[file] = curve_styles[i % total_styles]
                self._file_marker[file] = point_markers[i % total_markers]
                i += 1

    # -------------------------------------------------- Private Class Methods --------------------------------------------------
    def _get_data_along_transect(self, data_file_name: str, transect_points: list[list[float]], long_col_name: str, lat_col_name: str) -> pd.DataFrame:
        """
        Gets all data that was collected along the given transect and returns that data as a dataframe.
        Returns None if no data could be extracted with the given transect.

        Args:
            data_file_name (str): Name of the file containing data to extract for the time-series plot
            transect_points (list[list[float]]): List of coordinates for the transect's start (first item/list) and end (second item/list) points
                ^ [
                    [start point's longitude/easting, start point's latitude/northing],
                    [end point's longitude/easting, end point's latitude/northing]
                ]
            long_col_name (str): Name of the column containing the longitude/easting of each data point
            lat_col_name (str): Name of the column containing the latitude/northing of each data point
        """
        name, extension = os.path.splitext(data_file_name)
        extension = extension.lower()
        if extension == ".asc":
            geotiff_path = self._data_dir_path + "/" + self._data_converter.geodata_dir + "/" + name + ".tif"
            # Convert ASCII grid file into a new GeoTIFF (if not created yet).
            self._data_converter.convert_ascii_grid_data_into_geotiff(self._data_dir_path + "/" + data_file_name, geotiff_path)
            # Clip data collected along the clicked transect from the given data file.
            dataset = rxr.open_rasterio(geotiff_path)
            try:
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
            # Convert data into a DataFrame for easier plotting.
            clipped_dataset = clipped_dataset.squeeze().drop("spatial_ref").drop("band")
            clipped_dataset.name = self._data_col_name
            clipped_dataframe = clipped_dataset.to_dataframe().reset_index()
            no_data_val = clipped_dataset.attrs["_FillValue"]
            clipped_dataframe = clipped_dataframe[clipped_dataframe[self._data_col_name] != no_data_val]
            clipped_geodataframe = gpd.GeoDataFrame(
                data = clipped_dataframe,
                geometry = gpd.points_from_xy(
                    x = clipped_dataframe["x"],
                    y = clipped_dataframe["y"],
                    crs = self._data_converter.epsg
                )
            )
            # Calculate each point's distance from the transect's start point.
            transect_start_point = Point(transect_points[0])
            clipped_geodataframe[self._dist_col_name] = [point.distance(transect_start_point) for point in clipped_geodataframe.geometry]
            clipped_data_dataframe = clipped_geodataframe.drop(columns = "geometry").rename(
                columns = {
                    "x": long_col_name,
                    "y": lat_col_name
                }
            ).reset_index(drop = True)
            return clipped_data_dataframe
        # Return None if there's currently no implementation to extract data from the data file yet.
        print("Error extracting data along a transect from", data_file_name, ":", "Files with the", extension, "file format are not supported yet.")
        return None

    def _create_time_series_plot(self, data: dict = {}) -> hv.Overlay:
        """
        Creates a time-series plot for data collected along a clicked transect on the map.

        Args:
            data (dict): Dictionary mapping each data column (keys) to a list of values for that column (values)
        """
        # Remove informational key-value pairs that aren't part of the time-series plot.
        transect_file = data.pop("transect_file", None)
        num_transects = data.pop("num_clicked_transects", 0)
        long_col_name = data.pop("longitude_col_name", "Longitude")     # removes and returns the value with the specified key ("longitude_col_name"), else a default value ("Longitude") is returned
        lat_col_name = data.pop("latitude_col_name", "Latitude")
        # Assign the time-series plot's options.
        x_axis_col = self._dist_col_name
        y_axis_col = self._data_col_name
        other_val_cols = [long_col_name, lat_col_name]
        overlay_options = opts.Overlay(
            title = "Time-Series",
            xlabel = "Across-Shore Distance (m)",
            ylabel = "Elevation (m)",
            active_tools = ["pan", "wheel_zoom"],
            toolbar = None, show_legend = True,
            height = 500, responsive = True, padding = 0.1
        )
        if num_transects == 1:
            # Get the ID of the selected transect without the set's brackets.
            transect_id = str(set(data["Transect ID"]))[1:-1]
            # For each data file, plot its data collected along the clicked transect.
            plot = None
            for file in self._data_files:
                # Clip data along the selected transect for each data file.
                clipped_dataframe = self._get_data_along_transect(
                    data_file_name = file,
                    transect_points = list(zip(
                        data[long_col_name],
                        data[lat_col_name],
                        strict = True
                    )),
                    long_col_name = long_col_name,
                    lat_col_name = lat_col_name
                )
                if clipped_dataframe is not None:
                    # Plot clipped data.
                    clipped_data_curve_plot = hv.Curve(
                        data = clipped_dataframe,
                        kdims = x_axis_col,
                        vdims = y_axis_col,
                        label = file
                    ).opts(
                        color = self._file_color[file],
                        line_dash = self._file_line[file]
                    )
                    clipped_data_point_plot = hv.Points(
                        data = clipped_dataframe,
                        kdims = [x_axis_col, y_axis_col],
                        vdims = other_val_cols,
                        label = file
                    ).opts(
                        color = self._file_color[file],
                        marker = self._file_marker[file],
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
                    "Scroll to zoom in and out of the plot."
                ])
                # Return the overlay plot containing data collected along the transect for all data files.
                return plot.opts(overlay_options)
            else:
                self._modal_heading_pipe.event(data = [
                    "No Time-Series Available",
                    "Unfortunately, no data has been collected along your selected transect (Transect {} from {}). Please select another transect or create your own transect.".format(
                        transect_id, transect_file
                    )
                ])
        elif num_transects > 1:
            # Make the selected transects' IDs readable for the error message.
            ids = [str(id) for id in set(data["Transect ID"])]
            transect_ids = ""
            if len(ids) == 2: transect_ids = "{} and {}".format(*ids)
            else: transect_ids = ", ".join(ids[:-1]) + "and " + str(ids[-1])
            self._modal_heading_pipe.event(data = [
                "More than One Selected Transect",
                "Transects with IDs {} from {} have been selected. Please select only one transect because time-series of data can (currently) only be generated from one transect.".format(
                    transect_ids, transect_file
                )
            ])
        # Return an overlay plot with placeholder plots for each data file if exactly 1 transect has not been selected yet or there's no data overlapping the clicked transect.
        # ^ since DynamicMap requires callback to always return the same viewable element (in this case, Overlay)
        # ^ DynamicMap currently doesn't update plots properly when new plots are added to the initially returned plots, so placeholder/empty plots are created for each data file
        plot = None
        for file in self._data_files:
            empty_curve_plot = hv.Curve(data = [], kdims = x_axis_col, vdims = y_axis_col, label = file)
            empty_point_plot = hv.Points(data = [], kdims = [x_axis_col, y_axis_col], vdims = other_val_cols, label = file).opts(tools = ["hover"])
            placeholder_file_plots = empty_curve_plot * empty_point_plot
            if plot is None: plot = placeholder_file_plots
            else: plot = plot * placeholder_file_plots
        return plot.opts(overlay_options)

    def _create_clicked_transect_table(self, data: dict = {}) -> hv.Table:
        """
        Creates a table containing information about the clicked transect's start and end points.

        Args:
            data (dict): Dictionary mapping each data column (keys) to a list of values for that column (values)
        """
        # Add a new "Point Type" column to differentiate the transect points.
        data[self._point_type_col_name] = ["start", "end"]
        return hv.Table(
            data = data,
            kdims = [self._point_type_col_name],
            vdims = data.pop("table_cols", [])
        ).opts(
            title = "Selected Transect's Data",
            editable = False, fit_columns = True
        )

    def _update_heading_text(self, data: list[str] = ["", ""]) -> None:
        """
        Updates the heading text at the top of the popup modal by setting new values for the rendered Panel markdown objects.

        Args:
            data (list[str]): List of strings to display in the modal's heading
        """
        title = "#### {}".format(data[0])
        details = "##### {}".format(data[1])
        self._modal_heading[0].object = title
        self._modal_heading[1].object = details
    
    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @property
    def clicked_transect_pipe(self) -> hv.streams.Pipe:
        """
        Returns a pipe stream that contains information about the most recently clicked transect.
        Other objects can pass the data into this stream by triggering an event.
        """
        return self._clicked_transect_pipe

    @property
    def time_series_plot(self) -> hv.DynamicMap:
        """
        Returns a time-series plot for data collected near the selected transect on the map.
        """
        return self._time_series_plot
    
    # @param.depends("_create_clicked_transect_table")
    @property
    def clicked_transect_data(self) -> hv.DynamicMap:
        """
        Returns a table for the selected transect on the map.
        """
        # table_dynamic_map_vals = self._clicked_transect_table.values()
        # if len(table_dynamic_map_vals):
        #     holoviews_table = table_dynamic_map_vals[0]
        #     print(holoviews_table)
        #     # Create a pandas dataframe containing the clicked transect's data.
        #     col_names = holoviews_table.dimensions(
        #         selection = "all",
        #         label = "name"
        #     )
        #     print(col_names)
        #     clicked_transect_dataframe = holoviews_table.dframe(
        #         dimensions = col_names,
        #         # multi_index = True
        #     )
        #     # dict = holoviews_table.columns(dimensions = ["Point Type", "Easting (meters)", "Northing (meters)", "Transect ID"])
        #     # print(dict)
        #     # clicked_transect_dataframe = pd.DataFrame(dict)
        #     print("table dataframe", clicked_transect_dataframe)

        # # Return a customized Panel DataFrame widget.
        # return pn.widgets.DataFrame(
        #     self._clicked_transect_dataframe,
        #     name = "Selected Transect's Data",
        #     show_index = True, auto_edit = False, text_align = "center"
        # )
        
        # print(self._clicked_transect_table.last)
        # if self._clicked_transect_table.last is None:
        #     return self._clicked_transect_table
        # else:
        #     print("dframe", self._clicked_transect_table.last.dframe())
        #     widget = pn.widgets.DataFrame(
        #         value = self._clicked_transect_table.last.dframe(),
        #         name = "Selected Transect's Data",
        #         show_index = True, auto_edit = False, text_align = "center"
        #     )
        #     print(widget)
        return self._clicked_transect_table

    @property
    def content(self) -> pn.Column:
        """
        Returns a Panel column with components to display in the popup modal.
        """
        return pn.Column(
            objects = [
                *(self._modal_heading),
                self._time_series_plot,
                self._clicked_transect_table
            ],
            sizing_mode = "stretch_width"
        )