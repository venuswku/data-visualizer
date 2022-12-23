# Standard library imports
import os

# External dependencies imports
import param
import panel as pn
import holoviews as hv
# from holoviews import opts
import rioxarray as rxr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point, LineString
from bokeh.palettes import Set2
from .DataMap import DataMap

### PopupModal is used to display a time-series plot or any other data/message in the app's modal. ###
class PopupModal(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    user_selected_data_files = param.ListSelector(label = "Time-Series Data")
    show_time_series = param.Boolean(default = True)
    modal_content = param.List()

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, data_converter: DataMap, template: pn.template, time_series_data_col_names: list[str] = [], **params) -> None:
        """
        Creates a new instance of the PopupModal class with its instance variables.

        Args:
            data_converter (DataMap): Instance containing methods for converting data files to allow quicker loading onto a map
            template (panel.template): Data visualizer app's template
            time_series_data_col_names (list[str]): Optional list of column names for columns containing data for the time-series' y-axis
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._data_converter = data_converter
        self._app_template = template
        self._all_data_cols = time_series_data_col_names
        
        # _data_dir_path = path to the directory containing all the data for the time-series plot
        self._data_dir_path = "./data/Elwha/Time-Series Data"
        # _dist_col_name = name of the column that stores the x-axis values (distance from shore) for the time-series plot
        self._dist_col_name = "Across-Shore Distance (m)"
        # _default_y_axis_data_col_name = default name of the column that stores the y-axis values for the time-series plot (default is often used for data in ASCII grid files)
        self._default_y_axis_data_col_name = "Elevation (m)"
        # _point_type_col_name = name of the column that stores the type of transect point (either start or end)
        self._point_type_col_name = "Point Type"
        # The following list of constant variables are keys that appear in the dictionary that DataMap sends into PopupModal's _clicked_transects_pipe stream.
        # ^ When the _clicked_transects_pipe stream gets sent a new dictionary, the dictionary is passed into the _create_time_series_plot() callback as the `data` keyword argument.
        [self._clicked_transects_file, self._num_clicked_transects, self._clicked_transects_longitude_col,
        self._clicked_transects_latitude_col, self._clicked_transects_table_cols, self._clicked_transects_id_col] = data_converter.clicked_transects_info_keys

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _modal_heading_pipe = pipe stream that contains text for the popup modal's heading
        # ^ first string in list is the title of the modal
        # ^ second string in list is more detail about the modal's contents
        self._modal_heading_pipe = hv.streams.Pipe(data = [])
        # Update the _modal_heading property whenever data in the modal heading pipe changes.
        self._modal_heading_pipe.add_subscriber(self._update_heading_text)
        # _modal_heading = list containing markdown objects for the modal heading (title for first markdown, details for second markdown)
        self._modal_heading = [pn.pane.Markdown(object = ""), pn.pane.Markdown(object = "", margin = (-20, 0, 0, 5))]
        # _modal_content = vertical layout container holding all Panel objects that appear in the popup modal
        self._modal_content = pn.Column(objects = [], sizing_mode = "stretch_width")

        # _clicked_transects_pipe = pipe stream that contains info about the most recently clicked transect(s)
        self._clicked_transects_pipe = hv.streams.Pipe(data = {})
        # _y_axis_data_col_name = name of the column that stores the y-axis values for the time-series plot
        self._y_axis_data_col_name = self._default_y_axis_data_col_name
        # _time_series_plot = time-series plot for data collected along the most recently clicked transect (path)
        self._time_series_plot = hv.DynamicMap(self._create_time_series_plot, streams = [self._clicked_transects_pipe]).opts(
            title = "Time-Series",
            xlabel = self._dist_col_name,
            active_tools = ["pan", "wheel_zoom"],
            show_legend = True, toolbar = None,
            height = 500, responsive = True, padding = 0.1
        )
        self._time_series_panel_obj = pn.panel(self._time_series_plot)
        # _data_files_multichoice = custom widget that stores the user's selected data files for the time-series
        self._data_files_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.user_selected_data_files,
            placeholder = "Choose one or more data files to display in the time-series",
            solid = False, sizing_mode = "stretch_width"
        )
        # _clicked_transects_table = table containing information about the clicked transect(s)'s start and end points
        self._clicked_transects_table = pn.widgets.DataFrame(
            value = pd.DataFrame(),
            name = "Selected Transect(s) Data",
            show_index = True, auto_edit = False, text_align = "center",
            sizing_mode = "stretch_both", margin = (-20, 0, 0, 0)
        )
        
        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # Assign styles for each data file in the time-series plot.
        point_colors = list(Set2[8])
        curve_styles = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
        point_markers = ["o", "^", "s", "d", "x", ">", "*", "v", "+", "<"]
        total_colors, total_styles, total_markers = len(point_colors), len(curve_styles), len(point_markers)
        self._all_data_files, self._file_color, self._file_line, self._file_marker, i = [], {}, {}, {}, 0
        if os.path.isdir(self._data_dir_path):
            for file in os.listdir(self._data_dir_path):
                if os.path.isfile(os.path.join(self._data_dir_path, file)):
                    self._all_data_files.append(file)
                    self._file_color[file] = point_colors[i % total_colors]
                    self._file_line[file] = curve_styles[i % total_styles]
                    self._file_marker[file] = point_markers[i % total_markers]
                    i += 1
        
        # Set available options for the widget that lets the user choose what data to display in the time-series plot.
        self._data_files_multichoice.options = self._all_data_files
        self._data_files_multichoice.value = self._all_data_files

        self._loading =  pn.pane.Markdown(
            object = """
                <style>
                .loading-div {display: flex; justify-content: center; align-items: center; width: 100%; height: 100%}
                </style>
                <div class="loading-div"><h4>Loading Time-Series</h4></div>
            """,
            loading = True, visible = True,
            width = 500, height = 500, sizing_mode = "stretch_both"
        )
        # self.param.modal_content.default
        self._modal_content.objects = [
            # self._loading,
            *(self._modal_heading),
            pn.panel(self._time_series_plot, visible = True),
            pn.Row(
                pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
                self._data_files_multichoice,
                sizing_mode = "stretch_both"
            )
        ]

    # -------------------------------------------------- Private Class Methods --------------------------------------------------
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
            # Convert ASCII grid file into a new GeoTIFF (if not created yet).
            geotiff_path = self._data_dir_path + "/" + self._data_converter.geodata_dir + "/" + name + ".tif"
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
                    crs = self._data_converter.epsg
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
            ).reset_index(drop = True)
            return clipped_data_dataframe
        elif extension in [".csv", ".txt"]:
            # Convert CSV or TXT data file into a new GeoJSON (if not created yet).
            geojson_path = self._data_dir_path + "/" + self._data_converter.geodata_dir + "/" + name + ".geojson"
            self._data_converter.convert_csv_txt_data_into_geojson(self._data_dir_path + "/" + data_file_name, geojson_path)
            # Reproject the data file to match the transect's projection.
            data_geodataframe = gpd.read_file(filename = geojson_path).to_crs(crs = self._data_converter.epsg)
            # print("data_geodataframe", data_geodataframe.head())
            # Add buffer/padding to the clicked transect, which is created with the given transect's start and end point coordinates.
            # ^ Buffer allows data points within a certain distance from the clicked transect to be included in the time-series (since it's rare for data points to lie exactly on a transect).
            padded_transect = LineString(transect_points).buffer(1.0)
            # Create GeoDataFrame for the padded transect.
            clicked_transect_geodataframe = gpd.GeoDataFrame(
                data = {"geometry": [padded_transect]},
                geometry = "geometry",
                crs = self._data_converter.epsg
            )
            # print("clicked_transect_geodataframe", clicked_transect_geodataframe.head())
            # Clip data collected along the clicked transect from the given data file.
            clipped_geodataframe = data_geodataframe.clip(mask = clicked_transect_geodataframe)
            # print("clipped_geodataframe", clipped_geodataframe.head())
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
            clipped_data_dataframe = clipped_geodataframe.drop(columns = "geometry").reset_index(drop = True)
            # print("clipped_data_dataframe", clipped_data_dataframe.head())
            # Get name of the column with time-series' y-axis values.
            self._y_axis_data_col_name = self._get_data_col_name(list(clipped_data_dataframe.columns))
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
        # with pn.param.set_values(self._modal_content, loading = True):
        # with pn.param.set_values(self._loading, visible = True):
        # Get informational key-value pairs that aren't part of the time-series plot.
        transect_file = data.get(self._clicked_transects_file, None)
        num_transects = data.get(self._num_clicked_transects, 0)
        long_col_name = data.get(self._clicked_transects_longitude_col, "Longitude")
        lat_col_name = data.get(self._clicked_transects_latitude_col, "Latitude")
        transect_id_col_name = data.get(self._clicked_transects_id_col, "Transect ID")
        self._update_clicked_transects_table(info = data)
        if num_transects == 1:
            # Get the ID of the selected transect without the set's brackets.
            transect_id = str(set(data[transect_id_col_name]))[1:-1]
            # For each data file, plot its data collected along the clicked transect.
            plot = None
            for file in self._data_files_multichoice.value:
                # Clip data along the selected transect for each data file.
                clipped_dataframe = self._get_data_along_transect(
                    data_file_name = file,
                    transect_points = list(zip(data[long_col_name], data[lat_col_name], strict = True)),
                    long_col_name = long_col_name,
                    lat_col_name = lat_col_name
                )
                if clipped_dataframe is not None:
                    # Assign the time-series plot's options.
                    x_axis_col = self._dist_col_name
                    y_axis_col = self._y_axis_data_col_name
                    other_val_cols = [col for col in clipped_dataframe.columns if col not in [x_axis_col, y_axis_col]]
                    # print(self._time_series_plot.opts.info())
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
                    "Scroll on the axes or data area to zoom in and out of the plot."
                ])
                self.show_time_series = True
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
        plot = None
        for file in self._all_data_files:
            empty_curve_plot, empty_point_plot = hv.Curve(data = []), hv.Points(data = [])
            placeholder_file_plots = empty_curve_plot * empty_point_plot
            if plot is None: plot = placeholder_file_plots
            else: plot = plot * placeholder_file_plots
        self.show_time_series = False
        return plot
    
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
        num_transects = info.get(self._num_clicked_transects, 0)
        new_transects_data[self._point_type_col_name] = ["start", "end"] * num_transects
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
    
    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @property
    def clicked_transects_pipe(self) -> hv.streams.Pipe:
        """
        Returns a pipe stream that contains information about the most recently clicked transect(s).
        Other classes can pass the data into this stream by triggering an event.
        """
        return self._clicked_transects_pipe

    # @param.depends("_create_time_series_plot")#"show_time_series"
    # def _update_modal_content(self) -> None:
    #     # self._data_files_multichoice.visible = self.show_time_series
    #     # self._modal_content.objects = [
    #     #     *(self._modal_heading),
    #     #     pn.panel(self._time_series_plot, visible = self.show_time_series),
    #     #     pn.Row(
    #     #         pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
    #     #         self._data_files_multichoice
    #     #     )
    #     # ]
        
    #     # self._data_files_multichoice.visible = self.show_time_series
    #     # self._modal_content = pn.Column(
    #     #     objects = [
    #     #         *(self._modal_heading),
    #     #         pn.panel(self._time_series_plot, visible = self.show_time_series),
    #     #         pn.Row(
    #     #             pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
    #     #             self._data_files_multichoice
    #     #         )
    #     #     ],
    #     #     sizing_mode = "stretch_width"
    #     # )

    #     # self._time_series_plot.visible = self.show_time_series
    #     # self._data_files_multichoice.visible = self.show_time_series

    #     self._data_files_multichoice.visible = self.show_time_series
    #     self.modal_content = [
    #         # self._loading,
    #         *(self._modal_heading),
    #         pn.panel(self._time_series_plot, visible = self.show_time_series),
    #         pn.Row(
    #             pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
    #             self._data_files_multichoice,
    #             sizing_mode = "stretch_both"
    #         )
    #     ]

    # @property
    # @param.depends("modal_content")
    @param.depends("_create_time_series_plot")#"show_time_series"
    def content(self) -> pn.Column:
        """
        Returns a Panel column with components to display in the popup modal.
        """
        # with pn.param.set_values(self._modal_content, loading = True):
        # self._data_files_multichoice.visible = self.show_time_series
        # return pn.Column(
        #     objects = [
        #         *(self._modal_heading),
        #         pn.panel(self._time_series_plot, visible = self.show_time_series),
        #         pn.Row(
        #             pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
        #             self._data_files_multichoice
        #         )
        #     ],
        #     sizing_mode = "stretch_width"
        # )
        
        # self._data_files_multichoice.visible = self.show_time_series
        # self._modal_content.objects = [
        #     *(self._modal_heading),
        #     pn.panel(self._time_series_plot, visible = self.show_time_series),
        #     pn.Row(
        #         pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
        #         self._data_files_multichoice,
        #         sizing_mode = "stretch_both"
        #     )
        # ]
        
        # Open the app's modal to display info/error message about the selected transect(s).
        self._app_template.open_modal()
        # Update the content of all the modal components.
        self._data_files_multichoice.visible = self.show_time_series
        self._modal_content = pn.Column(
            objects = [
                # self._loading,
                *(self._modal_heading),
                pn.panel(self._time_series_plot, visible = self.show_time_series),
                pn.Row(
                    pn.Column("Selected Transect(s) Data", self._clicked_transects_table),
                    self._data_files_multichoice
                )
            ],
            sizing_mode = "stretch_width"
        )
        # Empty data in the clicked transects pipe without triggering an event
        # # (in case the user clicks on the same transect again and the pipe doesn't notice a value change).
        # self._clicked_transects_pipe.update(data = {})
        # Return the new modal content.
        return self._modal_content
        # return pn.Column(objects = self.modal_content, sizing_mode = "stretch_width")