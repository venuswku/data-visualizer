# Standard library imports
import os
import json
from datetime import datetime

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
import holoviews as hv
from holoviews.operation.datashader import regrid, datashade
import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import LineString
from bokeh.palettes import Bokeh, Set2
from io import BytesIO

### DataMap is used for displaying the inputted data files onto a map. ###
class DataMap(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    basemap = param.Selector(label = "Basemap")
    collection = param.Selector(label = "Collection")
    transects = param.ListSelector(label = "Transects")
    clicked_transects_info = param.Dict(default = {}, label = "Information About the Recently Clicked Transect(s)")
    view_user_transect_time_series = param.Event(label = "Indicator for Displaying the Time-Series for Data Along the User-Drawn Transect")
    data_file_path = param.Path(default = None, label = "Path to the Most Recently Selected Data File from PopupModal's _data_files_checkbox_group Widget", allow_None = True)

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, **params) -> None:
        """
        Creates a new instance of the DataMap class with its instance variables.
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        # _root_data_dir_path = path to the root directory that contains all available datasets/collections for the app
        self._root_data_dir_path = os.path.abspath("./data")
        # _default_crs = default coordinate reference system for the user-drawn transect and other plots
        self._default_crs = ccrs.PlateCarree()
        # _app_main_color = theme color used for all the Panel widgets in this app
        self._app_main_color = "#2196f3"
        
        # _all_basemaps = dictionary mapping each basemap name (keys) to a basemap WMTS (web mapping tile source) layer (values)
        self._all_basemaps = {
            "Default": gts.OSM,
            "Satellite": gts.EsriImagery,
            "Topographic": gts.OpenTopoMap,
            "Black & White": gts.StamenToner,
            "Dark": gts.CartoDark
        }
        # _all_collections = list of provided collections (subfolders in the root data directory)
        self._all_collections = {}
        for file in os.listdir(self._root_data_dir_path):
            file_path = os.path.join(self._root_data_dir_path, file)
            if os.path.isdir(file_path):
                collection_info_json_path = os.path.join(file_path, "collection_info.json")
                json_file = open(collection_info_json_path)
                collection_info = json.load(json_file)
                if file in collection_info:
                    collection_title = collection_info[file]
                    self._all_collections[collection_title] = file
                else:
                    self._all_collections[file] = file

        # _create_own_transect_option = Name of the option for the user to create their own transect
        self._create_own_transect_option = "Draw My Own Transect"
        # _clicked_transects_file_key = clicked_transects_info parameter's dictionary key that corresponds to the file that contains the clicked transect(s)
        self._clicked_transects_file_key = "transect_file"
        # _num_clicked_transects_key = clicked_transects_info parameter's dictionary key that corresponds to the number of clicked transect(s)
        self._num_clicked_transects_key = "num_clicked_transects"
        # _clicked_transects_crs_key = clicked_transects_info parameter's dictionary key that corresponds to the clicked transect(s)'s coordinate reference system
        self._clicked_transects_crs_key = "crs"
        # _clicked_transects_longitude_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the longitude/easting of the clicked transects' start and end points
        self._clicked_transects_longitude_key = "longitude_col_name"
        # _clicked_transects_latitude_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the latitude/northing of the clicked transects' start and end points
        self._clicked_transects_latitude_key = "latitude_col_name"
        # _clicked_transects_data_cols_key = clicked_transects_info parameter's dictionary key that corresponds to the names of the columns to display in the popup modal's clicked transect(s) data table
        self._clicked_transects_data_cols_key = "table_cols"
        # _clicked_transects_id_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the ID of the clicked transects' start and end points
        self._clicked_transects_id_key = "id_col_name"
        # _transects_folder_name = Name of the folder containing files with transect data
        # ^ should be same as `transects_subdir_name` in utils/preprocess_data.py
        self._transects_folder_name = "Transects"
        # _transects_id_col_name = Name of the column containing the ID of each clicked transect
        # ^ should be same as `transect_geojson_id_property` in utils/preprocess_data.py because it was used to assign the ID property for each transect in the outputted GeoJSON
        self._transects_id_col_name = "Transect ID"
        # _transect_start_point_prop_name = Name of the GeoJSON property containing the start point of a transect
        # ^ should be same as `transect_geojson_start_point_property` in utils/preprocess_data.py because it was used to assign the start point property for each transect in the outputted GeoJSON
        self._transect_start_point_prop_name = "Start Point"
        # _transect_end_point_prop_name = Name of the GeoJSON property containing the end point of a transect
        # ^ should be same as `transect_geojson_end_point_property` in utils/preprocess_data.py because it was used to assign the end point property for each transect in the outputted GeoJSON
        self._transect_end_point_prop_name = "End Point"

        # Colors and markers for data in the map and time-series plots in the popup modal:
        self._palette1_colors = Bokeh[8]
        self._palette2_colors = list(Set2[8])
        self._markers = ["o", "^", "s", "d", "*", "+"]
        self._curve_styles = ["solid", "dashed", "dotted", "dotdash", "dashdot"]
        self._total_palette1_colors = len(self._palette1_colors)
        self._total_palette2_colors = len(self._palette2_colors)
        self._total_markers = len(self._markers)
        self._total_line_styles = len(self._curve_styles)
        
        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _data_map_plot = overlay plot containing the selected basemap and all the data (categories, transects, etc.) plots
        self._data_map_plot = pn.pane.HoloViews(object = None, sizing_mode = "scale_both")
        # _created_plots = dictionary mapping each file's path (keys) to its created plot (values)
        self._created_plots = {}
        
        # _selected_basemap_plot = WMTS (web mapping tile source) layer containing the user's selected basemap
        self._selected_basemap_plot = list(self._all_basemaps.values())[0]
        
        # _collection_dir_path = path to the selected collection's directory
        self._collection_dir_path = None
        # _collection_crs = coordinate reference system of the chosen collection
        self._collection_crs = None
        # _selected_collection_info = information about the selected collection, which is loaded from its collection_info.json file
        # ^ name of the JSON file should be same as `outputted_collection_json_name` in utils/preprocess_data.py
        self._selected_collection_info = {}
        # _data_file_options_dict = dictionary mapping the option name (key) of each data file in the selected collection to the data file's path (value)
        self._data_file_options_dict = {}

        # _data_file_color = directory mapping each data file path (key) to a color (value) for the map and time-series plots
        self._data_file_color = {}
        # _data_file_marker = directory mapping each data file path (key) to a marker (value) for the map and time-series plots
        self._data_file_marker = {}
        # _data_file_line_style = directory mapping each data file path (key) to a line style (value) for the time-series plot
        self._data_file_line_style = {}
        # _selected_data_plot = point or image plot for the most recently selected data file
        # ^ None if the no data file was selected by the user
        self._selected_data_plot = None
        
        # _all_transect_files = list of files containing transects to display on the map
        self._all_transect_files = []
        # _transect_colors = dictionary mapping each transect file (key) to a color (value), which will be used for the color of its path plots
        self._transect_colors = {}
        # _selected_transects_plot = overlay of path plots if the user selected one or more transect files to display on the map
        # ^ None if the user didn't provide any transect files or the user didn't select to display a transect file
        self._selected_transects_plot = None
        # _tapped_data_streams = dictionary mapping each transect file's path (key) to a selection stream (value), which saves the file's most recently clicked data element (path) on the map
        self._tapped_data_streams = {}

        # _user_transect_plot = path plot used when the user wants to create their own transect to display on the map
        self._user_transect_plot = gv.Path(data = [], crs = self._default_crs, label = self._create_own_transect_option)#.opts(projection = self._default_crs)
        # _edit_user_transect_stream = stream that allows user to add and move the start and end points of their own transect
        self._edit_user_transect_stream = hv.streams.PolyDraw(
            source = self._user_transect_plot,
            num_objects = 1,
            styles = {"line_color": ["black"], "line_width": [5]},
            vertex_style = {"fill_color": "black"},
            show_vertices = True, drag = True
        )
        self._edit_user_transect_stream.add_subscriber(self._set_ability_to_save_user_transect)

        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # Set basemap widget's options.
        self.param.basemap.objects = self._all_basemaps.keys()

        # Set collection widget's options using the _all_collections dictionary.
        self._collection_select = pn.widgets.Select.from_param(
            parameter = self.param.collection,
            options = self._all_collections
        )

        # Set transect widget's options.
        self._transects_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.transects,
            options = self._all_transect_files + [self._create_own_transect_option],
            placeholder = "Choose one or more transect files to display",
            solid = False
        )
        self._update_collection_objects()
        # Show an error popup if there are any errors that occurred while creating plots for the data map or popup modal.
        self._error_messages = []
        self._error_popup_text = pn.widgets.TextInput(value = "", visible = False)
        self._error_popup_text.jscallback(
            value = """
            if (cb_obj.value) {
                window.alert(cb_obj.value);
                cb_obj.value = ""
            }
            """
        )

        # Create a button for viewing the time-series for data along the user-drawn transect.
        self._view_user_transect_time_series_button = pn.widgets.Button.from_param(
            parameter = self.param.view_user_transect_time_series,
            name = "View Time-Series for Drawn Transect",
            visible = False, disabled = True, loading = False, button_type = "primary"
        )
        # Create a button for downloading the user-drawn transect as a GeoJSON.
        self._user_drawn_transect_download_button = pn.widgets.FileDownload(
            filename = "drawn_transect.geojson",
            callback = self._get_user_drawn_transect_geojson,
            label = "Save Drawn Transect",
            visible = False, disabled = True, loading = False, button_type = "default"
        )
        # Create a tips section and a dropdown for instructions on how to use the PolyDraw tool.
        self._drawing_user_transect_instructions = pn.Column(
            pn.pane.Alert(object = """
                <b>Tips</b>
                <ul>
                    <li>Repeat the <b>Add</b> steps below to simultaneously delete an existing transect and add a new one.</li>
                    <li>Scroll down towards the bottom of the map to view all the available tools.
                        <ul>
                            <li>The Pan tool <img src="https://raw.githubusercontent.com/venuswku/data-visualizer/main/assets/PanTool.png" alt="Pan tool" width="16"/> allows you to move around the map by holding down your mouse/trackpad and dragging.</li>
                            <li>
                                The Pan tool is automatically disabled whenever the Polygon Draw tool <img src="https://raw.githubusercontent.com/venuswku/data-visualizer/main/assets/PolygonDrawTool.png" alt="Polygon Draw tool" width="16"/> is enabled (allowing you to draw your own transect), and vice versa. 
                                Click on either tool to toggle between them.
                            </li>
                        </ul>
                    </li>
                </ul>
            """, alert_type = "primary"),
            pn.Accordion((
                "Draw Your Own Transects",
                pn.pane.Markdown(object = """
                    <b>Add</b>
                    <ul>
                        <li>Double click to add the start point.</li>
                        <li>(Optional) Single click to add each subsequent point.</li>
                        <li>If you want to restart adding transect points, press the ESC key.</li>
                        <li>Double click to add the end point.</li>
                    </ul>
                    <b>Move</b>
                    <ul>
                        <li>Click to select an existing transect.</li>
                        <li>Then drag the transect to move it.</li>
                        <li>Transect points will be moved once you let go of the mouse/trackpad.</li>
                    </ul>
                    <b>Delete</b>
                    <ul>
                        <li>Click to select an existing transect.</li>
                        <li>Then press the BACKSPACE (Windows) or DELETE (Mac) key while the cursor is within the map area.</li>
                    </ul>
                """, sizing_mode = "stretch_width")
            )), visible = False, margin = 0
        )        

    # -------------------------------------------------- Private Class Methods --------------------------------------------------    
    def _plot_geojson_points(self, data_file_option: str) -> gv.Points:
        """
        Creates a point plot from a GeoJSON file containing Points.

        Args:
            data_file_option (str): Option name of the most recently selected data file from PopupModal's _data_files_checkbox_group widget
        """
        # Read the GeoJSON as a GeoDataFrame.
        geodataframe = gpd.read_file(self.data_file_path)
        latitude_col, longitude_col, non_lat_long_cols = None, None, []
        for col in geodataframe.columns:
            col_name = col.lower()
            if "lat" in col_name: latitude_col = col
            elif "lon" in col_name: longitude_col = col
            elif col_name != "geometry": non_lat_long_cols.append(col)
        # Create a point plot with the GeoDataFrame.datashade()
        point_plot = gv.Points(
            data = geodataframe,
            kdims = [longitude_col, latitude_col],
            vdims = non_lat_long_cols,
            label = data_file_option
        ).opts(
            color = self._data_file_color[self.data_file_path],
            marker = self._data_file_marker[self.data_file_path],
            hover_color = self._app_main_color,
            tools = ["hover"],# responsive = True,
            size = 10, muted_alpha = 0.01
        )
        return point_plot
    
    def _plot_geojson_linestrings(self, geojson_file_path: str, filename: str) -> gv.Path:
        """
        Creates a path plot from a GeoJSON file containing LineStrings. Returns None if the a transect is invalid (e.g. less than 2 points).

        Args:
            geojson_file_path (str): Path to the GeoJSON file containing LineStrings
            filename (str): Name of the transect file that corresponds to the returned path plot
        """
        # Convert the CRS from the GeoJSON file into GeoView's default CRS (Plate Carree), if necessary.
        geodataframe = gpd.read_file(geojson_file_path)
        geojson_crs = geodataframe.crs
        if geojson_crs is not None:
            geojson_epsg_code = ccrs.CRS(geojson_crs).to_epsg()
            # Only use projected coordinate systems, not geodetic coordinate systems like EPSG:4326/WGS-84 (https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.epsg.html).
            if geojson_epsg_code not in [4326]: geodataframe = geodataframe.to_crs(crs = self._default_crs)
        # Check if LineStrings/transects contain 2 or more points.
        if any(map(lambda transect: len(transect.coords) < 2, geodataframe.geometry)):
            # Set the message of the error popup.
            error = "Error displaying {} as a transect plot: Found an invalid transect containing less than two points. Transects are expected to have two or more points.".format(filename)
            self._error_messages.append(error)
            # Print an error message if the given file contains an invalid transect.
            print(error)
            return None
        elif any(map(lambda transect: len(transect.coords) > 2, geodataframe.geometry)):
            # If any of the file's transects contains more than 2 points, add a color column for the transect plot's lines.
            transect_color = self._transect_colors[filename]
            color_col_name = "color"
            geodataframe.insert(
                loc = len(geodataframe.columns),
                column = color_col_name,
                value = [transect_color] * len(geodataframe.index)
            )
            # Create a contour plot in order to avoid the transect from being split into sub-geometries when plotted as a path
            # (sub-geometries/segments give _get_clicked_transect_info() the wrong index when only a segment is clicked).
            return gv.Contours(
                data = geodataframe,
                # vdims = [col for col in list(geodataframe.columns) if col != color_col_name],
                crs = self._default_crs,
                # label = "{}: {}".format(self._transects_folder_name, filename)
            ).opts(
                color = color_col_name,
                hover_color = self._app_main_color,
                selection_color = self._app_main_color,
                nonselection_color = self._transect_colors[filename],
                nonselection_alpha = 1, selection_alpha = 1,
                tools = ["hover", "tap"]
            )
        else:
            # Create a path plot with the correct CRS.
            return gv.Path(
                data = geodataframe,
                crs = self._default_crs,
                label = "{}: {}".format(self._transects_folder_name, filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
            ).opts(
                color = self._transect_colors[filename],
                hover_color = self._app_main_color,
                selection_color = self._app_main_color,
                nonselection_color = self._transect_colors[filename],
                nonselection_alpha = 1, selection_alpha = 1,
                tools = ["hover", "tap"]
            )
    
    def _create_data_plot(self) -> None:
        """
        Creates a point/image plot containing the given file's data.
        """
        # Read the file and create a plot from it.
        subdir_path, filename = os.path.split(self.data_file_path)
        _, extension = os.path.splitext(filename)
        extension = extension.lower()
        plot = None
        if extension == ".geojson":
            # Create a point plot with the GeoJSON.
            subdir_name = os.path.basename(subdir_path)
            plot = self._plot_geojson_points(data_file_option = "{}: {}".format(subdir_name, filename))
        elif extension in [".tif", ".tiff"]:
            # Create an image plot with the GeoTIFF.regrid()
            plot = gv.load_tiff(
                self.data_file_path,
                vdims = "Elevation (meters)",
                nan_nodata = True
            ).opts(
                cmap = "Turbo",
                tools = ["hover"],
                alpha = 0.5,
                # responsive = True
            )
        if plot is None:
            print("Error displaying", filename, "as a point/image plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            # Save the created plot.
            self._created_plots[self.data_file_path] = plot
    
    def _create_path_plot(self, filename: str) -> None:
        """
        Creates a path/contour plot containing the given file's transects.

        Args:
            filename (str): Name of the file containing transects
        """
        # Read the given file.
        file_path = os.path.join(self._collection_dir_path, self._transects_folder_name, filename)
        _, extension = os.path.splitext(filename)
        extension = extension.lower()
        # Create a path/contour plot from the given file.
        plot = None
        if extension == ".geojson":
            plot = self._plot_geojson_linestrings(file_path, filename)
        else:
            print("Error displaying", filename, "as a transect plot:", "Input files with the", extension, "file format are not supported yet.")
        # Save the transect plot, if created.
        if plot is not None: self._created_plots[file_path] = plot

    def _get_clicked_transect_info(self, **params: dict) -> None:
        """
        Gets information about the most recently clicked transect on the map, which is used to update the popup modal's contents (time-series plot and transect data table).

        Args:
            params (dict): Dictionary mapping each transect file's path (keys) to a list containing the indices of selected/clicked/tapped transects (values) from its transect file
        """
        # print("Selection1D stream's parameter:", params)
        with pn.param.set_values(self._data_map_plot, loading = True):
            # Set names for the longitude and latitude data columns in the popup modal's data table.
            long_col_name = "Easting (meters)"
            lat_col_name = "Northing (meters)"
            # Save information about the recently clicked transect(s) in a dictionary.
            clicked_transects_info_dict = {}
            # Find the user's clicked/selected transect(s).
            for file_path in params:
                _, filename = os.path.split(file_path)
                if (filename in self._all_transect_files) and params[file_path]:
                    clicked_transect_indices = params[file_path]
                    num_clicked_transects = len(clicked_transect_indices)
                    # Transform the transect's coordinates into a CRS with meters as a unit.
                    transect_file_geodataframe = gpd.read_file(filename = file_path)
                    transect_crs, transect_geodataframe_crs = self._collection_crs, transect_file_geodataframe.crs
                    if transect_geodataframe_crs is not None:
                        geojson_epsg_code = ccrs.CRS(transect_geodataframe_crs).to_epsg()
                        if geojson_epsg_code == 4326: transect_file_geodataframe = transect_file_geodataframe.to_crs(crs = self._collection_crs)
                        else: transect_crs = ccrs.epsg(geojson_epsg_code)
                    else:
                        transect_file_geodataframe = transect_file_geodataframe.set_crs(crs = self._collection_crs)
                    # Reset transect file's Selection1D stream parameter to its default value (empty list []).
                    self._tapped_data_streams[file_path].reset()
                    # Add information about the clicked transect(s).
                    clicked_transects_info_dict[self._clicked_transects_file_key] = filename
                    clicked_transects_info_dict[self._num_clicked_transects_key] = num_clicked_transects
                    clicked_transects_info_dict[self._clicked_transects_crs_key] = transect_crs
                    clicked_transects_info_dict[self._clicked_transects_longitude_key] = long_col_name
                    clicked_transects_info_dict[self._clicked_transects_latitude_key] = lat_col_name
                    clicked_transects_info_dict[self._clicked_transects_id_key] = self._transects_id_col_name
                    # Specify the names of columns to display in the popup modal's data table.
                    clicked_transects_info_dict[self._clicked_transects_data_cols_key] = [long_col_name, lat_col_name, self._transects_id_col_name]
                    clicked_transects_info_dict[long_col_name] = []
                    clicked_transects_info_dict[lat_col_name] = []
                    clicked_transects_info_dict[self._transects_id_col_name] = []
                    # Get all the transects/paths from the clicked transect's file.
                    transects_file_plot = self._created_plots[file_path]
                    transect_file_paths = transects_file_plot.split()
                    # Get data for each of the user's clicked transect(s).
                    for transect_index in clicked_transect_indices:
                        path_plot = transect_file_paths[transect_index]
                        path_info = path_plot.columns(dimensions = [self._transects_id_col_name])
                        transect_id_col_vals = path_info[self._transects_id_col_name].tolist()
                        clicked_transects_info_dict[self._transects_id_col_name].extend(transect_id_col_vals)
                        transect_id = list(set(transect_id_col_vals))[0]
                        # Make sure to get each transect's easting and northing (meters) coordinates because the time-series calculations only work with non-negative values.
                        transect_geodataframe = transect_file_geodataframe[transect_file_geodataframe[self._transects_id_col_name] == transect_id]
                        transect_geojson = json.loads(transect_geodataframe.to_json())
                        transect_points = transect_geojson["features"][0]["geometry"]["coordinates"]
                        clicked_transects_info_dict[long_col_name].extend([point[0] for point in transect_points])
                        clicked_transects_info_dict[lat_col_name].extend([point[1] for point in transect_points])
                    # Stop iterating through all the transect files once a clicked transect is found.
                    break
            # Update the clicked_transects_info parameter in order to update the time-series plot, transect data table, or error message in the popup modal.
            self.clicked_transects_info = clicked_transects_info_dict
            # Reset the clicked_transects_info parameter in case the user clicked on the same transect again.
            # ^ If the parameter isn't reset, then the parameter value stays the same if the same transect is consecutively clicked more than once,
            #   meaning the info won't be sent to the modal (by Application class's _update_clicked_transects_info() method) and the modal won't open.
            self.clicked_transects_info = {}

    @param.depends("view_user_transect_time_series", watch = True)
    def _get_user_drawn_transect_info_as_clicked_transect(self) -> None:
        """
        Gets information about the user-drawn transect, which is used to update the popup modal's contents (time-series plot and transect data table).
        The view_user_transect_time_series event parameter triggers this method, acting as if the user-drawn transect was clicked.
        Bokeh doesn't allow the Tap and PolyDraw tools to be active at the same time, so an event button will be used to call this method instead of a click from the user.
        """
        with pn.param.set_values(self._data_map_plot, loading = True):
            # Save points of the most recently drawn user transect.
            longitude_col_vals = self._edit_user_transect_stream.data["xs"][0]
            latitude_col_vals = self._edit_user_transect_stream.data["ys"][0]
            longitude_col_name = "Longitude"
            latitude_col_name = "Latitude"
            self._user_transect_plot.data = [{longitude_col_name: longitude_col_vals, latitude_col_name: latitude_col_vals}]
            # Get information about the user-drawn transect as if the user-drawn transect was clicked.
            user_transect_info_dict = {
                self._clicked_transects_file_key: "User-Drawn Transect",
                self._num_clicked_transects_key: 1,
                self._clicked_transects_longitude_key: longitude_col_name,
                self._clicked_transects_latitude_key: latitude_col_name,
                self._clicked_transects_data_cols_key: [self._transects_id_col_name, longitude_col_name, latitude_col_name],
                longitude_col_name: longitude_col_vals,
                latitude_col_name: latitude_col_vals
            }
            num_points_in_user_transect = len(longitude_col_vals)
            user_transect_info_dict[self._transects_id_col_name] = [0] * num_points_in_user_transect
            # Update the clicked_transects_info parameter in order to update the time-series plot, transect data table, or error message in the popup modal.
            self.clicked_transects_info = user_transect_info_dict
            # Reset the clicked_transects_info parameter in case the user wants to view the time-series for the user-drawn transect again.
            # ^ If the parameter isn't reset, then the parameter value stays the same if the "View Time-Series for Drawn Transect" button is consecutively clicked more than once,
            #   meaning the info won't be sent to the modal (by Application class's _update_clicked_transects_info() method) and the modal won't open.
            self.clicked_transects_info = {}
    
    def _get_user_drawn_transect_geojson(self) -> dict:
        """
        Gets the user-drawn transect's points by storing them in a GeoJSON file object.
        """
        with pn.param.set_values(self._user_drawn_transect_download_button, loading = True):
            transect_points = list(zip(
                self._edit_user_transect_stream.data["xs"][0],
                self._edit_user_transect_stream.data["ys"][0],
                strict = True
            ))
            start_point = transect_points[0]
            end_point = transect_points[-1]
            # Create a GeoDataFrame from the user-drawn points stored in the PolyDraw stream.
            geodataframe = gpd.GeoDataFrame(
                data = {
                    "Type": ["User-Drawn Transect"],
                    "Date Created": [datetime.now().strftime("%B %d, %Y %H:%M:%S")],
                    self._transect_start_point_prop_name: ["({}, {})".format(start_point[0], start_point[1])],
                    self._transect_end_point_prop_name: ["({}, {})".format(end_point[0], end_point[1])],
                    self._transects_id_col_name: [0],
                    "geometry": [LineString(transect_points)]
                },
                crs = ccrs.PlateCarree()
            )
            # Convert the GeoDataFrame into a GeoJSON string with the CRS specified.
            geojson_dict = json.loads(geodataframe.to_json())
            del geojson_dict["features"][0]["id"]
            geojson_str = json.dumps(geojson_dict)
            self._user_drawn_transect_download_button.disabled = True
            # Convert the GeoJSON string into a file object and return it for downloading.
            return BytesIO(geojson_str.encode())

    def _set_ability_to_save_user_transect(self, data: dict) -> None:
        """
        Either enables or disables the button for downloading the user-drawn transect depending on the number of transect points added.

        Args:
            data (dict): Dictionary containing lists of coordinates for the "xs"/longitude and "ys"/latitude of each user-drawn transect point
        """
        # Only enable the download button when the user plotted at least 2 transect points on the map.
        if (len(data["xs"]) and len(set(data["xs"][0])) >= 2) and (len(data["ys"]) and len(set(data["ys"][0])) >= 2):
            self._user_drawn_transect_download_button.disabled = False
            self._view_user_transect_time_series_button.disabled = False
        else:
            self._user_drawn_transect_download_button.disabled = True
            self._view_user_transect_time_series_button.disabled = True
            if (not len(data["xs"])) and (not len(data["ys"])): self._user_transect_plot.data = []

    @param.depends("basemap", watch = True)
    def _update_basemap_plot(self) -> None:
        """
        Creates basemap WMTS (web mapping tile source) plot whenever the selected basemap name changes.
        """
        # Get the name of the newly selected basemap.
        if self.basemap is None:
            selected_basemap_name = list(self._all_basemaps.keys())[0]
        else:
            selected_basemap_name = self.basemap
        # Create the plot containing the basemap.
        new_basemap_plot = self._all_basemaps[selected_basemap_name]
        # Save basemap plot.
        self._selected_basemap_plot = new_basemap_plot

    @param.depends("collection", watch = True)
    def _update_collection_objects(self) -> None:
        """
        Updates everything (internal class properties, widgets, plots, etc.) that depends on the collection parameter whenever the parameter changes.
        """
        if self.collection is None: self.collection = list(self._all_collections.values())[0]
        self._collection_dir_path = os.path.join(self._root_data_dir_path, self.collection)
        collection_info_json_path = os.path.join(self._collection_dir_path, "collection_info.json")  # name of the JSON file should be same as `outputted_collection_json_name` in utils/preprocess_data.py
        if os.path.exists(collection_info_json_path):
            # Get the CRS of the new collection.
            json_file = open(collection_info_json_path)
            self._selected_collection_info = json.load(json_file)
            collection_epsg_code = self._selected_collection_info.get("epsg", 4326)   # name of the key should be same as `collection_epsg_property` in utils/preprocess_data.py
            self._collection_crs = ccrs.epsg(collection_epsg_code)
            # Get all data files' widget option names (i.e. data file names) from collection directory.
            i, self._data_file_options_dict = 0, {}
            collection_subdirs = [file for file in os.listdir(self._collection_dir_path) if os.path.isdir(os.path.join(self._collection_dir_path, file)) and (file != self._transects_folder_name)]
            for subdir in collection_subdirs:
                subdir_path = os.path.join(self._collection_dir_path, subdir)
                for file in [file for file in os.listdir(subdir_path) if os.path.isfile(os.path.join(subdir_path, file))]:
                    data_file_path = os.path.join(subdir_path, file)
                    self._data_file_options_dict[file] = data_file_path
                    # Set styles for each data file's plot.
                    self._data_file_color[data_file_path] = self._palette2_colors[i % self._total_palette2_colors]
                    self._data_file_line_style[data_file_path] = self._curve_styles[i % self._total_line_styles]
                    self._data_file_marker[data_file_path] = self._markers[i % self._total_markers]
                    i += 1
            # Get the transect widget's new options.
            transects_dir_path = os.path.join(self._collection_dir_path, self._transects_folder_name)
            if os.path.isdir(transects_dir_path):
                self._all_transect_files = [file for file in os.listdir(transects_dir_path) if os.path.isfile(os.path.join(transects_dir_path, file))]
            else:
                self._all_transect_files = []
            self._transects_multichoice.options = self._all_transect_files + [self._create_own_transect_option]
            self._transects_multichoice.value = []
            # Reassign colors and tap streams for the new collection's transects.
            if self._all_transect_files:
                for i, transect_option in enumerate(self._transects_multichoice.options):
                    self._transect_colors[transect_option] = self._palette1_colors[i % self._total_palette1_colors]
                for file in self._all_transect_files:
                    transect_file_path = os.path.join(transects_dir_path, file)
                    self._tapped_data_streams[transect_file_path] = hv.streams.Selection1D(source = None, rename = {"index": transect_file_path})
                    # Specify a callable subscriber function that gets called whenever any transect from the file is clicked/tapped.
                    self._tapped_data_streams[transect_file_path].add_subscriber(self._get_clicked_transect_info)
            # Reset the map's plots.
            self._selected_data_plot = None
            self._selected_transects_plot = None
        else:
            self._selected_collection_info = {}
            print("Error with collection {}: Please preprocess the chosen collection with `preprocess_data.py`.".format(self.collection))

    @param.depends("data_file_path", watch = True)
    def _update_selected_data_plot(self) -> None:
        """
        Creates a point or image plot whenever the most recently selected time-series data file changes.
        """
        # Only when a time-series data file is selected...
        if self.data_file_path is not None:
            # Create a plot with data from the file.
            new_data_plot = None
            if self.data_file_path not in self._created_plots: self._create_data_plot()
            # Display the data file's plot if it was created.
            # ^ plots aren't created for unsupported files -> e.g. PNG
            if self.data_file_path in self._created_plots: new_data_plot = self._created_plots[self.data_file_path]
            # Save the data plot.
            self._selected_data_plot = new_data_plot
        else:
            self._selected_data_plot = None

    @param.depends("transects", watch = True)
    def _update_selected_transects_plot(self) -> None:
        """
        Creates an overlay of path plots whenever the selected transect files change.
        """
        # Only when the widget is initialized and at least one transect file is selected...
        if self.transects is not None:
            # Create an overlay of path plots with transects from each selected transect file.
            new_transects_plot = None
            # If the user wants to create their own transect, display buttons related to the user-drawn transect.
            if self._create_own_transect_option in self.transects:
                self._view_user_transect_time_series_button.visible = True
                self._user_drawn_transect_download_button.visible = True
                self._drawing_user_transect_instructions.visible = True
            else:
                self._view_user_transect_time_series_button.visible = False
                self._user_drawn_transect_download_button.visible = False
                self._drawing_user_transect_instructions.visible = False
            for file in self.transects:
                file_path = os.path.join(self._collection_dir_path, self._transects_folder_name, file)
                # Allow user to draw start and end points when they selected to draw their own transect.
                if file == self._create_own_transect_option:
                    # Display an editable path plot for the user to modify their transect's start and end points.
                    if new_transects_plot is None:
                        new_transects_plot = self._user_transect_plot
                    else:
                        new_transects_plot = (new_transects_plot * self._user_transect_plot)
                else:
                    # Create the selected transect file's path plot if we never read the file before.
                    if file_path not in self._created_plots:
                        self._create_path_plot(file)
                        # Save the new plot as a source for the transect file's Selection1D stream.
                        if file_path in self._created_plots: self._tapped_data_streams[file_path].source = self._created_plots[file_path]
                    # Display the transect file's path plot if it was created.
                    # ^ plots aren't created for unsupported files
                    if file_path in self._created_plots:
                        if new_transects_plot is None:
                            new_transects_plot = self._created_plots[file_path]
                        else:
                            new_transects_plot = (new_transects_plot * self._created_plots[file_path])
            # Save overlaid transect plots.
            self._selected_transects_plot = new_transects_plot
    
    def _update_map_data_ranges(self, plot: any, element: any) -> None:
        """
        Fixes the map's data range when the user selected the option to create their own transect
        because an empty path plot causes the map to automatically zoom in to the middle of the map.

        Args:
            plot (any): HoloViews object rendering the plot; this hook/method is applied after the plot is rendered
            element (any): Element rendered in the plot
        """
        # print(plot, element)
        # print('plot.state:   ', plot.state)
        # print('plot.handles: ', sorted(plot.handles.keys()))
        # print('plot.handles.x_range: ', plot.handles['x_range'])
        # print("plot.handles.x_range dict:", plot.handles['x_range'].__dict__)
        # print(plot.handles['x_range'].start)
        # print(plot.handles['x_range'].end)
        # print("plot.handles.y_range dict:", plot.handles['y_range'].__dict__)
        if (self.transects is not None) and (len(self.transects) == 1) and (self._create_own_transect_option in self.transects) and (not len(self._user_transect_plot.data)) and (self.data_file_path is None):
            plot.handles["x_range"].start = plot.handles["x_range"].reset_start = -20037508.342789244
            plot.handles["x_range"].end = plot.handles["x_range"].reset_end = 20037508.342789244
            plot.handles["y_range"].start = plot.handles["y_range"].reset_start = -20037508.342789248
            plot.handles["y_range"].end = plot.handles["y_range"].reset_end = 20037508.342789248

    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @param.depends("_update_basemap_plot", "_update_collection_objects", "_update_selected_data_plot", "_update_selected_transects_plot", "_get_clicked_transect_info")
    def plot(self) -> gv.Overlay:
        """
        Returns the selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        # Overlay the selected plots.
        current_active_tools = ["pan", "wheel_zoom"]
        new_plot = self._selected_basemap_plot
        if self._selected_data_plot is not None:
            new_plot = (new_plot * self._selected_data_plot)
        if self._selected_transects_plot is not None:
            new_plot = (new_plot * self._selected_transects_plot)
            if self._create_own_transect_option in self.transects:
                current_active_tools.append("poly_draw")
        # Save the overlaid plots.
        self._data_map_plot.object = new_plot.opts(
            xaxis = None, yaxis = None,
            tools = ["zoom_in", "zoom_out", "tap"],
            active_tools = current_active_tools,
            toolbar = "below",#None,"above"
            title = "", #show_legend = True,
            hooks = [self._update_map_data_ranges]
        )
        # Display browser popup for any errors that occurred while updating the data map.
        if self._error_messages:
            self._error_popup_text.value = "\n".join(self._error_messages)
            # Reset the error messages to an empty list in order to indicate that there are no errors by default.
            self._error_messages = []
        return self._data_map_plot

    @property
    def param_widgets(self) -> list[any]:
        """
        Returns a list of parameters (will have default widget) or custom Panel widgets for parameters used in the app.
        """
        widgets = [
            self.param.basemap,
            self._collection_select,
            self._transects_multichoice,
            self._error_popup_text,
            self._view_user_transect_time_series_button,
            self._user_drawn_transect_download_button,
            self._drawing_user_transect_instructions
        ]
        return widgets

    @property
    def app_main_color(self) -> str:
        """
        Returns the main color of the application (used for navigation bar, hovered and clicked transects, etc.).
        """
        return self._app_main_color

    @property
    def selected_basemap(self) -> any:
        """
        Returns the user's selected basemap.
        """
        return self._selected_basemap_plot
    
    @property
    def selected_collection_dir_path(self) -> str:
        """
        Returns the directory path to the user's selected collection.
        """
        return self._collection_dir_path
    
    @property
    def selected_collection_crs(self) -> ccrs:
        """
        Returns the selected collection's coordinate reference system.
        """
        return self._collection_crs

    @property
    def selected_collection_json_info(self) -> ccrs:
        """
        Returns the selected collection's information from its `collection_info.json` file.
        """
        return self._selected_collection_info
    
    @property
    def data_file_options(self) -> ccrs:
        """
        Returns the dictionary mapping the option name (key) of each data file in the selected collection to the data file's path (value).
        """
        return self._data_file_options_dict
    
    @property
    def transects_dir_name(self) -> str:
        """
        Returns the name of the directory that is reserved for storing transects.
        """
        return self._transects_folder_name
    
    @property
    def map_default_crs(self) -> ccrs:
        """
        Returns the data map's default coordinate reference system.
        """
        return self._default_crs

    @property
    def data_file_color(self) -> dict:
        """
        Returns the dictionary that maps the path to each data file in the selected collection to a plot color.
        """
        return self._data_file_color
    
    @property
    def data_file_marker(self) -> dict:
        """
        Returns the dictionary that maps the path to each data file in the selected collection to a plot marker.
        """
        return self._data_file_marker
    
    @property
    def data_file_line_style(self) -> dict:
        """
        Returns the dictionary that maps the path to each data file in the selected collection to a plot line style.
        """
        return self._data_file_line_style
    
    @property
    def error_messages(self) -> list:
        """
        Returns the list containing error messages to display in a browser popup window if there's any errors while creating plots for the data map or popup modal.
        """
        return self._error_messages

    @property
    def clicked_transects_info_keys(self) -> list[str]:
        """
        Returns a list of keys that appear in the clicked_transects_info parameter's dictionary,
        which is used by the popup modal to display information about the user's clicked transect(s) from this map.
        """
        return [
            self._clicked_transects_file_key,
            self._num_clicked_transects_key,
            self._clicked_transects_crs_key,
            self._clicked_transects_longitude_key,
            self._clicked_transects_latitude_key,
            self._clicked_transects_data_cols_key,
            self._clicked_transects_id_key
        ]