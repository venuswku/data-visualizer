# Standard library importspn
import os
import json
from datetime import datetime

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
import holoviews as hv
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import LineString
import rioxarray as rxr
from bokeh.palettes import Bokeh
from io import BytesIO

### DataMap is used for displaying the inputted data files onto a map. ###
class DataMap(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    basemap = param.Selector(label = "Basemap")
    categories = param.ListSelector(label = "Data Categories")
    transects = param.ListSelector(label = "Transects")
    clicked_transects_info = param.Dict(default = {}, label = "Information About the Recently Clicked Transect(s)")
    view_user_transect_time_series = param.Event(label = "Indicator for Displaying the Time-Series for Data Along the User-Drawn Transect")

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, data_dir_path: str, latitude_col_names: list[str], longitude_col_names: list[str], colors: dict = {}, basemap_options: dict = {"Default": gts.OSM}, **params) -> None:
        """
        Creates a new instance of the DataMap class with its instance variables.

        Args:
            data_dir_path (str): Path to the directory containing all the data category subfolders and their data files
            latitude_col_names (list[str]): Possible names of the column containing the latitude of each data point
            longitude_col_names (list[str]): Possible names of the column containing the longitude of each data point
            colors (dict): Optional dictionary mapping each data category name (keys) to a color (values), which will be the color of its data points
            basemap_options (dict): Optional dictionary mapping each basemap name (keys) to a basemap WMTS (web mapping tile source) layer (values)
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        self._data_dir_path = data_dir_path
        self._all_lat_cols = latitude_col_names
        self._all_long_cols = longitude_col_names
        self._app_main_color = "#2196f3"
        
        # _transects_folder_name = Name of the folder containing files with transect data
        self._transects_folder_name = "Transects"
        # _create_own_transect_option = Name of the option for the user to create their own transect
        self._create_own_transect_option = "Draw My Own Transect"
        # _geodata_folder_name = Name of the folder containing GeoJSON/GeoTIFF files that were created by georeferencing data files (txt, csv, asc)
        # ^ allows data to load faster onto the map
        self._geodata_folder_name = "GeoData"
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
        # _transects_id_col_name = Name of the column containing the ID of each clicked transect
        # ^ initially created in _convert_transect_data_into_geojson() when assigning an ID property for each transect in the outputted GeoJSON
        self._transects_id_col_name = "Transect ID"
        # _transect_start_point_prop_name = Name of the GeoJSON property containing the start point of a transect (in meters)
        self._transect_start_point_prop_name = "Start Point"
        # _transect_end_point_prop_name = Name of the GeoJSON property containing the end point of a transect (in meters)
        self._transect_end_point_prop_name = "End Point"
        
        # _default_crs = default coordinate reference system for the user-drawn transect and other plots
        self._default_crs = ccrs.PlateCarree()
        # _crs = custom coordinate reference system for the projected data (e.g. GeoTIFFs)
        # ^ can be created from a dictionary of PROJ parameters
        # ^ https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.CRS.html#cartopy.crs.CRS.__init__
        self._crs = ccrs.CRS({
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
        # _epsg_code = number code that corresponds to the _crs specified above
        self._epsg_code = 32148
        # _epsg = publicly registered coordinate system for the projected data
        # ^ should be close, if not equivalent, to the custom CRS defined above (_crs)
        self._epsg = ccrs.epsg(self._epsg_code)

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _data_map_plot = overlay plot containing the selected basemap and all the data (categories, transects, etc.) plots
        self._data_map_plot = pn.pane.HoloViews(object = None, sizing_mode = "scale_both")
        # _created_plots = dictionary mapping each file's path (keys) to its created plot (values)
        self._created_plots = {}
        
        # _all_basemaps = dictionary mapping each basemap name (keys) to a basemap WMTS (web mapping tile source) layer (values)
        self._all_basemaps = basemap_options
        # _selected_basemap_plot = WMTS (web mapping tile source) layer containing the user's selected basemap
        self._selected_basemap_plot = basemap_options[list(basemap_options.keys())[0]]
        
        # _all_categories = list of data categories (subfolders in the root data directory -> excludes transects)
        self._all_categories = [file for file in os.listdir(data_dir_path) if os.path.isdir(data_dir_path + "/" + file) and (file != self._transects_folder_name)]
        # _category_colors = dictionary mapping each data category (keys) to a color (values), which will be used for the color of its point plots
        self._category_colors = {}
        # _category_markers = dictionary mapping each data category (keys) to a marker (values), which will be used for the marker of its point plots
        self._category_markers = {}
        # _selected_categories_plot = overlay of point plots for each data category selected by the user
        # ^ None if the no data files were provided by the user or the user didn't select any data categories
        self._selected_categories_plot = None
        
        # _all_transect_files = list of files containing transects to display on the map
        self._all_transect_files = []
        transects_dir_path = data_dir_path + "/" + self._transects_folder_name
        if os.path.isdir(transects_dir_path):
            self._all_transect_files = [file for file in os.listdir(transects_dir_path) if os.path.isfile(os.path.join(transects_dir_path, file))]
        # _transect_colors = dictionary mapping each transect file (keys) to a color (values), which will be used for the color of its path plots
        self._transect_colors = {}
        # _selected_transects_plot = overlay of path plots if the user selected one or more transect files to display on the map
        # ^ None if the user didn't provide any transect files or the user didn't select to display a transect file
        self._selected_transects_plot = None

        # _tapped_data_streams = dictionary mapping each transect file's path (keys) to a selection stream (values), which saves the file's most recently clicked data element (path) on the map
        self._tapped_data_streams = {}
        for file in self._all_transect_files:
            file_path = transects_dir_path + "/" + file
            self._tapped_data_streams[file_path] = hv.streams.Selection1D(source = None, rename = {"index": file_path})
            # Specify a callable subscriber function that gets called whenever any transect from the file is clicked/tapped.
            self._tapped_data_streams[file_path].add_subscriber(self._get_clicked_transect_info)

        # _user_transect_plot = path plot used when the user wants to create their own transect to display on the map
        self._user_transect_plot = gv.Path(data = [], crs = self._default_crs)#.opts(projection = self._default_crs)
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
        self.param.basemap.objects = basemap_options.keys()

        # Set data category widget's options.
        self._categories_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.categories,
            options = self._all_categories,
            placeholder = "Choose one or more data categories to display",
            solid = False
        )

        # Set transect widget's options.
        self._transects_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.transects,
            options = self._all_transect_files + [self._create_own_transect_option],
            placeholder = "Choose one or more transect files to display",
            solid = False
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
                            <li>The Pan tool <img src="https://raw.githubusercontent.com/venuswku/data-visualizer/draw-own-transect/assets/PanTool.png" alt="Pan tool" width="16"/> allows you to move around the map by holding down your mouse/trackpad and dragging.</li>
                            <li>
                                The Pan tool is automatically disabled whenever the Polygon Draw tool <img src="https://raw.githubusercontent.com/venuswku/data-visualizer/draw-own-transect/assets/PolygonDrawTool.png" alt="Polygon Draw tool" width="16"/> is enabled (allowing you to draw your own transect), and vice versa. 
                                Click on either tool to toggle between them.
                            </li>
                        </ul>
                    </li>
                </ul>
            """, alert_type = "primary"),
            pn.pane.Markdown(object = """
                <details>
                    <summary><b>How to Draw Your Own Transects</b></summary>
                    <div style="padding-left:14px">
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
                    </div>
                </details>
            """, sizing_mode = "stretch_width"
            ), visible = False, margin = (0, 5, 5, 10)
        )

        # Set color and marker for each data category.
        palette_colors = Bokeh[8]
        markers = ["o", "^", "s", "d", "x", ">", "*", "v", "+", "<"]
        total_palette_colors, total_markers = len(palette_colors), len(markers)
        for i, category in enumerate(self._all_categories):
            # Assign the color that the user chose, if provided.
            if category in colors:
                self._category_colors[category] = colors[category]
            # Else assign a color from the Bokeh palette.
            else:
                self._category_colors[category] = palette_colors[i % total_palette_colors]
            # Assign a marker.
            self._category_markers[category] = markers[i % total_markers]
        # Set color for each transect option.
        for i, transect_option in enumerate(self._transects_multichoice.options):
            self._transect_colors[transect_option] = palette_colors[(len(category) + i) % total_palette_colors]

    # -------------------------------------------------- Private Class Methods --------------------------------------------------    
    def _plot_geojson_points(self, geojson_file_path: str, data_category: str) -> gv.Points:
        """
        Creates a point plot from a GeoJSON file containing Points.

        Args:
            geojson_file_path (str): Path to the GeoJSON file containing Points
            data_category (str): Name of the data category that the data file belongs to
        """
        # Read the GeoJSON as a GeoDataFrame.
        geodataframe = gpd.read_file(geojson_file_path)
        latitude_col, longitude_col, non_lat_long_cols = None, None, []
        for col in geodataframe.columns:
            if col in self._all_lat_cols: latitude_col = col
            elif col in self._all_long_cols: longitude_col = col
            elif col != "geometry": non_lat_long_cols.append(col)
        # Create a point plot with the GeoDataFrame.
        point_plot = gv.Points(
            data = geodataframe,
            kdims = [longitude_col, latitude_col],
            vdims = non_lat_long_cols,
            label = data_category
        ).opts(
            color = self._category_colors[data_category],
            marker = self._category_markers[data_category],
            hover_color = self._app_main_color,
            tools = ["hover"],
            size = 10, muted_alpha = 0.01
        )
        return point_plot
    
    def _plot_geojson_linestrings(self, geojson_file_path: str, filename: str) -> gv.Path:
        """
        Creates a path plot from a GeoJSON file containing LineStrings.

        Args:
            geojson_file_path (str): Path to the GeoJSON file containing LineStrings
            filename (str): Name of the transect file that corresponds to the returned path plot
        """
        # plot_crs = self._default_crs
        # Convert the CRS from the GeoJSON file into GeoView's default CRS (Plate Carree), if necessary.
        geodataframe = gpd.read_file(geojson_file_path)
        geojson_crs = geodataframe.crs
        if geojson_crs is not None:
            geojson_epsg_code = ccrs.CRS(geojson_crs).to_epsg()
            # Only use projected coordinate systems, not geodetic coordinate systems like EPSG:4326/WGS-84 (https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.epsg.html).
            if geojson_epsg_code not in [4326]:
                # plot_crs = ccrs.epsg(geojson_epsg_code)
                geodataframe = geodataframe.to_crs(crs = self._default_crs)
        # Create a path plot with the correct CRS.
        return gv.Path(
            data = geodataframe,
            crs = self._default_crs,#plot_crs,
            label = "{}: {}".format(self._transects_folder_name, filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
        ).opts(
            color = self._transect_colors[filename],
            hover_color = self._app_main_color,
            selection_color = self._app_main_color,
            nonselection_color = self._transect_colors[filename],
            nonselection_alpha = 1, selection_alpha = 1,
            tools = ["hover", "tap"]
        )
    
    def _create_data_plot(self, filename: str, category: str) -> None:
        """
        Creates a point/image plot containing the given file's data.

        Args:
            filename (str): Name of the file containing data
            category (str): Name of the data category that the file belongs to
        """
        # Read the file and create a plot from it.
        file_path = self._data_dir_path + "/" + category + "/" + filename
        [name, extension] = os.path.splitext(filename)
        extension = extension.lower()
        category_geodata_dir_path = self._data_dir_path + "/" + category + "/" + self._geodata_folder_name
        plot = None
        if extension in [".csv", ".txt"]:
            # Convert the data file into a new GeoJSON (if not created yet).
            geojson_path = category_geodata_dir_path + "/" + name + ".geojson"
            self.convert_csv_txt_data_into_geojson(file_path, geojson_path)
            # Create a point plot with the GeoJSON.
            plot = self._plot_geojson_points(geojson_path, category)
        elif extension == ".asc":
            # Convert ASCII grid file into a new GeoTIFF (if not created yet).
            geotiff_path = category_geodata_dir_path + "/" + name + ".tif"
            self.convert_ascii_grid_data_into_geotiff(file_path, geotiff_path)
            # Create an image plot with the GeoTIFF.
            plot = gv.load_tiff(
                geotiff_path,
                vdims = "Elevation (meters)",
                nan_nodata = True
            ).opts(
                cmap = "Turbo",
                tools = ["hover"],
                alpha = 0.5
            )
        if plot is None:
            print("Error displaying", name + extension, "as a point/image plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            # Save the created plot.
            self._created_plots[file_path] = plot
    
    def _create_path_plot(self, filename: str) -> None:
        """
        Creates a path plot containing the given file's paths.

        Args:
            filename (str): Name of the file containing paths
        """
        # Read the given file.
        file_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + filename
        [name, extension] = os.path.splitext(filename)
        extension = extension.lower()
        transects_geodata_dir_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + self._geodata_folder_name
        # Create a path plot from the given file.
        plot = None
        if extension == ".txt":
            geojson_path = transects_geodata_dir_path + "/" + name + ".geojson"
            # Get the GeoJSON for the given path file.
            self._convert_transect_data_into_geojson(file_path, geojson_path)
            # Create a path plot with the path GeoJSON.
            plot = self._plot_geojson_linestrings(geojson_path, filename)
        elif extension == ".geojson":
            plot = self._plot_geojson_linestrings(file_path, filename)
        # Save the path plot, if created.
        if plot is None:
            print("Error displaying", name + extension, "as a path plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            self._created_plots[file_path] = plot

    def _convert_transect_data_into_geojson(self, file_path: str, geojson_path: str) -> None:
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
                    [name, extension] = os.path.splitext(filename)
                    if extension == ".geojson": transect_file_geodataframe = gpd.read_file(filename = file_path)
                    else: transect_file_geodataframe = gpd.read_file(filename = self._data_dir_path + "/" + self._transects_folder_name + "/" + self._geodata_folder_name + "/" + name + ".geojson")
                    transect_crs, transect_geodataframe_crs = self._epsg, transect_file_geodataframe.crs
                    if transect_geodataframe_crs is not None:
                        geojson_epsg_code = ccrs.CRS(transect_geodataframe_crs).to_epsg()
                        if geojson_epsg_code == 4326: transect_file_geodataframe = transect_file_geodataframe.set_crs(crs = self._default_crs, allow_override = True).to_crs(crs = self._epsg)
                        else: transect_crs = ccrs.epsg(geojson_epsg_code)
                    else:
                        transect_file_geodataframe = transect_file_geodataframe.set_crs(crs = self._epsg)
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
                        path_info = transect_file_paths[transect_index].columns(dimensions = [self._transects_id_col_name])
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
            self._user_transect_plot.data = [{"Longitude": longitude_col_vals, "Latitude": latitude_col_vals}]
            # Transform the transect's coordinates into a CRS with meters as a unit.
            easting_col_name = "Easting (meters)"
            northing_col_name = "Northing (meters)"
            # transformed_points = self._epsg.transform_points(src_crs = self._default_crs, x = longitude_col_vals, y = latitude_col_vals)
            # easting_col_vals = [point[0] for point in transformed_points]
            # northing_col_vals = [point[1] for point in transformed_points]
            # Get information about the user-drawn transect as if the user-drawn transect was clicked.
            user_transect_info_dict = {
                self._clicked_transects_file_key: "User-Drawn Transect",
                self._num_clicked_transects_key: 1,
                # self._clicked_transects_crs_key: self._epsg,
                self._clicked_transects_longitude_key: easting_col_name,
                self._clicked_transects_latitude_key: northing_col_name,
                self._clicked_transects_data_cols_key: [self._transects_id_col_name, easting_col_name, northing_col_name],
                easting_col_name: longitude_col_vals,#easting_col_vals
                northing_col_name: latitude_col_vals#northing_col_vals
            }
            num_points_in_user_transect = len(longitude_col_vals)
            user_transect_info_dict[self._transects_id_col_name] = [0] * num_points_in_user_transect
            # Update the clicked_transects_info parameter in order to update the time-series plot, transect data table, or error message in the popup modal.
            self.clicked_transects_info = user_transect_info_dict
            # print(user_transect_info_dict)
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

    @param.depends("categories", watch = True)
    def _update_selected_categories_plot(self) -> None:
        """
        Creates an overlay of data plots whenever the selected data categories change.
        """
        # Only when the widget is initialized and at least one data category is selected...
        if self.categories is not None:
            # Create a plot with data from each selected category.
            new_data_plot = None
            selected_category_names = self.categories
            for category in self._all_categories:
                category_dir_path = self._data_dir_path + "/" + category
                category_files = [file for file in os.listdir(category_dir_path) if os.path.isfile(os.path.join(category_dir_path, file))]
                for file in category_files:
                    if category in selected_category_names:
                        # Create the selected data's point plot if we never read the file before.
                        file_path = category_dir_path + "/" + file
                        if file_path not in self._created_plots:
                            self._create_data_plot(file, category)
                        # Display the data file's point plot if it was created.
                        # ^ plots aren't created for unsupported files -> e.g. png files don't have data points
                        if file_path in self._created_plots:
                            if new_data_plot is None:
                                new_data_plot = self._created_plots[file_path]
                            else:
                                new_data_plot = (new_data_plot * self._created_plots[file_path])
            # Save overlaid category plots.
            self._selected_categories_plot = new_data_plot

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
                file_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + file
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
                        self._tapped_data_streams[file_path].source = self._created_plots[file_path]
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
        if (self.transects is not None) and (len(self.transects) == 1) and (self._create_own_transect_option in self.transects) and (not len(self._user_transect_plot.data)):
            if (self.categories is None) or ((self.categories is not None) and (len(self.categories) == 0)):
                plot.handles["x_range"].start = plot.handles["x_range"].reset_start = -20037508.342789244
                plot.handles["x_range"].end = plot.handles["x_range"].reset_end = 20037508.342789244
                plot.handles["y_range"].start = plot.handles["y_range"].reset_start = -20037508.342789248
                plot.handles["y_range"].end = plot.handles["y_range"].reset_end = 20037508.342789248

    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @param.depends("_update_basemap_plot", "_update_selected_categories_plot", "_update_selected_transects_plot", "_get_clicked_transect_info")
    def plot(self) -> gv.Overlay:
        """
        Returns the selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        # Overlay the selected plots.
        current_active_tools = ["pan", "wheel_zoom"]
        new_plot = self._selected_basemap_plot
        if self._selected_categories_plot is not None:
            new_plot = (new_plot * self._selected_categories_plot)
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
            title = "", show_legend = True,
            hooks = [self._update_map_data_ranges]
        )
        return self._data_map_plot

    @property
    def param_widgets(self) -> list[any]:
        """
        Returns a list of parameters (will have default widget) or custom Panel widgets for parameters used in the app.
        """
        widgets = [
            self.param.basemap,
            self._categories_multichoice,
            self._transects_multichoice,
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
    def geodata_dir(self) -> str:
        """
        Returns name of the directory containing GeoJSON/GeoTIFF files that were created by georeferencing data files (txt, csv, asc).
        """
        return self._geodata_folder_name

    @property
    def selected_basemap(self) -> any:
        """
        Returns the user's selected basemap.
        """
        return self._selected_basemap_plot
    
    @property
    def map_default_crs(self) -> ccrs:
        """
        Returns the data map's default coordinate reference system.
        """
        return self._default_crs

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
    
    def convert_csv_txt_data_into_geojson(self, file_path: str, geojson_path: str) -> None:
        """
        Creates and saves a GeoJSON file containing Points for each data point in the given dataframe.
        
        Args:
            file_path (str): Path to the file containing data points
            geojson_path (str): Path to the newly created GeoJSON file, which is a FeatureCollection of Points
        """
        if not os.path.exists(geojson_path):
            # Create the GeoData folder if it doesn't exist yet.
            geodata_dir_path, _ = os.path.split(geojson_path)
            if not os.path.isdir(geodata_dir_path): os.makedirs(geodata_dir_path)
            # Read the data file as a DataFrame and replace any NaN values with "N/A".
            dataframe = pd.read_csv(file_path).fillna("N/A")
            # Ignore any unnamed columns.
            dataframe = dataframe.loc[:, ~dataframe.columns.str.match("Unnamed")]
            [latitude_col] = [col for col in dataframe.columns if col in self._all_lat_cols]
            [longitude_col] = [col for col in dataframe.columns if col in self._all_long_cols]
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

    def convert_ascii_grid_data_into_geotiff(self, file_path: str, geotiff_path: str) -> None:
        """
        Converts an ASCII grid file into a GeoTIFF file.

        Args:
            file_path (str): Path to the ASCII grid file
            geotiff_path (str): Path to the newly created GeoTIFF file
        """
        if not os.path.exists(geotiff_path):
            # Create the GeoData folder if it doesn't exist yet.
            geodata_dir_path, _ = os.path.split(geotiff_path)
            if not os.path.isdir(geodata_dir_path): os.makedirs(geodata_dir_path)
            # Add custom projection to the data file based on Elwha data's metadata.
            dataset = rxr.open_rasterio(file_path)
            dataset.rio.write_crs(self._epsg, inplace = True)
            # Save the data as a GeoTIFF.
            dataset.rio.to_raster(
                raster_path = geotiff_path,
                driver = "GTiff"
            )