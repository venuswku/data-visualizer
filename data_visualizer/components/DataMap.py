# Standard library imports
import os
import json
import asyncio
from datetime import datetime
import time

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
import holoviews as hv
from holoviews.operation.datashader import dynspread, rasterize, inspect_points, datashade
import dask_geopandas
import spatialpandas as spd
import geopandas as gpd
import rioxarray as rxr
import cartopy.crs as ccrs
from shapely.geometry import LineString
from bokeh.models import HoverTool
from bokeh.palettes import Bokeh
from io import BytesIO

### DataMap is used for displaying the inputted data files onto a map. ###
class DataMap(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    basemap = param.Selector(label = "Basemap")
    collection = param.Selector(label = "Collection")
    transects = param.ListSelector(label = "Transects")
    clicked_transects_info = param.Dict(default = {}, label = "Information About the Recently Clicked Transect(s)")
    data_file_paths = param.List(default = [], label = "List of Paths to Data Files to Display on the Map")
    
    view_user_transect_time_series = param.Event(label = "Indicator for Displaying the Time-Series for Data Along the User-Drawn Transect")
    update_accordion_section = param.Event(label = "Indicator for Updating the DataMap's Accordion Sections")

    # -------------------------------------------------- Constructor --------------------------------------------------
    def __init__(self, time_series_data: list[str] = [], **params) -> None:
        """
        Creates a new instance of the DataMap class with its instance variables.

        Args:
            time_series_data (list[str]): List of column names for columns containing data for the time-series' y-axis
        """
        super().__init__(**params)

        # -------------------------------------------------- Constants --------------------------------------------------
        # _all_time_series_data = dictionary mapping each collection (key) to a list of column names (value) containing data for the collection's time-series
        self._all_time_series_data = time_series_data
        # _root_data_dir_path = path to the root directory that contains all available datasets/collections for the app
        self._root_data_dir_path = os.path.relpath("./data")
        # _default_crs = default coordinate reference system for the user-drawn transect and other plots
        self._default_crs = ccrs.PlateCarree()
        # _app_main_color = theme color used for all the Panel widgets in this app
        self._app_main_color = "#2196f3"
        
        # _all_basemaps = dictionary mapping each basemap name (key) to a basemap WMTS (web mapping tile source) layer (value)
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
        # _palette_colors = list of palette colors for transects in the map
        self._palette_colors = Bokeh[8]
        self._total_palette_colors = len(self._palette_colors)
        
        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
        # _data_map_plot = overlay plot containing the selected basemap and all the data (categories, transects, etc.) plots
        self._data_map_plot = None
        # _created_plots = dictionary mapping each file's path (key) to its created plot (value)
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

        # _data_file_options_dict = dictionary mapping the option name (key) of each data file in the selected collection to the data file's path (value)
        self._data_file_options_dict = {}
        # _selected_data_plot = point or image plot for the most recently selected data file
        # ^ None if the no data file was selected by the user
        self._selected_data_plot = None

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
        # Show an error popup if there are any errors that occurred while creating plots for the data map.
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
        # Create an accordion section for tips and instructions on how to use the PolyDraw tool to draw a custom transect.
        self._display_user_drawn_transect_instructions = False
        self._drawing_user_transect_accordion_section = (
            "Draw Your Own Transects",
            pn.Column(
                pn.pane.Alert(object = """
                    <p style="font-weight: bold; margin: 10px 0px 0px 0px">Tips</p>
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
                """, alert_type = "primary", margin = (0, 10, 15, 10)),
                pn.pane.Markdown(object = """
                    <style>
                    .instruction-subheading { font-weight: bold; margin: 10px 0px 0px 0px }
                    </style>
                    <p class="instruction-subheading">Add</p>
                    <ul>
                        <li>Double click to add the start point.</li>
                        <li>(Optional) Single click to add each subsequent point.</li>
                        <li>If you want to restart adding transect points, press the ESC key.</li>
                        <li>Double click to add the end point.</li>
                    </ul>
                    <p class="instruction-subheading">Move</p>
                    <ul>
                        <li>Click to select an existing transect.</li>
                        <li>Then drag the transect to move it.</li>
                        <li>Transect points will be moved once you let go of the mouse/trackpad.</li>
                    </ul>
                    <p class="instruction-subheading">Delete</p>
                    <ul>
                        <li>Click to select an existing transect.</li>
                        <li>Then press the BACKSPACE (Windows) or DELETE (Mac) key while the cursor is within the map area.</li>
                    </ul>
                """, margin = (0, 10), sizing_mode = "stretch_width")
            )
        )

    # -------------------------------------------------- Private Class Methods --------------------------------------------------    
    def _plot_parquet_points(self, data_file_path: str, data_file_option: str) -> gv.Points:
        """
        Creates a point plot from Parquet partition files containing Points.

        Args:
            data_file_path (str): Path to the directory containing data to plot
            data_file_option (str): Option name of the most recently selected data file from PopupModal's _data_files_checkbox_group widget
        """
        start_time = time.time()
        collection_time_series_data_cols = self._all_time_series_data[self.collection]
        # Read the Parquet file as a geopandas GeoDataFrame.
        geodataframe = dask_geopandas.read_parquet(data_file_path).compute()
        latitude_col, longitude_col, time_series_data_col, non_lat_long_cols = None, None, None, []
        for col in geodataframe.columns:
            col_name = col.lower()
            if "lat" in col_name: latitude_col = col
            elif "lon" in col_name: longitude_col = col
            elif col in collection_time_series_data_cols: time_series_data_col = col
            elif col_name != "geometry": non_lat_long_cols.append(col)
        geodataframe = geodataframe.drop(columns = [latitude_col, longitude_col])
        mid_time = time.time()
        print("Reading parquet data took {} seconds.".format(mid_time - start_time))
        # Convert the geopandas GeoDataFrame into a spatialpandas GeoDataFrame for the geometry values to be compatible with GeoViews.
        spatial_pandas_geodataframe = spd.GeoDataFrame(geodataframe)
        # Create a point plot with the spatialpandas GeoDataFrame.
        custom_hover_tool = HoverTool(tooltips = [
            ("Longitude", "$x"),
            ("Latitude", "$y"),
            (time_series_data_col, "@image")
        ])
        point_plot = dynspread(
            rasterize(
                gv.Points(
                    data = spatial_pandas_geodataframe,
                    kdims = [longitude_col, latitude_col],
                    vdims = [time_series_data_col] + non_lat_long_cols,
                    # label = data_file_option
                ),
            ).opts(
                cmap = "Turbo", tools = [custom_hover_tool],
                cnorm = "eq_hist", responsive = True
            ),
            max_px = 5
        )
        end_time = time.time()
        print("Creating point plot took {} seconds.".format(end_time - mid_time))
        # Return the point plot.
        return point_plot
    
    def _plot_geojson_linestrings(self, geojson_file_path: str) -> gv.Path | None:
        """
        Creates a path plot from a GeoJSON file containing LineStrings. Returns None if the a transect is invalid (e.g. less than 2 points).

        Args:
            geojson_file_path (str): Path to the GeoJSON file containing LineStrings
        """
        filename = os.path.basename(geojson_file_path)
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
            self._error_messages.append("⚠️" + error)
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
                crs = self._default_crs,
                label = "Transects: {}".format(filename)
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
                label = "Transects: {}".format(filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
            ).opts(
                color = self._transect_colors[filename],
                hover_color = self._app_main_color,
                selection_color = self._app_main_color,
                nonselection_color = self._transect_colors[filename],
                nonselection_alpha = 1, selection_alpha = 1,
                tools = ["hover", "tap"]
            )
    
    def _plot_parquet_linestrings(self, parquet_file_path: str) -> gv.Path | None:
        """
        Creates a path plot from a Parquet directory of partition files containing LineStrings. Returns None if the a transect is invalid (e.g. less than 2 points).

        Args:
            parquet_file_path (str): Path to the Parquet directory containing LineStrings
        """
        filename = os.path.basename(parquet_file_path)
        # Read the Parquet file as a geopandas GeoDataFrame.
        gpd_geodataframe = dask_geopandas.read_parquet(parquet_file_path).compute()
        transect_info_cols = [col for col in gpd_geodataframe.columns if col != "geometry"]
        # Return the path plot.
        return gv.Path(
            data = spd.GeoDataFrame(gpd_geodataframe),
            vdims = transect_info_cols,
            crs = self._collection_crs,
            label = "Transects: {}".format(filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
        ).opts(
            hover_color = self._app_main_color,
            selection_color = self._app_main_color,
            nonselection_color = self._transect_colors[filename],
            nonselection_alpha = 1, selection_alpha = 1,
            tools = ["hover", "tap"], responsive = True
        )

    def _create_data_plot(self, data_file_path: str) -> None:
        """
        Creates and returns a point/image plot containing the given file's data.

        Args:
            data_file_path (str): Path to the file containing data to plot
        """
        start_time = time.time()
        # Read the file and create a plot from it.
        subdir_path, filename = os.path.split(data_file_path)
        subdir_name = os.path.basename(subdir_path)
        _, extension = os.path.splitext(filename)
        extension = extension.lower()
        plot = None
        if extension in [".parq", ".parquet"]:
            # Create a point plot with the Parquet partition files.
            plot = self._plot_parquet_points(
                data_file_path = data_file_path,
                data_file_option = self._selected_collection_info.get(data_file_path, "{}: {}".format(subdir_name, filename))
            )
        elif extension in [".tif", ".tiff"]:
            # Create an image plot with the GeoTIFF.
            dataset = rxr.open_rasterio(filename = data_file_path)#.drop_vars(names = ["spatial_ref"], errors = "ignore")
            print(dataset)
            # thing = gv.Dataset(
            #     dataset,
            #     kdims = list(dataset.dims),
            #     vdims = "Elevation (meters)",
            #     nan_nodata = True,
            #     # label = self._selected_collection_info.get(data_file_path, "{}: {}".format(subdir_name, filename))
            # ).to(gv.Image)
            thing = gv.util.from_xarray(
                da = dataset,
                # kims = list(dataset.dims),
                vdims = "Elevation (meters)",
                nan_nodata = True,
                # label = self._selected_collection_info.get(data_file_path, "{}: {}".format(subdir_name, filename))
            )
            print(thing)
            plot = rasterize(
                # gv.load_tiff(
                #     data_file_path,
                #     vdims = "Elevation (meters)",
                #     nan_nodata = True,
                #     # label = self._selected_collection_info.get(data_file_path, "{}: {}".format(subdir_name, filename))
                # )
                thing
            ).opts(
                cmap = "Turbo",
                tools = ["hover"],
                alpha = 0.5,
                responsive = True
            )
        if plot is None:
            print("Error displaying", filename, "as a point/image plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            # Save the created plot.
            self._created_plots[data_file_path] = plot
        end_time = time.time()
        print("Creating data plot for {} took {} seconds.".format(data_file_path, end_time - start_time))
    
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
            plot = self._plot_geojson_linestrings(file_path)
        elif extension in [".parq", ".parquet"]:
            plot = self._plot_parquet_linestrings(file_path)
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
        # with pn.param.set_values(self._data_map_plot, loading = True):
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
                if filename.endswith(".geojson"): transect_file_geodataframe = gpd.read_file(filename = file_path)
                elif filename.endswith(".parq"): transect_file_geodataframe = dask_geopandas.read_parquet(file_path).compute()
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
        # with pn.param.set_values(self._data_map_plot, loading = True):
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
            self._data_file_options_dict = {}
            collection_subdirs = [file for file in os.listdir(self._collection_dir_path) if os.path.isdir(os.path.join(self._collection_dir_path, file)) and (file != self._transects_folder_name)]
            for subdir in collection_subdirs:
                subdir_path = os.path.join(self._collection_dir_path, subdir)
                for file in [file for file in os.listdir(subdir_path) if os.path.isfile(os.path.join(subdir_path, file)) or file.endswith(".parq") or file.endswith(".parquet")]:
                    data_file_path = os.path.join(subdir_path, file)
                    self._data_file_options_dict[file] = data_file_path
            # Get the transect widget's new options.
            transects_dir_path = os.path.join(self._collection_dir_path, self._transects_folder_name)
            if os.path.isdir(transects_dir_path):
                self._all_transect_files = [file for file in os.listdir(transects_dir_path) if os.path.isfile(os.path.join(transects_dir_path, file)) or file.endswith(".parq") or file.endswith(".parquet")]
            else:
                self._all_transect_files = []
            self._transects_multichoice.options = self._all_transect_files + [self._create_own_transect_option]
            self._transects_multichoice.value = []
            # Reassign colors and tap streams for the new collection's transects.
            if self._all_transect_files:
                for i, transect_option in enumerate(self._transects_multichoice.options):
                    self._transect_colors[transect_option] = self._palette_colors[i % self._total_palette_colors]
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

    @param.depends("transects", watch = True)
    def _update_selected_transects_plot(self) -> None:
        """
        Creates an overlay of path plots whenever the selected transect files change.
        """
        # Only when the widget is initialized and at least one transect file is selected...
        if self.transects is not None:
            # Create an overlay of path plots with transects from each selected transect file.
            new_transects_plot = None
            # If the user wants to create their own transect, display widgets related to the user-drawn transect.
            if self._create_own_transect_option in self.transects:
                self._view_user_transect_time_series_button.visible = True
                self._user_drawn_transect_download_button.visible = True
                self._display_user_drawn_transect_instructions = True
            else:
                self._view_user_transect_time_series_button.visible = False
                self._user_drawn_transect_download_button.visible = False
                self._display_user_drawn_transect_instructions = False
            # Update the accordion sections that contain the widgets related to the user-drawn transect.
            self.update_accordion_section = True
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
    
    @param.depends("data_file_paths", watch = True)
    def _update_selected_data_plots(self) -> None:
        """
        Creates an overlay of point or image plots whenever the list of paths for time-series data changes.
        """
        # with pn.param.set_values(self._data_map_plot, loading = True):
        print("_update_selected_data_plots", self.data_file_paths)
        # Only when the list of time-series data files is initiated...
        if self.data_file_paths is not None:
            # Overlay all data files' plots.
            start_time = time.time()
            new_data_plot = None
            for file_path in self.data_file_paths:
                # Create the selected data file's plot if we never read the file before.
                if file_path not in self._created_plots: self._create_data_plot(file_path)
                # Display the data file's plot if it was created.
                # ^ plots aren't created for unsupported files
                if file_path in self._created_plots:
                    if new_data_plot is None:
                        new_data_plot = self._created_plots[file_path]
                    else:
                        new_data_plot = (new_data_plot * self._created_plots[file_path])            
            end_time = time.time()
            print("Overlaying all data plots took {} seconds.".format(end_time - start_time))
            # Save the new data plot.
            self._selected_data_plot = new_data_plot
        else:
            self._selected_data_plot = None

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
        if (self.transects is not None) and (len(self.transects) == 1) and (self._create_own_transect_option in self.transects) and (not len(self._user_transect_plot.data)) and (not self.data_file_paths):
            plot.handles["x_range"].start = plot.handles["x_range"].reset_start = -20037508.342789244
            plot.handles["x_range"].end = plot.handles["x_range"].reset_end = 20037508.342789244
            plot.handles["y_range"].start = plot.handles["y_range"].reset_start = -20037508.342789248
            plot.handles["y_range"].end = plot.handles["y_range"].reset_end = 20037508.342789248

    # -------------------------------------------------- Public Class Properties & Methods --------------------------------------------------
    @param.depends("_update_basemap_plot", "_update_collection_objects", "_update_selected_transects_plot", "_update_selected_data_plots", "_get_clicked_transect_info")
    def plot(self) -> gv.Overlay:
        """
        Returns the selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        start_time = time.time()
        # Overlay the selected plots.
        current_active_tools = ["pan", "wheel_zoom"]
        new_plot = self._selected_basemap_plot
        if self._selected_data_plot is not None:
            new_plot = (new_plot * self._selected_data_plot)
        if self._selected_transects_plot is not None:
            new_plot = (new_plot * self._selected_transects_plot)
            if self._create_own_transect_option in self.transects:
                current_active_tools.append("poly_draw")
        mid_time = time.time()
        print("Plotting all selected plots on the map took {} seconds.".format(mid_time - start_time))
        # Save the overlaid plots.
        self._data_map_plot = new_plot.opts(
            xaxis = None, yaxis = None,
            tools = ["zoom_in", "zoom_out", "tap"],
            active_tools = current_active_tools,
            toolbar = "below",
            title = "", show_legend = True,
            hooks = [self._update_map_data_ranges]
        )
        end_time = time.time()
        print("Rendering new data map on browser took {} seconds.".format(end_time - mid_time))
        # Display browser popup for any errors that occurred while updating the data map.
        if self._error_messages:
            self._error_popup_text.value = "\n".join(self._error_messages)
            # Reset the error messages to an empty list in order to indicate that there are no errors by default.
            self._error_messages = []
        return self._data_map_plot

    def get_sidebar_widgets(self) -> list[any]:
        """
        Returns a list of parameters (will have default widget) or custom Panel widgets for parameters used in the app.
        """
        widgets = [
            self.param.basemap,
            self._collection_select,
            self._transects_multichoice,
            self._error_popup_text,
            self._view_user_transect_time_series_button,
            self._user_drawn_transect_download_button
        ]
        return widgets
    
    def get_accordion_sections(self) -> list:
        """
        Returns a list of tuples, each containing the name of the accordion section and its content.
        """
        sections = []
        if self._display_user_drawn_transect_instructions: sections.append(self._drawing_user_transect_accordion_section)
        return sections

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
    def error_message(self) -> pn.widgets.TextInput:
        """
        Returns the widget containing error messages to display in a browser popup window if there's any errors while creating plots for the data map or popup modal.
        """
        return self._error_popup_text

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
    
    @property
    def all_data_cols(self) -> list[str]:
        """
        Returns a list of column names containing data for the selected collection's time-series.
        """
        return self._all_time_series_data[self.collection]