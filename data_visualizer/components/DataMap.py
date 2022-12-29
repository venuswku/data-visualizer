# Standard library importspn
import os

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import holoviews as hv
from bokeh.palettes import Bokeh

### DataMap is used for displaying the inputted data files onto a map. ###
class DataMap(param.Parameterized):
    # -------------------------------------------------- Parameters --------------------------------------------------
    basemap = param.Selector(label = "Basemap")
    categories = param.ListSelector(label = "Data Categories")
    transects = param.ListSelector(label = "Transects")
    clicked_transects_info = param.Dict(default = {})

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
        self._create_own_transect_option = "Create My Own Transect"
        # _geodata_folder_name = Name of the folder containing GeoJSON/GeoTIFF files that were created by georeferencing data files (txt, csv, asc)
        # ^ allows data to load faster onto the map
        self._geodata_folder_name = "GeoData"
        # _clicked_transects_file_key = clicked_transects_info parameter's dictionary key that corresponds to the file that contains the clicked transect(s)
        self._clicked_transects_file_key = "transect_file"
        # _num_clicked_transects_key = clicked_transects_info parameter's dictionary key that corresponds to the number of clicked transect(s)
        self._num_clicked_transects_key = "num_clicked_transects"
        # _clicked_transects_longitude_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the longitude/easting of the clicked transects' start and end points
        self._clicked_transects_longitude_key = "longitude_col_name"
        # _clicked_transects_latitude_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the latitude/northing of the clicked transects' start and end points
        self._clicked_transects_latitude_key = "latitude_col_name"
        # _clicked_transects_data_cols_key = clicked_transects_info parameter's dictionary key that corresponds to the names of the columns to display in the popup modal's clicked transect(s) data table
        self._clicked_transects_data_cols_key = "table_cols"
        # _clicked_transects_id_key = clicked_transects_info parameter's dictionary key that corresponds to the name of the column containing the ID of the clicked transects' start and end points
        self._clicked_transects_id_key = "id_col_name"
        # _transects_id_col_name = Name of the column containing the ID of each clicked transect
        # ^ initially created in _create_transects_geojson() when assigning an ID property for each transect in the outputted GeoJSON
        self._transects_id_col_name = "Transect ID"
        
        # _crs = custom coordinate reference system for the projected data
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
        # _epsg = publicly registered coordinate system for the projected data
        # ^ should be close, if not equivalent, to the custom CRS defined above (_crs)
        self._epsg = ccrs.epsg(32148)

        # -------------------------------------------------- Internal Class Properties --------------------------------------------------
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
        self._categories_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.categories,
            options = self._all_categories,
            placeholder = "Choose one or more data categories to display",
            solid = False
        )
        # _selected_categories_plot = overlay of point plots for each data category selected by the user
        # ^ empty overlay plot if the no data files were provided by the user or the user didn't select any data categories
        self._selected_categories_plot = hv.DynamicMap(
            self._update_selected_categories_plot,
            streams = dict(selected_data_categories = self._categories_multichoice.param.value)
        ).opts(hooks = [self.plot], responsive = True)#
        
        # _all_transect_files = list of files containing transects to display on the map
        self._all_transect_files = []
        transects_dir_path = data_dir_path + "/" + self._transects_folder_name
        if os.path.isdir(transects_dir_path):
            self._all_transect_files = [file for file in os.listdir(transects_dir_path) if os.path.isfile(os.path.join(transects_dir_path, file))]
        # _transect_colors = dictionary mapping each transect file (keys) to a color (values), which will be used for the color of its path plots
        self._transect_colors = {}
        self._transects_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.transects,
            options = self._all_transect_files + [self._create_own_transect_option],
            placeholder = "Choose one or more transect files to display",
            solid = False
        )
        # _selected_transects_plot = overlay of path plots if the user selected one or more transect files to display on the map
        # ^ empty overlay plot if the user didn't provide any transect files or the user didn't select to display a transect file
        self._selected_transects_plot = hv.DynamicMap(
            self._update_selected_transects_plot,
            streams = dict(selected_transect_files = self._transects_multichoice.param.value)
        ).opts(hooks = [self.plot], responsive = True)#

        # _tapped_data_streams = dictionary mapping each transect file's path (keys) to a selection stream (values), which saves the file's most recently clicked data element (path) on the map
        self._tapped_data_streams = {}
        for file in self._all_transect_files:
            file_path = transects_dir_path + "/" + file
            self._tapped_data_streams[file_path] = hv.streams.Selection1D(source = None, rename = {"index": file_path})
            # Specify a callable subscriber function that gets called whenever any transect from the file is clicked/tapped.
            self._tapped_data_streams[file_path].add_subscriber(self._get_clicked_transects_info)

        # _user_transect_plot = predefined path plot if the user wanted to create their own transect to display on the map
        self._user_transect_plot = gv.Path(
            data = [[(296856.9100, 131388.7700), (296416.5400, 132035.8500)]],#[[]],
            crs = self._epsg
        ).opts(active_tools = ["poly_draw"])
        # self._user_transect_plot = hv.Curve(data = np.array([[(-123.5688, 48.1523), (-123.5626, 48.1476)]])).opts(active_tools = ["point_draw"])
        # self._user_transect_plot = hv.Points(
        #     data = np.array([[-123.5688, 48.1523], [-123.5626, 48.1476]])
        # ).opts(active_tools = ["point_draw"], color = "black")
        # self._user_transect_plot = hv.Curve([]).opts(
        #     active_tools = ["point_draw"],
        #     color = "black"
        # )
        # _edit_user_transect_stream = stream that allows user to move the start and end points of their own transect
        self._edit_user_transect_stream = hv.streams.PolyDraw(
            source = self._user_transect_plot,
            num_objects = 1,
            drag = True,
            styles = {"line_color": ["black"], "line_width": [5]},
            show_vertices = True,
            vertex_style = {"fill_color": "black"}
        )
        # self._edit_user_transect_stream = hv.streams.PointDraw(source = self._user_transect_plot, num_objects = 2)
        # self._edit_user_transect_stream = hv.streams.CurveEdit(
        #     # data = self._user_transect_plot.columns(),
        #     source = self._user_transect_plot,
        #     num_objects = 2,
        #     add = False,
        #     style = {"color": "black", "size": 10}
        # )

        # -------------------------------------------------- Widget and Plot Options --------------------------------------------------
        # Set basemap widget's options.
        self.param.basemap.objects = basemap_options.keys()

        # Set data category widget's options.
        # self._categories_multichoice = pn.widgets.MultiChoice.from_param(
        #     parameter = self.param.categories,
        #     options = self._all_categories,
        #     placeholder = "Choose one or more data categories to display",
        #     solid = False
        # )
        print("categories", self._categories_multichoice.param.value)

        # Set transect widget's options.
        # self._transects_multichoice = pn.widgets.MultiChoice.from_param(
        #     parameter = self.param.transects,
        #     options = self._all_transect_files + [self._create_own_transect_option],
        #     placeholder = "Choose one or more transect files to display",
        #     solid = False
        # )
        print("transects", self._transects_multichoice.param.value)

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
        geodataframe = gpd.read_file(geojson_file_path)
        path_plot = gv.Path(
            data = geodataframe,
            crs = self._epsg,
            label = "{}: {}".format(self._transects_folder_name, filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
        ).opts(
            color = self._transect_colors[filename],
            hover_color = self._app_main_color,
            selection_color = self._app_main_color,
            nonselection_color = self._transect_colors[filename],
            nonselection_alpha = 1, selection_alpha = 1,
            tools = ["hover", "tap"], selected = []
        )
        return path_plot

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
            self._create_transects_geojson(file_path, geojson_path)
            # Create a path plot with the path GeoJSON.
            plot = self._plot_geojson_linestrings(geojson_path, filename)
        # Save the path plot, if created.
        if plot is None:
            print("Error displaying", name + extension, "as a path plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            self._created_plots[file_path] = plot

    def _create_transects_geojson(self, file_path: str, geojson_path: str) -> None:
        """
        Creates and saves a GeoJSON file containing LineStrings for each transect in the given transect file.

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
                        transect_feature["properties"]["Start Point"] = "({}, {})".format(x, y)
                    else:
                        # Add the transect's end point.
                        transect_feature["geometry"]["coordinates"].append(point)
                        transect_feature["properties"]["End Point"] = "({}, {})".format(x, y)
                        # Save the transect to the FeatureCollection.
                        features_list.append(transect_feature)
                        # Reset the feature for the next transect.
                        transect_feature = None
            # Convert the FeatureCollection into a GeoJSON.
            geodataframe = gpd.GeoDataFrame.from_features(
                {"type": "FeatureCollection", "features": features_list},
                crs = self._crs
            )
            # Save the GeoJSON file to skip converting the data file again.
            geodataframe.to_file(geojson_path, driver = "GeoJSON")

    def _get_clicked_transects_info(self, **params: dict) -> None:
        """
        Gets information about the clicked transect on the map, which is used to update the popup modal's contents (time-series plot and transect data table).

        Args:
            params (dict): Dictionary mapping each transect file's path (keys) to a list containing the indices of selected/clicked/tapped transects (values) from its transect file
        """
        # print("Selection1D stream's parameter:", params)
        # Set custom names for the latitude and longitude data columns, or set them to None to keep the default column names ("Longitude", "Latitude").
        custom_long_col_name = "Easting (meters)"
        custom_lat_col_name = "Northing (meters)"
        # Save information about the clicked transect(s) in a dictionary.
        clicked_transects_info_dict = {}
        # Find the user's clicked/selected transect(s).
        for file_path in params:
            _, filename = os.path.split(file_path)
            if (filename in self._all_transect_files) and params[file_path]:
                clicked_transect_indices = params[file_path]
                num_clicked_transects = len(clicked_transect_indices)
                # Reset transect file's Selection1D stream parameter to its default value (empty list []).
                self._tapped_data_streams[file_path].reset()
                # self._tapped_data_streams[file_path].event(index = [])
                # Add information about the clicked transect(s).
                clicked_transects_info_dict[self._clicked_transects_file_key] = filename
                clicked_transects_info_dict[self._num_clicked_transects_key] = num_clicked_transects
                clicked_transects_info_dict[self._clicked_transects_id_key] = self._transects_id_col_name
                # Specify the names of columns to display in the popup modal's data table.
                clicked_transects_info_dict[self._clicked_transects_data_cols_key] = []
                # Get all the transects/paths from the clicked transect's file.
                transects_file_plot = self._created_plots[file_path]
                # transects_file_plot.opts(selected = [])
                transect_file_paths = transects_file_plot.split()
                # Get data for each of the user's clicked transect(s).
                for transect_index in clicked_transect_indices:
                    clicked_transects_info = transect_file_paths[transect_index].columns(dimensions = ["Longitude", "Latitude", self._transects_id_col_name])
                    # Rename GeoViews' default coordinate column names if custom column names were provided.
                    for col, values in clicked_transects_info.items():
                        if custom_long_col_name and (col == "Longitude"):
                            col = custom_long_col_name
                            clicked_transects_info_dict[self._clicked_transects_longitude_key] = custom_long_col_name
                        elif custom_lat_col_name and (col == "Latitude"):
                            col = custom_lat_col_name
                            clicked_transects_info_dict[self._clicked_transects_latitude_key] = custom_lat_col_name
                        # Save the column's name if it isn't in the list of data columns to display in the popup modal's data table.
                        if col not in clicked_transects_info_dict[self._clicked_transects_data_cols_key]:
                            clicked_transects_info_dict[self._clicked_transects_data_cols_key].append(col)
                        # Convert the numpy array of column values into a Python list and save to make it easier to iterate over the values.
                        curr_transect_col_vals = values.tolist()
                        # Save column values for the clicked transect.
                        prev_transects_col_vals = clicked_transects_info_dict.get(col, [])
                        clicked_transects_info_dict[col] = prev_transects_col_vals + curr_transect_col_vals
                # Stop iterating through all the transect files once a clicked transect is found.
                break
        # Update the clicked_transects_info parameter in order to update the time-series plot, transect data table, or error message in the popup modal.
        self.clicked_transects_info = clicked_transects_info_dict
        # # Reset the clicked_transects_info parameter in case the user clicks on the same transect again.
        # # ^ If the parameter isn't reset, then the parameter value stays the same, meaning the info won't be sent to the modal and the modal won't open.
        # self.clicked_transects_info = {}

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

    # @param.depends("categories")
    def _update_selected_categories_plot(self, selected_data_categories: list[str]) -> hv.Overlay:#, **params: dict
        """
        Creates an overlay of data plots whenever the selected data categories change.
        """
        print("_update_selected_categories_plot", self.categories, selected_data_categories)
        # print("_update_selected_categories_plot", params)
        # Only when the widget is initialized and at least one data category is selected...
        if self.categories is not None:
            # Create a plot with data from each selected category.
            new_data_plot = None
            for category in self._all_categories:
                category_dir_path = self._data_dir_path + "/" + category
                category_files = [file for file in os.listdir(category_dir_path) if os.path.isfile(os.path.join(category_dir_path, file))]
                for file in category_files:
                    if category in self.categories:
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
            print(new_data_plot)
            # Return overlaid category plots.
            if new_data_plot is not None: return new_data_plot
        # Return an empty/placeholder overlay plot if no data categories were selected.
        return gv.Points(data = []) * gv.Points(data = [])

    @param.depends("_get_clicked_transects_info")#"transects", 
    def _update_selected_transects_plot(self, selected_transect_files: list[str]) -> hv.Overlay:#, **params: dict
        """
        Creates an overlay of path plots whenever the selected transect files change or clicked transects need to be unselected after getting all their info.
        """
        print("_update_selected_transects_plot", self.transects, selected_transect_files)
        # print("_update_selected_transects_plot", params)
        # Only when the widget is initialized and at least one transect file is selected...
        if selected_transect_files:
            # Create an overlay of path plots with transects from each selected transect file.
            new_transects_plot = None
            for file in selected_transect_files:
                file_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + file
                # Allow user to draw start and end points when they selected to draw their own transect.
                if file == self._create_own_transect_option:
                    # Display an editable curve plot for the user to modify their transect's start and end points.
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
            # Return overlaid transect plots.
            print(new_transects_plot)
            if len(selected_transect_files) > 1: return new_transects_plot
            else: return new_transects_plot * gv.Path(data = [])
        else:
            # Return an empty/placeholder overlay plot if no transects were selected.
            return gv.Path(data = []) * gv.Path(data = [])

    # -------------------------------------------------- Public Class Methods --------------------------------------------------
    @param.depends("_update_basemap_plot")#, "_update_selected_categories_plot", "_update_selected_transects_plot"
    def plot(self) -> gv.Overlay:
        """
        Returns selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        default_active_tools = ["pan", "wheel_zoom"]
        # Overlay the selected plots.
        # new_plot = self._selected_basemap_plot# * self._user_transect_plot
        print(self._selected_categories_plot.info, self._selected_transects_plot.info)
        # if self._selected_categories_plot.dimensions():
        #     new_plot = new_plot * self._selected_categories_plot
        # if self._selected_transects_plot.dimensions():
        #     new_plot = new_plot * self._selected_transects_plot

        # if self._create_own_transect_option in self.transects:
        #     default_active_tools.append("poly_draw")
        # Return the overlaid plots.
        # return new_plot.opts(
        #     xaxis = None, yaxis = None,
        #     tools = ["zoom_in", "zoom_out", "poly_draw", "save"],
        #     active_tools = default_active_tools,
        #     toolbar = "above",
        #     # toolbar = None,
        #     title = "", show_legend = True
        # )
        return pn.panel(
            (self._selected_basemap_plot * self._selected_categories_plot * self._selected_transects_plot).opts(
                xaxis = None, yaxis = None,
                tools = ["zoom_in", "zoom_out", "poly_draw", "save"],
                active_tools = default_active_tools,
                toolbar = "above",
                # toolbar = None,
                title = "", show_legend = True
            ),
            sizing_mode = "scale_both"
        )

    @property
    def param_widgets(self) -> list[any]:
        """
        Returns a list of parameters (will have default widget) or custom Panel widgets for parameters used in the app.
        """
        widgets = [
            self.param.basemap,
            self._categories_multichoice
        ]
        # If the user provided file(s) containing transects in a subfolder along with the data categories, then display transect widget.
        if len(self._all_transect_files): widgets.append(self._transects_multichoice)
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
    def epsg(self) -> ccrs.epsg:
        """
        Returns the data's projection for a projected coordinate system corresponding to the stated EPSG code.
        """
        return self._epsg

    @property
    def selected_basemap(self) -> any:
        """
        Returns the user's selected basemap.
        """
        return self._selected_basemap_plot

    @property
    def clicked_transects_info_keys(self) -> list[str]:
        """
        Returns a list of keys that appear in the clicked_transects_info parameter's dictionary,
        which is used by the popup modal to display information about the user's clicked transect(s) from this map.
        """
        return [
            self._clicked_transects_file_key,
            self._num_clicked_transects_key,
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
                    crs = "EPSG:4326"
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
            dataset.rio.write_crs(self._crs, inplace = True)
            # Save the data as a GeoTIFF.
            dataset.rio.to_raster(
                raster_path = geotiff_path,
                driver = "GTiff"
            )