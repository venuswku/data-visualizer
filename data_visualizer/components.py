# Standard library imports
import os

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
from geoviews import opts
import rioxarray as rxr
import pandas as pd
import geopandas as gpd
import cartopy.crs as ccrs
import holoviews as hv
from bokeh.palettes import Bokeh

class DataMap(param.Parameterized):
    # Parameters with GUI widgets
    basemap = param.Selector(label = "Basemap")
    categories = param.ListSelector(label = "Data Categories")
    transects = param.ListSelector(label = "Transects")

    def __init__(self, data_dir_path: str, latitude_col_names: list[str], longitude_col_names: list[str], template: pn.template, colors: dict = {}, basemap_options: dict = {"Default": gts.OSM}, **params) -> None:
        """
        Creates a new instance of the DataMap class with its instance variables.

        Args:
            data_dir_path (str): Path to the directory containing all the data category subfolders and their data files
            latitude_col_names (list[str]): Possible names of the column containing the latitude of each data point
            longitude_col_names (list[str]): Possible names of the column containing the longitude of each data point
            template (panel.template): Data visualizer app's template
            colors (dict): Optional dictionary mapping names of data categories (keys) to their colors (values)
            basemap_options (dict): Optional dictionary mapping basemap names (keys) to basemap WMTS (web mapping tile source) layers (values)
            legend_name (str): Optional name for the map legend, default empty string means that no title will be displayed in the legend
        """
        super().__init__(**params)

        # Set constants.
        self._data_dir_path = data_dir_path
        self._transects_folder_name = "Transects"
        self._all_lat_cols = latitude_col_names
        self._all_long_cols = longitude_col_names
        self._app_template = template
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

        # Initialize internal class properties.
        # _all_basemaps = dictionary mapping basemap names (keys) to basemap WMTS (web mapping tile source) layers (values)
        self._all_basemaps = basemap_options
        # _selected_basemap_plot = WMTS (web mapping tile source) layer containing the user's selected basemap
        self._selected_basemap_plot = basemap_options[list(basemap_options.keys())[0]]
        # _all_categories = list of data categories (subfolders in the root data directory -> excludes transects)
        self._all_categories = [file for file in os.listdir(data_dir_path) if os.path.isdir(data_dir_path + "/" + file) and (file != self._transects_folder_name)]
        # _category_colors = dictionary mapping data categories (keys) to their color (values) in data plots
        self._category_colors = {}
        # _category_markers = dictionary mapping data categories (keys) to their marker (values) in point plots
        self._category_markers = {}
        # _selected_categories_plot = overlay of point plots for each data category selected by the user
        # ^ None if the no data files were provided by the user or the user didn't select any data categories
        self._selected_categories_plot = None
        # _all_transect_files = list of files containing transects to display on the map
        self._all_transect_files = []
        if os.path.isdir(data_dir_path + "/" + self._transects_folder_name):
            self._all_transect_files = [file for file in os.listdir(data_dir_path + "/" + self._transects_folder_name)]
        # _transect_colors = dictionary mapping transect files (keys) to their color (values) in path plots
        self._transect_colors = {}
        # _selected_transects_plot = overlay of path plots for each transect file selected by the user
        # ^ None if the transects file isn't provided by the user or the user didn't select a transects file
        self._selected_transects_plot = None
        # _tapped_data_stream = stream that saves the most recently clicked data element (point, path, etc.) on the map
        self._tapped_data_stream = hv.streams.Selection1D(source = self._selected_transects_plot)
        # _timeseries_plot = timeseries plot for the most recently clicked data element
        self._time_series_plot = gv.DynamicMap(self._create_time_series, streams = [self._tapped_data_stream])
        # _created_plots = dictionary mapping filenames (keys) to their created plots (values)
        self._created_plots = {}

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
            options = self._all_transect_files,
            placeholder = "Choose one or more transect files to display",
            solid = False
        )

        # Set color and marker for each data category.
        palette_colors = Bokeh[8]
        total_palette_colors = len(palette_colors)
        markers = ["o", "^", "s", "d", "x", ">", "*", "v", "+", "<"]
        total_markers = len(markers)
        for i, category in enumerate(self._all_categories):
            # Assign the color that the user chose, if provided.
            if category in colors:
                self._category_colors[category] = colors[category]
            # Else assign a color from the Bokeh palette.
            else:
                self._category_colors[category] = palette_colors[i % total_palette_colors]
            # Assign a marker.
            self._category_markers[category] = markers[i % total_markers]
        # Set color for each transect file.
        for i, transect_file in enumerate(self._all_transect_files):
            self._transect_colors[transect_file] = palette_colors[(len(category) + i) % total_palette_colors]

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
        plot = None
        if extension in [".csv", ".txt"]:
            # Convert the data file into a GeoViews Dataset.
            dataframe = pd.read_csv(file_path)
            non_lat_long_cols, latitude_col, longitude_col = [], None, None
            for col in dataframe.columns:
                if col in self._all_lat_cols: latitude_col = col
                elif col in self._all_long_cols: longitude_col = col
                else: non_lat_long_cols.append(col)
            data = gv.Dataset(
                dataframe,
                kdims = non_lat_long_cols
            )
            # Create a point plot with the GeoViews Dataset.
            plot = data.to(
                gv.Points,
                kdims = [longitude_col, latitude_col],
                vdims = non_lat_long_cols,
                label = category
            ).opts(
                opts.Points(
                    color = self._category_colors[category],
                    marker = self._category_markers[category],
                    tools = ["hover"],
                    size = 10, muted_alpha = 0.01
                )
            )
        elif extension == ".asc":
            geotiff_path = self._data_dir_path + "/" + category + "/" + name + ".tif"
            # Convert ASCII grid file into a new GeoTIFF (if not created yet).
            if not os.path.exists(geotiff_path):
                dataset = rxr.open_rasterio(file_path)
                # Add custom projection based on the Elwha data's metadata.
                dataset.rio.write_crs(self._crs, inplace = True)
                # Save the data as a GeoTIFF.
                dataset.rio.to_raster(
                    raster_path = geotiff_path,
                    driver = "GTiff"
                )
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
            self._created_plots[filename] = plot

    def _create_path_plot(self, filename: str) -> None:
        """
        Creates a path plot containing the given file's paths.

        Args:
            filename (str): Name of the file containing paths
        """
        # Read the file and create a plot from it.
        file_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + filename
        [name, extension] = os.path.splitext(filename)
        extension = extension.lower()
        plot = None
        if extension == ".txt":
            geojson_path = self._data_dir_path + "/" + self._transects_folder_name + "/" + name + ".geojson"
            # Convert the data file into a new GeoJSON (if not created yet).
            if not os.path.exists(geojson_path):
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
                                "properties": {"Transect ID": id},
                                "geometry": {
                                    "type": "LineString",
                                    "coordinates": []
                                }
                            }
                            # Add the transect's start point.
                            transect_feature["geometry"]["coordinates"].append(point)
                            transect_feature["properties"]["Start Point (meters)"] = "({}, {})".format(point[0], point[1])
                        else:
                            # Add the transect's end point.
                            transect_feature["geometry"]["coordinates"].append(point)
                            transect_feature["properties"]["End Point (meters)"] = "({}, {})".format(point[0], point[1])
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
            # Create a path plot with the GeoJSON.
            geodataframe = gpd.read_file(geojson_path)
            plot = gv.Path(
                geodataframe,
                crs = self._epsg,
                label = "{}: {}".format(self._transects_folder_name, filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
            ).opts(
                opts.Path(
                    color = self._transect_colors[filename],
                    tools = ["hover", "tap"],
                    active_tools = ["tap"]
                )
            )
        elif extension == ".geojson":
            geodataframe = gpd.read_file(file_path)
            plot = gv.Path(
                geodataframe,
                crs = self._epsg,
                label = "{}: {}".format(self._transects_folder_name, filename)    # HoloViews 2.0: Paths will be in legend by default when a label is specified (https://github.com/holoviz/holoviews/issues/2601)
            ).opts(
                opts.Path(
                    color = self._transect_colors[filename],
                    tools = ["hover", "tap"],
                    active_tools = ["tap"]
                )
            )
        if plot is None:
            print("Error displaying", name + extension, "as a path plot:", "Input files with the", extension, "file format are not supported yet.")
        else:
            # Save the created plot.
            self._created_plots[filename] = plot
    
    def _create_time_series(self, index: int) -> any:
        """
        Creates a time-series plot for the selected data element on the map.

        Args:
            index (int): Index of the selected/clicked/tapped data element
        """
        plot = gv.Points([])
        # If a data element is selected...
        if len(index):
            # Open the app's modal to display the time-series plot.
            self._app_template.open_modal()
            print(index)
            # Get user's selected transect.
            [transect_index] = index
            for displayed_transect_file in self.transects:
                transect_plot = self._created_plots[displayed_transect_file]
                clicked_transect = transect_plot.iloc[transect_index]
                print(clicked_transect)
            # Create time-series plot for all data collected along the selected transect.
            plot = gv.Points([])
        return plot

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
        if (self.categories is not None) and (len(self.categories) > 0):
            # Create a plot with data from each selected category that is within the selected datetime range.
            new_data_plot = None
            selected_category_names = self.categories
            for category in self._all_categories:
                category_files = os.listdir(self._data_dir_path + "/" + category)
                for file in category_files:
                    if (category in selected_category_names):# and self._data_within_date_range(file):
                        # Create the selected data's point plot if we never read the file before.
                        if file not in self._created_plots:
                            self._create_data_plot(file, category)
                        # Display the data file's point plot if it was created.
                        # ^ plots aren't created for unsupported files -> e.g. png files don't have data points
                        if file in self._created_plots:
                            if new_data_plot is None:
                                new_data_plot = self._created_plots[file]
                            else:
                                new_data_plot = (new_data_plot * self._created_plots[file])
            # Save overlaid category plots.
            self._selected_categories_plot = new_data_plot

    @param.depends("transects", watch = True)
    def _update_selected_transects_plot(self) -> None:
        """
        Creates an overlay of path plots whenever the selected transect files change.
        """
        # Only when the widget is initialized and at least one transect file is selected...
        if (self.transects is not None) and (len(self.transects) > 0):
            # Create an overlay of path plots with transects from each selected transect file.
            new_transects_plot = None
            for file in self.transects:
                # Create the selected transect file's path plot if we never read the file before.
                if file not in self._created_plots:
                    self._create_path_plot(file)
                # Display the transect file's path plot if it was created.
                # ^ plots aren't created for unsupported files
                if file in self._created_plots:
                    if new_transects_plot is None:
                        new_transects_plot = self._created_plots[file]
                    else:
                        new_transects_plot = (new_transects_plot * self._created_plots[file])
            # Set tap stream to new transects overlay.
            self._tapped_data_stream.source = new_transects_plot
            # Save overlaid transect plots.
            self._selected_transects_plot = new_transects_plot

    @param.depends("_update_basemap_plot", "_update_selected_categories_plot", "_update_selected_transects_plot")
    def plot(self) -> gv.Overlay:
        """
        Returns selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        # Overlay the selected plots.
        new_plot = self._selected_basemap_plot
        if self._selected_categories_plot is not None:
            new_plot = (new_plot * self._selected_categories_plot)
        if self._selected_transects_plot is not None:
            new_plot = (new_plot * self._selected_transects_plot)
        # Return the overlaid plots.
        return new_plot.opts(
            xaxis = None, yaxis = None,
            active_tools = ["pan", "wheel_zoom"],
            toolbar = None, title = "", show_legend = True
        )

    @property
    def param_widgets(self) -> list:
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
    def time_series_plot(self) -> gv.DynamicMap:
        """
        Returns a time-series plot for the selected data element on the map.
        """
        return self._time_series_plot

class Application(param.Parameterized):
    # Main components
    data_map = param.ClassSelector(class_ = DataMap, is_instance = True)

    # Parameters with GUI widgets
    # view_time_series = param.Event(label = "View Time-Series Along Transect")

    def __init__(self, **params):
        """
        Creates a new instance of the Application class with its instance variables.
        """
        super().__init__(**params)
    
    # @param.depends("view_time_series", watch = True)
    # def _view_time_series(self):
    #     # Open the app's modal to display the time-series plot.
    #     self._app_template.open_modal()