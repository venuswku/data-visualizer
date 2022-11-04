# Standard library imports
import os

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
from geoviews import opts
import pandas as pd
import rioxarray as rxr
import cartopy
from bokeh.palettes import Bokeh

class DataMap(param.Parameterized):
    # Parameters with GUI widgets
    basemap = param.Selector(label = "Basemap")
    categories = param.ListSelector(label = "Data Categories")

    def __init__(self, data_dir_path: str, latitude_col_names: list[str], longitude_col_names: list[str], map_center: tuple = (0, 0), colors: dict = {}, data_details_button: "ipywidgets.Button" = None, basemap_options: dict = {"Default": gts.OSM}, **params) -> None:
        """
        Creates a new instance of the DataVisualizer class with its instance variables.

        Args:
            data_dir_path (str): Path to the directory containing all the data category subfolders and their data files
            latitude_col_names (list[str]): Possible names of the column containing the latitude of each data point
            longitude_col_names (list[str]): Possible names of the column containing the longitude of each data point
            map_center (tuple): Optional (Latitude, Longitude) tuple specifying the center of the map
            colors (dict): Optional dictionary mapping names of data categories (keys) to their colors (values)
            data_details_button ("ipywidgets.Button"): Optional button displayed in the popup of a hovered/clicked data point
            basemap_options (dict): Optional dictionary mapping basemap names (keys) to basemap WMTS (web mapping tile source) layers (values)
            legend_name (str): Optional name for the map legend, default empty string means that no title will be displayed in the legend
        """
        super().__init__(**params)

        # Set constants.
        self._data_dir_path = data_dir_path
        self._all_lat_cols = latitude_col_names
        self._all_long_cols = longitude_col_names
        # self._default_geojson_hover_color = "#2196f3"

        # Initialize internal class properties.
        # _all_basemaps = dictionary mapping basemap names (keys) to basemap WMTS (web mapping tile source) layers (values)
        self._all_basemaps = basemap_options
        # _selected_basemap_plot = WMTS (web mapping tile source) layer containing the user's selected basemap
        self._selected_basemap_plot = basemap_options[list(basemap_options.keys())[0]]
        # _all_categories = list of data categories (subfolders in the root data directory)
        self._all_categories = [file for file in os.listdir(data_dir_path) if os.path.isdir(data_dir_path + "/" + file)]
        # _category_colors = dictionary mapping data categories (keys) to their color (values) in data plots
        self._category_colors = {}
        # _category_markers = dictionary mapping data categories (keys) to their marker (values) in data plots
        self._category_markers = {}
        # _selected_categories_plot = overlay of data point plots for the user's selected data categories
        self._selected_categories_plot = None
        # _created_plots = dictionary mapping data filenames (keys) to their created data plots (values)
        self._created_plots = {}
        # _displayed_plots = set of unique data filenames that are displayed in the map
        self._displayed_plots = set()

        # Set basemap widget's options.
        self.param.basemap.objects = basemap_options.keys()

        # Set data category widget's options.
        self._categories_multichoice = pn.widgets.MultiChoice.from_param(
            parameter = self.param.categories,
            options = self._all_categories,
            placeholder = "Choose one or more data categories to display",
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
        
        # _plotter = instance of the DataPlotter class, which creates plots with given data
        # self._plotter = DataPlotter(
        #     data_dir_path = data_dir_path,
        #     category_colors = self._category_colors
        # )

    def _create_data_plot(self, filename: str, category: str, plots: gv.Overlay) -> gv.Overlay:
        """
        Creates and displays a plot containing the given file's data points.

        Args:
            filename (str): Name of the file containing data
            category (str): Name of the data category that the file belongs to
            plots (gv.Overlay): Overlaid plots that have been accumulated so far for displaying on the map
        """
        # print("create", filename)
        # Read the file and create a point plot from it.
        file_path = self._data_dir_path + "/" + category + "/" + filename
        [name, extension] = os.path.splitext(filename)
        extension = extension.lower()
        plot = None
        if extension in [".csv", ".txt"]:
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
            plot = data.to(
                gv.Points,
                kdims = [longitude_col, latitude_col],
                vdims = non_lat_long_cols,
                label = category
            ).options(
                opts.Points(
                    color = self._category_colors[category],
                    marker = self._category_markers[category],
                    tools = ["hover"],
                    size = 10
                )
            )
        elif extension == ".asc":
            geotiff_path = self._data_dir_path + "/" + category + "/" + name + ".tif"
            # Convert ASCII grid file into a new GeoTIFF (if not created yet).
            if not os.path.exists(geotiff_path):
                dataset = rxr.open_rasterio(file_path)
                dataset.rio.to_raster(
                    raster_path = geotiff_path,
                    driver = "GTiff"
                )
            # Create an image plot with the GeoTIFF.
            plot = gv.load_tiff(
                geotiff_path,
                crs = cartopy.crs.epsg(26914),
                nan_nodata = True
            )
        
        if plot is None:
            # Return the given overlay of plots if no new plot was created.
            print("Input files with the", extension, "file format is not supported yet.")
            return plots
        else:
            # Save the created plot.
            self._created_plots[filename] = plot
            # Display the created plot.
            new_plots = self._display_data_plot(filename, plots)
            return new_plots

    def _display_data_plot(self, filename: str, plots: gv.Overlay) -> gv.Overlay:
        """
        Displays a data plot on the map.

        Args:
            filename (str): Name of the file containing data
            plots (gv.Overlay): Overlaid plots that have been accumulated so far for displaying on the map
        """
        # print("display", filename)
        data_plot = self._created_plots[filename]
        # Assign the data plot to plots if there's no plots accumulated yet.
        if plots is None:
            new_plots = data_plot
        # Else overlay the data plot with the other accumulated data plots.
        else:
            new_plots = plots * data_plot
        # Add the data filename to the set of displayed plots.
        self._displayed_plots.add(filename)
        return new_plots

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
        new_plot = self._all_basemaps[selected_basemap_name]
        # Save basemap plot.
        self._selected_basemap_plot = new_plot

    @param.depends("categories", watch = True)
    def _update_category_plots(self) -> None:
        """
        Creates a plot with data from the selected categories whenever the selected data categories change.
        """
        # Get the selected data categories.
        if self.categories is None:
            return None
        else:
            # Create a plot with data from each selected category that is within the selected datetime range.
            new_displayed_plots = None
            selected_category_names = self.categories
            for category in self._all_categories:
                category_files = os.listdir(self._data_dir_path + "/" + category)
                for file in category_files:
                    if (category in selected_category_names):# and self._data_within_date_range(file):
                        # Create and display the selected data if we never read the file before.
                        if file not in self._created_plots:
                            new_displayed_plots = self._create_data_plot(file, category, new_displayed_plots)
                        # Display the selected data if it isn't in map yet.
                        else:
                            new_displayed_plots = self._display_data_plot(file, new_displayed_plots)
                    # Else hide the data if user didn't select to display it.
                    else:
                        self._displayed_plots.discard(file)
            # Save overlaid category plots.
            self._selected_categories_plot = new_displayed_plots

    @param.depends("_update_basemap_plot", "_update_category_plots")
    def plot(self) -> gv.Overlay:
        """
        Returns selected basemap and data plots as an overlay whenever any of the plots are updated.
        """
        # Overlay the selected plots.
        if self._selected_categories_plot is None:
            new_plot = self._selected_basemap_plot
        else:
            new_plot = (self._selected_basemap_plot * self._selected_categories_plot)
        # Return the overlaid plots.
        return new_plot.options(
            xaxis = None, yaxis = None,
            active_tools = ["pan", "wheel_zoom"],
            tools = ["zoom_in", "zoom_out", "save"],
            toolbar = "above",
            title = "",
            show_legend = True
        )

    @property
    def param_widgets(self):
        """
        Returns a list of parameters (will have default widget) or custom Panel widgets for parameters used in the app.
        """
        return [
            self.param.basemap,
            self._categories_multichoice
        ]

class DataPlotter(param.Parameterized):
    def __init__(self, **params):
        super().__init__(**params)
        self.original_dataset = gv.DynamicMap()

    # @property
    # def plot(self):
    #    return self.original_dataset

class Application(param.Parameterized):
    # Main components
    data_map = param.ClassSelector(class_ = DataMap, is_instance = True)

    # Parameters with GUI widgets
    

    def __init__(self, **params):
        super().__init__(**params)