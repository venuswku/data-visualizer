# Standard library imports
import os

# External dependencies imports
import param
import panel as pn
import geoviews as gv
import geoviews.tile_sources as gts
from geoviews import opts
import pandas as pd
from bokeh.palettes import Bokeh

class DataMap(param.Parameterized):
    # UI elements
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

        # _map = map containing data that user wants to visualize
        self._map = gv.DynamicMap(self.plot)

        # # popup = popup that displays information about a data point
        # popup_children = (HTML(),)
        # if data_details_button is not None: popup_children += (data_details_button,)
        # self.popup = Popup(
        #     child = VBox(children=(popup_children)),
        #     min_width = 300, max_width = 500,
        #     auto_close = False, name = "Popup"
        # )
        # self._map.add_layer(self.popup)


        # # _selected_geojson_data = dictionary with details (file path, feature with popup info, etc.) about the hovered/clicked GeoJSON feature
        # self._selected_geojson_data = {}
        
        # all_layers = {name1: layer1, name2: layer2, ...} dictionary to store all possible layers that could be on the map
        # ^ e.g. {
        #     "ew15_july_topo.txt": GeoJSON layer 1,
        #     "ew16_july_topo.txt": GeoJSON layer 2,
        #     "Default": default tile layer,
        #     "Satellite": satellite tile layer,
        #     ...
        # }
        # self.all_layers = defaultdict(lambda: None)
        
        # # Add placeholder layers for all data files to the map (initially no features) since new map layers currently can't be added once map is rendered on Panel app.
        # # ^ Will modify GeoJSON layer's `data` attribute when its data needs to be displayed.
        # legend_colors, palette_colors = {}, Bokeh[8]
        # data_categories = [file for file in os.listdir(data_dir_path) if os.path.isdir(data_dir_path + "/" + file)]
        # category_idx, total_palette_colors = 0, len(palette_colors)
        # for category in data_categories:
        #     category_path = data_dir_path + "/" + category
        #     category_files = [file for file in os.listdir(category_path)]
        #     # Assign category to a default color.
        #     default_category_color = palette_colors[category_idx % total_palette_colors]
        #     legend_colors[category] = default_category_color
        #     for file in category_files:
        #         placeholder_geojson = GeoJSON(data = self._geojsons[self._empty_geojson_name], name = file)
        #         # Initially set GeoJSON layer to a color from the Bokeh palette.
        #         placeholder_geojson.point_style = {
        #             "color": default_category_color,
        #             "opacity": 0.5,
        #             "fillColor": default_category_color,
        #             "fillOpacity": 0.3,
        #             "radius": 8,
        #             "weight": 1,
        #             "dashArray": 2
        #         }
        #         placeholder_geojson.hover_style = {
        #             "color": self._default_geojson_hover_color,
        #             "fillColor": self._default_geojson_hover_color,
        #             "weight": 3
        #         }
        #         # Set any custom styles.
        #         geojson_style_attributes = ["style", "point_style", "hover_style"]
        #         if category in category_styles:
        #             category_custom_style = category_styles[category]
        #             for style_attr, style_val in category_custom_style.items():
        #                 if style_attr in geojson_style_attributes:
        #                     if style_attr == "style":
        #                         placeholder_geojson.style = style_val
        #                         if "color" in style_val: legend_colors[category] = style_val["color"]
        #                     elif style_attr == "point_style":
        #                         placeholder_geojson.point_style = style_val
        #                         if "color" in style_val: legend_colors[category] = style_val["color"]
        #                     elif style_attr == "hover_style":
        #                         placeholder_geojson.hover_style = style_val
        #         # Add the placeholder GeoJSON layer to the _map.
        #         self._map.add_layer(placeholder_geojson)
        #     category_idx += 1

        # # Add a map legend if the GeoJSON data layers have different styling.
        # if len(legend_colors) > 1:
        #     self._map.add_control(
        #         LegendControl(
        #             name = legend_name,
        #             legend = legend_colors,
        #             position = "bottomright"
        #         )
        #     )
        
        # plotter = instance of the DataPlotter class, which creates plots with given data
        # self.plotter = DataPlotter(data_dir_path=data_dir_path, category_colors=legend_colors)

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
        dataframe = pd.read_csv(self._data_dir_path + "/" + category + "/" + filename)
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
                    if (category in selected_category_names): # and data_within_date_range(file):
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
    def __init__(self):
        self.original_dataset = gv.DynamicMap()

    # @property
    # def plot(self):
    #    return self.original_dataset

class Application(param.Parameterized):
    # Main components
    data_map = param.ClassSelector(class_ = DataMap, is_instance = True)

    # UI elements
    

    def __init__(self, **params):
        # self.map = pn.pane.HoloViews()
        super().__init__(**params)