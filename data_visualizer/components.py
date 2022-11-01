# Standard library imports
import os
import json
from collections import defaultdict

# External dependencies imports
import panel as pn
import param
import geoviews as gv
import geoviews.tile_sources as gts
from ipyleaflet import Map, basemap_to_tiles, GeoJSON, Popup, LayersControl, FullScreenControl, LegendControl
from ipywidgets import Layout, HTML, VBox
from bokeh.palettes import Bokeh

class DataMap(param.Parameterized):
    # UI elements
    basemap = param.Selector(label="Basemap")

    def __init__(self, data_dir_path: str, map_center: tuple = (0, 0), category_styles: dict = {}, data_details_button: "ipywidgets.Button" = None, basemap_options: dict = {"Default": gts.OSM}, legend_name: str = "", **params) -> None:
        """
        Creates a new instance of the DataVisualizer class with its instance variables.

        Args:
        data_dir_path (str): Path to the directory containing all the data category subfolders and their data files
        map_center (tuple): Optional (Latitude, Longitude) tuple specifying the center of the map
        category_styles (dict): Optional dictionary mapping names of data categories (keys) to their optional styling on a GeoJSON layer (values)
            ^ e.g. {
                "category1": {
                    "point_style": {"opacity": 0.5, "radius": 1},
                    "hover_style": {"color": "blue", "radius": 2}
                },
                "category2": {
                    "style": {"color": "black", "opacity": 0.5, "radius": 1},
                    "hover_style": {"color": "red", "opacity": 1, "radius": 2}
                }
            }
            ^ Make sure keys in data[data_category] match a GeoJSON attribute (https://ipyleaflet.readthedocs.io/en/latest/layers/geo_json.html#attributes-and-methods) and values are dictionaries {}.
        data_details_button ("ipywidgets.Button"): Optional button displayed in the popup of a hovered/clicked data point
        basemap_options (dict): Optional dictionary mapping basemap names (keys) to basemap layers (values)
        legend_name (str): Optional name for the map legend, default empty string means that no title will be displayed in the legend
        """
        super().__init__(**params)

        self._default_geojson_hover_color = "#2196f3"
        self._empty_geojson_name = "no_data"

        # Set basemap selector's options.
        self.param.basemap.objects = basemap_options.keys()
        self._all_basemaps = basemap_options

        # _map = map containing data that user wants to visualize
        # self._map = Map(
        #     center = map_center,
        #     zoom = 15, max_zoom = 18,
        #     layout = Layout(height="calc(100vh - 94px)")
        # )
        self._map = gv.DynamicMap(self.plot)
    
        # self._map.add_control(FullScreenControl())
        # self._map.add_control(LayersControl(position="topright"))

        # # popup = popup that displays information about a data point
        # popup_children = (HTML(),)
        # if data_details_button is not None: popup_children += (data_details_button,)
        # self.popup = Popup(
        #     child = VBox(children=(popup_children)),
        #     min_width = 300, max_width = 500,
        #     auto_close = False, name = "Popup"
        # )
        # self._map.add_layer(self.popup)

        # _selected_geojson_data = dictionary with details (file path, feature with popup info, etc.) about the hovered/clicked GeoJSON feature
        self._selected_geojson_data = {}

        # _geojsons = {name1: GeoJSON1, name2: GeoJSON2, ...} dictionary to store all data that was read from data files
        self._geojsons = {
            self._empty_geojson_name: {"type": "FeatureCollection", "features": []}
        }
        
        # all_layers = {name1: layer1, name2: layer2, ...} dictionary to store all possible layers that could be on the map
        # ^ e.g. {
        #     "ew15_july_topo.txt": GeoJSON layer 1,
        #     "ew16_july_topo.txt": GeoJSON layer 2,
        #     "Default": default tile layer,
        #     "Satellite": satellite tile layer,
        #     ...
        # }
        self.all_layers = defaultdict(lambda: None)
        
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

    # def update_basemap(self) -> None:
    #     """
    #     Updates the visibility of all basemap tile layers in order to display the newly selected basemap.
    #     """
    #     print("update", self.param.basemap)
    #     # for basemap_name in self.basemap.objects:
    #     #     basemap_tile_layer = self.all_layers[basemap_name]
    #     #     if basemap_tile_layer.name == self.basemap: basemap_tile_layer.visible = True
    #     #     else: basemap_tile_layer.visible = False

    @param.depends("basemap")
    def plot(self):
        print("plot", self.basemap)
        if self.basemap is None:
            selected_basemap_name = list(self._all_basemaps.keys())[0]
        else:
            selected_basemap_name = self.basemap
        return self._all_basemaps[selected_basemap_name]

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