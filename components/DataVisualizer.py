# Standard library imports
import os
import json
from collections import defaultdict
import math

# External dependencies imports
import pandas as pd
import geopandas
from ipyleaflet import Map, basemaps, basemap_to_tiles, GeoJSON, Popup, LayersControl, FullScreenControl, LegendControl
from ipywidgets import Layout, HTML, VBox
from bokeh.palettes import Bokeh
from components.DataPlotter import DataPlotter

# Constants
default_geojson_hover_color = "#2196f3"
empty_geojson_name = "no_data"

class DataVisualizer:
  def __init__(self, data_dir_path: str, map_center: tuple = (0, 0), category_styles: dict = {}, data_details_button: "ipywidgets.Button" = None, basemap_options: dict = {"Default": basemaps.OpenStreetMap.Mapnik}, legend_name: str = "") -> None:
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
    # map = map containing data that user wants to visualize
    self.map = Map(
      center = map_center,
      zoom = 15, max_zoom = 18, layout = Layout(height="calc(100vh - 94px)")
    )
  
    self.map.add_control(FullScreenControl())
    self.map.add_control(LayersControl(position="topright"))

    # popup = popup that displays information about a data point
    popup_children = (HTML(),)
    if data_details_button is not None: popup_children += (data_details_button,)
    self.popup = Popup(
      child = VBox(children=(popup_children)),
      min_width = 300, max_width = 500,
      auto_close = False, name = "Popup"
    )
    self.map.add_layer(self.popup)

    # selected_geojson_data = dictionary with details (file path, feature with popup info, etc.) about the hovered/clicked GeoJSON feature
    self.selected_geojson_data = {}

    # geojsons = {name1: GeoJSON1, name2: GeoJSON2, ...} dictionary to store all data that was read from data files
    self.geojsons = {
      empty_geojson_name: {"type": "FeatureCollection", "features": []}
    }
    
    # all_layers = {name1: layer1, name2: layer2, ...} dictionary to store all possible layers that could be on the map
    # ^ e.g. {
    #   "ew15_july_topo.txt": GeoJSON layer 1,
    #   "ew16_july_topo.txt": GeoJSON layer 2,
    #   "Default": default tile layer,
    #   "Satellite": satellite tile layer,
    #   ...
    # }
    self.all_layers = defaultdict(lambda: None)

    # basemaps = list containing all basemap names that the user could choose from
    self.basemaps = basemap_options.keys()
    # Add all basemaps to the map first in order to update the visibility of their tile layers even after the map is rendered.
    for name, basemap in basemap_options.items():
      tile_layer = basemap_to_tiles(basemap)
      tile_layer.show_loading, tile_layer.name = True, name
      tile_layer.visible = False
      self.map.add_layer(tile_layer)
      self.all_layers[name] = tile_layer
    
    # Add placeholder layers for all data files to the map (initially no features) since new map layers currently can't be added once map is rendered on Panel app.
    # ^ Will modify GeoJSON layer's `data` attribute when its data needs to be displayed.
    legend_colors, palette_colors = {}, Bokeh[8]
    data_categories = [file for file in os.listdir(data_dir_path) if os.path.isdir(data_dir_path + "/" + file)]
    category_idx, total_palette_colors = 0, len(palette_colors)
    for category in data_categories:
      category_path = data_dir_path + "/" + category
      category_files = [file for file in os.listdir(category_path)]
      # Assign category to a default color.
      default_category_color = palette_colors[category_idx % total_palette_colors]
      legend_colors[category] = default_category_color
      for file in category_files:
        placeholder_geojson = GeoJSON(data = self.geojsons[empty_geojson_name], name = file)
        # Initially set GeoJSON layer to a color from the Bokeh palette.
        placeholder_geojson.point_style = {"color": default_category_color, "opacity": 0.5, "fillColor": default_category_color, "fillOpacity": 0.3, "radius": 8, "weight": 1, "dashArray": 2}
        placeholder_geojson.hover_style = {"color": default_geojson_hover_color, "fillColor": default_geojson_hover_color, "weight": 3}
        # Set any custom styles.
        geojson_style_attributes = ["style", "point_style", "hover_style"]
        if category in category_styles:
          category_custom_style = category_styles[category]
          for style_attr, style_val in category_custom_style.items():
            if style_attr in geojson_style_attributes:
              if style_attr == "style":
                placeholder_geojson.style = style_val
                if "color" in style_val: legend_colors[category] = style_val["color"]
              elif style_attr == "point_style":
                placeholder_geojson.point_style = style_val
                if "color" in style_val: legend_colors[category] = style_val["color"]
              elif style_attr == "hover_style":
                placeholder_geojson.hover_style = style_val
        # Add the placeholder GeoJSON layer to the map.
        self.map.add_layer(placeholder_geojson)
      category_idx += 1
    
    # Add a map legend if the GeoJSON data layers have different styling.
    if len(legend_colors) > 1:
      self.map.add_control(
        LegendControl(
          name = legend_name,
          legend = legend_colors,
          position = "bottomright"
        )
      )
    
    # plotter = instance of the DataPlotter class, which creates plots with given data
    self.plotter = DataPlotter(data_dir_path=data_dir_path, category_colors=legend_colors)

  def get_existing_property(self, possible_prop_names_and_units: dict, feature_info: dict) -> str:
    """
    Gets the value and optional unit of an existing property from a GeoJSON feature.

    Args:
      possible_prop_names_and_units (dict): Dictionary mapping possible names for the existing property (keys) to optional units corresponding to the property (values)
      feature_info (dict): Dictionary mapping all properties for a GeoJSON feature to their values

    Returns:
      str: Value of the existing property and its optional unit
    """
    for prop_name, prop_unit in possible_prop_names_and_units.items():
      if prop_name in feature_info:
        return str(feature_info[prop_name]) + " " + prop_unit
    return "N/A"

  def get_label_vals(self, values: list[any], feature_info: dict) -> str:
    """
    Gets label values that are displayed in the popup.

    Args:
      values (list[any]): List containing dataframe column names or units that match the popup label
        ^ units or any other text are placed in lists []
        ^ values that have different column names in different files are listed in dictionaries {}, where keys are possible column names and values are units for the value (empty string "" if you don't want to include units)
        ^ column names for the needed data are strings ""
      feature_info (dict): Dictionary mapping all properties for a GeoJSON feature to their values

    Returns:
      str: Values for a label in the popup
    """
    label_values = ""
    for val in values:
      val_type = type(val)
      # If the value is a list, simply add the list value(s) to the returned result because it contains a unit or some literal text.
      if val_type is list:
        label_values += "".join(val)
      # Else if the value is a dictionary, which contains all the possible feature properties and corresponding units for a label, so find and add the value of the existing property.
      elif val_type is dict:
        label_values += self.get_existing_property(val, feature_info)
      # Else the value is a string, which contains the name of an existing dataframe column, so append feature_info[val]. 
      else:
        label_values += str(feature_info[val])
    return label_values

  def get_dataframe_col(self, possible_col_names: list[str], dataframe: "pandas.DataFrame") -> pd.DataFrame:
    """
    Gets the specified column of a dataframe.

    Args:
      possible_col_names (list[str]): List of possible names for the specified column
      dataframe (pandas.DataFrame): Two-dimensional (N columns by N rows) data structure that the specified column should be stored in

    Returns:
      pandas.DataFrame: Specified data column
    """
    for col_name in possible_col_names:
      if col_name in dataframe:
        return dataframe[col_name]
  
  def display_popup_info(self, popup_content: dict, feature: "geojson.Feature", data_file_path: str) -> None:
    """
    Opens the popup at the location of the hovered/clicked GeoJSON feature.

    Args:
      popup_content (dict): Dictionary mapping labels that appear on the popup (keys) to lists containing dataframe column names or units that match the label (values)
        ^ e.g. {
          "Date & Time Collected": ["Date", "Time", [" UTC"]],
          "Orthometric Height": [{"Ortho_Ht_m": "meters", "Ortho_ht_km": "kilometers", "ortho_ht_m": "meters"}]
        }
      feature (geojson.Feature): GeoJSON feature for the data point that had a mouse event
      data_file_path (str): Path to the file containing the hovered/clicked GeoJSON feature
    """
    # Save information about hovered/clicked GeoJSON feature.
    self.selected_geojson_data["path"] = data_file_path
    self.selected_geojson_data["feature"] = feature

    # Create HTML for popup.
    self.popup.location = list(reversed(feature["geometry"]["coordinates"]))
    popup_html = self.popup.child.children[0]
    popup_html.value = ""
    for label, values in popup_content.items():
      popup_html.value += "<b>{}</b> {}<br>".format(label, self.get_label_vals(values, feature["properties"]))
    self.popup.open_popup(location=self.popup.location)
  
  def create_geojson(self, data_path: str, name: str, popup_content: dict, longitude_col_names: list[str], latitude_col_names: list[str]) -> None:
    """
    Creates and displays a GeoJSON layer containing data points (at most 200 points per layer) on the map.

    Args:
      data_path (str): Path to the file that contains the layer's data points
      name (str): Name of the newly added GeoJSON layer
      popup_content (dict): Content displayed in a popup when hovering or clicking on a data point
      longitude_col_names (list[str]): Possible names of the column containing the longitude of each data point
      latitude_col_names (list[str]): Possible names of the column containing the latitude of each data point
    """
    self.popup.close_popup()
    for layer in self.map.layers:
      if layer.name == name:
        # Convert the data file into a GeoJSON.
        dataframe = pd.read_csv(data_path)
        # Get a random sample of 200 data points because large datasets lead to low performance and overcrowded data points.
        max_data_points = 200
        if len(dataframe.index) > max_data_points:
          dataframe = dataframe.sample(max_data_points)
        geodataframe = geopandas.GeoDataFrame(
          dataframe,
          geometry = geopandas.points_from_xy(
            self.get_dataframe_col(longitude_col_names, dataframe),
            self.get_dataframe_col(latitude_col_names, dataframe),
            crs = "EPSG:4326"
          )
        )
        geojson_str = geodataframe.to_json()
        geojson = json.loads(geojson_str)
        # Assign the new GeoJSON data to its corresponding layer in order to display it on the map.
        layer.data = geojson
        # Add mouse event handlers.
        layer.on_click(
          lambda feature, **kwargs: self.display_popup_info(
            popup_content = popup_content,
            feature = feature,
            data_file_path = data_path
          )
        )
        layer.on_hover(
          lambda feature, **kwargs: self.display_popup_info(
            popup_content = popup_content,
            feature = feature,
            data_file_path = data_path
          )
        )
        # Add GeoJSON layer to map and save it.
        self.all_layers[name] = layer
        self.geojsons[name] = geojson
        break
  
  def display_geojson(self, layer_name: str) -> None:
    """
    Displays data for a GeoJSON layer on the map.

    Args:
      layer_name (str): Name of a layer to display data on the map
    """
    self.popup.close_popup()
    if layer_name in self.geojsons:
      layer = self.all_layers[layer_name]
      layer.data = self.geojsons[layer_name]

  def hide_geojson(self, layer_name: str) -> None:
    """
    Hides data for a GeoJSON layer on the map.

    Args:
      layer_name (str): Name of a layer to hide data on the map
    """
    self.popup.close_popup()
    if layer_name in self.all_layers:
      layer = self.all_layers[layer_name]
      layer.data = self.geojsons[empty_geojson_name]

  def update_basemap(self, event: dict) -> None:
    """
    Updates the visibility of all basemap tile layers in order to display the newly selected basemap.

    Args:
      event (dict): information on an event that gets fired when the selected basemap value changes
    """
    newly_selected_basemap = event.new
    for basemap_name in self.basemaps:
      basemap_tile_layer = self.all_layers[basemap_name]
      if basemap_tile_layer.name == newly_selected_basemap: basemap_tile_layer.visible = True
      else: basemap_tile_layer.visible = False